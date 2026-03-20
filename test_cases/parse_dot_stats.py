#!/usr/bin/env python3
"""Parse Hedgehog dot file execution stats and produce a per-node timing report.

Extracts timing info from state manager node labels in the dot file,
categorizes nodes, and computes overhead metrics.

All visible timed nodes in the dot file are state managers (states wrapped
in hh::StateManager). Multi-threaded kernel tasks are inside sub-graphs
and only show queue stats, not execution timing. Kernel time is estimated
as: Graph Time - Sum(all state D+E times).
"""

import re
import sys
import os
from collections import defaultdict
from pathlib import Path
import glob as globmod


def parse_dot_file(path):
    """Parse a Hedgehog dot file and extract node timing stats."""
    with open(path) as f:
        content = f.read()

    # Extract graph-level stats
    graph_match = re.search(r'Execution duration:([\d.]+)([mu]?s)', content)
    graph_time = 0.0
    if graph_match:
        val = float(graph_match.group(1))
        unit = graph_match.group(2)
        if unit == 'ms': val /= 1000.0
        elif unit == 'us': val /= 1e6
        graph_time = val

    nodes = []

    # All timed nodes are state managers with this format (all on one line in the dot file):
    # x0x... [label="Name\nElements: N\nWait: Xs\nLock state: ...\nEmpty ready list: ...\nDequeue+Exec: Xms\nExec/Element: ...\n<extra>",shape=...];
    # Process line-by-line to avoid cross-line false matches.
    nodes = []
    sm_pattern = re.compile(
        r'\[label="([^"\\]+)\\n'         # node name (no quotes or backslashes)
        r'Elements:\s*(\d+)\\n'          # element count
        r'Wait:\s*([\d.]+)([mu]?s)\\n'   # wait time
        r'Lock state:\s*([\d.]+)([mu]?s)\\n'
        r'Empty ready list:\s*([\d.]+)([mu]?s)\\n'
        r'Dequeue\+Exec:\s*([\d.]+)([mu]?s)\\n'  # dequeue+exec time
        r'Exec/Element:\s*([\d.]+)([nmu]?s)'      # exec per element
        r'(?:\\n(.+?))?'                             # optional extra info
        r'"'
    )

    def to_ms(val, unit):
        if unit == 's': return val * 1000.0
        if unit == 'ms': return val
        if unit == 'us': return val / 1000.0
        if unit == 'ns': return val / 1e6
        return val

    for line in content.splitlines():
        m = sm_pattern.search(line)
        if not m:
            continue
        name = m.group(1)
        elements = int(m.group(2))
        wait_ms = to_ms(float(m.group(3)), m.group(4))
        lock_ms = to_ms(float(m.group(5)), m.group(6))
        empty_ms = to_ms(float(m.group(7)), m.group(8))
        dequeue_exec_ms = to_ms(float(m.group(9)), m.group(10))
        exec_per_elem_ms = to_ms(float(m.group(11)), m.group(12))
        extra = m.group(13) or ""
        nodes.append({
            'name': name,
            'elements': elements, 'wait_ms': wait_ms,
            'lock_ms': lock_ms, 'empty_ms': empty_ms,
            'dequeue_exec_ms': dequeue_exec_ms,
            'exec_per_elem_ms': exec_per_elem_ms,
            'exec_ms': exec_per_elem_ms * elements,
            'extra': extra.replace('\\n', '\n'),
        })

    return graph_time, nodes


def categorize_nodes(nodes):
    """Categorize nodes by role in the pipeline."""
    categories = {
        'barrier_states': [],   # Merged BarrierState nodes (collect N, barrier fn, scatter N)
        'orchestrators': [],    # Orchestrator/Decompose/Fork states (dispatch to parallel kernels)
        'collectors': [],       # Collector/Reassemble/Join/Coll states (gather from parallel kernels)
        'loop_control': [],     # TimestepLoop, RetryLoop, TerminationSink, RetryMomDivCollector
        'other': [],
    }

    barrier_state_names = {
        'MeshExchange(1)', 'MeshExchange(3)', 'MeshExchange(4)',
        'MeshExchange(6a)', 'MeshExchange(6b)', 'MeshExchange(7)',
        'Hvac+InitDiv', 'Soot+Hvac', 'RemoveMove',
        'PredDivExchange', 'CorrDivExchange',
        'PredPressure', 'CorrPressure',
        'TimestepDump',
    }

    loop_names = {
        'TimestepLoop', 'TerminationSink', 'RetryLoop', 'RetryMomDivCollector',
        'PressurePostLoop',
    }

    for n in nodes:
        name = n['name']
        if name in barrier_state_names:
            categories['barrier_states'].append(n)
        elif name in loop_names:
            categories['loop_control'].append(n)
        elif ('Collector' in name or 'Reassemble' in name or
              'Coll' in name or 'Mid' in name or
              name.startswith('Join') or name.startswith('PredJoin')):
            categories['collectors'].append(n)
        elif ('Orch' in name or 'Decompose' in name or
              name.startswith('Fork') or name.startswith('PredFork') or
              name == 'PredStep1Orch'):
            categories['orchestrators'].append(n)
        else:
            categories['other'].append(n)

    return categories


def print_report(test_name, graph_time, nodes):
    """Print a formatted timing report."""
    cats = categorize_nodes(nodes)
    graph_ms = graph_time * 1000.0 if graph_time > 0 else 0.001

    print(f"\n{'='*90}")
    print(f"  {test_name}")
    print(f"  Graph execution time: {graph_time:.3f}s ({graph_ms:.1f}ms)")
    print(f"  Total timed nodes: {len(nodes)}")
    print(f"{'='*90}")

    # Barrier States (the merged nodes)
    if cats['barrier_states']:
        print(f"\n  BARRIER STATES ({len(cats['barrier_states'])} nodes)")
        print(f"  {'Name':<25} {'D+E(ms)':>10} {'Elements':>8} {'Avg(ms)':>10} {'Routines'}")
        print(f"  {'-'*25} {'-'*10} {'-'*8} {'-'*10} {'-'*40}")
        total_de = 0
        for n in sorted(cats['barrier_states'], key=lambda x: -x['dequeue_exec_ms']):
            lines = [l for l in n['extra'].split('\n') if l.strip()]
            # Filter: routine names are uppercase, timing lines start with digits
            routines = [l for l in lines if l and not l[0].isdigit()]
            routines_str = ', '.join(routines)
            total_de += n['dequeue_exec_ms']
            print(f"  {n['name']:<25} {n['dequeue_exec_ms']:>10.3f} {n['elements']:>8} {n['exec_per_elem_ms']:>10.4f}   {routines_str[:55]}")
        print(f"  {'TOTAL':<25} {total_de:>10.3f} {'':>8} {'':>10}   {total_de/graph_ms*100:.1f}% of graph")

    # Orchestrators (dispatch overhead for parallel kernels)
    if cats['orchestrators']:
        print(f"\n  ORCHESTRATORS ({len(cats['orchestrators'])} nodes)")
        print(f"  {'Name':<35} {'D+E(ms)':>10} {'Elements':>8} {'EmptyRL(ms)':>12}")
        print(f"  {'-'*35} {'-'*10} {'-'*8} {'-'*12}")
        total_de = 0
        for n in sorted(cats['orchestrators'], key=lambda x: -x['dequeue_exec_ms']):
            total_de += n['dequeue_exec_ms']
            print(f"  {n['name']:<35} {n['dequeue_exec_ms']:>10.3f} {n['elements']:>8} {n['empty_ms']:>12.3f}")
        print(f"  {'TOTAL':<35} {total_de:>10.3f} {'':>8} {'':>12}  {total_de/graph_ms*100:.1f}% of graph")

    # Collectors
    if cats['collectors']:
        print(f"\n  COLLECTORS ({len(cats['collectors'])} nodes)")
        print(f"  {'Name':<35} {'D+E(ms)':>10} {'Elements':>8} {'EmptyRL(ms)':>12}")
        print(f"  {'-'*35} {'-'*10} {'-'*8} {'-'*12}")
        total_de = 0
        for n in sorted(cats['collectors'], key=lambda x: -x['dequeue_exec_ms']):
            total_de += n['dequeue_exec_ms']
            print(f"  {n['name']:<35} {n['dequeue_exec_ms']:>10.3f} {n['elements']:>8} {n['empty_ms']:>12.3f}")
        print(f"  {'TOTAL':<35} {total_de:>10.3f} {'':>8} {'':>12}  {total_de/graph_ms*100:.1f}% of graph")

    # Loop control
    if cats['loop_control']:
        print(f"\n  LOOP CONTROL ({len(cats['loop_control'])} nodes)")
        print(f"  {'Name':<35} {'D+E(ms)':>10} {'Elements':>8}")
        print(f"  {'-'*35} {'-'*10} {'-'*8}")
        for n in sorted(cats['loop_control'], key=lambda x: -x['dequeue_exec_ms']):
            print(f"  {n['name']:<35} {n['dequeue_exec_ms']:>10.3f} {n['elements']:>8}")

    # Other (uncategorized)
    if cats['other']:
        print(f"\n  OTHER ({len(cats['other'])} nodes)")
        print(f"  {'Name':<35} {'D+E(ms)':>10} {'Elements':>8}")
        print(f"  {'-'*35} {'-'*10} {'-'*8}")
        for n in sorted(cats['other'], key=lambda x: -x['dequeue_exec_ms']):
            print(f"  {n['name']:<35} {n['dequeue_exec_ms']:>10.3f} {n['elements']:>8}")

    # Overhead summary
    barrier_de = sum(n['dequeue_exec_ms'] for n in cats['barrier_states'])
    orch_de = sum(n['dequeue_exec_ms'] for n in cats['orchestrators'])
    coll_de = sum(n['dequeue_exec_ms'] for n in cats['collectors'])
    loop_de = sum(n['dequeue_exec_ms'] for n in cats['loop_control'])
    other_de = sum(n['dequeue_exec_ms'] for n in cats['other'])
    total_state_de = sum(n['dequeue_exec_ms'] for n in nodes)
    kernel_est = graph_ms - total_state_de  # estimated parallel kernel time

    print(f"\n  OVERHEAD SUMMARY")
    print(f"  {'Category':<35} {'D+E(ms)':>12} {'% of Graph':>12}")
    print(f"  {'-'*35} {'-'*12} {'-'*12}")
    print(f"  {'Barrier States (sequential)':<35} {barrier_de:>12.3f} {barrier_de/graph_ms*100:>11.1f}%")
    print(f"  {'Orchestrators (dispatch)':<35} {orch_de:>12.3f} {orch_de/graph_ms*100:>11.1f}%")
    print(f"  {'Collectors (gather)':<35} {coll_de:>12.3f} {coll_de/graph_ms*100:>11.1f}%")
    print(f"  {'Loop Control':<35} {loop_de:>12.3f} {loop_de/graph_ms*100:>11.1f}%")
    if cats['other']:
        print(f"  {'Other':<35} {other_de:>12.3f} {other_de/graph_ms*100:>11.1f}%")
    print(f"  {'-'*35} {'-'*12} {'-'*12}")
    print(f"  {'Total State Overhead':<35} {total_state_de:>12.3f} {total_state_de/graph_ms*100:>11.1f}%")
    print(f"  {'Kernel Execution (estimated)':<35} {kernel_est:>12.3f} {kernel_est/graph_ms*100:>11.1f}%")
    print(f"  {'Graph Total':<35} {graph_ms:>12.3f} {'100.0%':>12}")

    # Top-15 nodes by D+E
    print(f"\n  TOP 15 NODES BY DEQUEUE+EXEC TIME")
    print(f"  {'Name':<35} {'D+E(ms)':>10} {'% of Graph':>12} {'Category'}")
    print(f"  {'-'*35} {'-'*10} {'-'*12} {'-'*20}")
    for n in sorted(nodes, key=lambda x: -x['dequeue_exec_ms'])[:15]:
        cat = 'barrier'
        if n in cats['orchestrators']: cat = 'orchestrator'
        elif n in cats['collectors']: cat = 'collector'
        elif n in cats['loop_control']: cat = 'loop'
        elif n in cats['other']: cat = 'other'
        print(f"  {n['name']:<35} {n['dequeue_exec_ms']:>10.3f} {n['dequeue_exec_ms']/graph_ms*100:>11.1f}%   {cat}")


def find_dot_files(base_dir):
    """Find all test directories with dot files."""
    results = []
    run_dir = os.path.join(base_dir, 'run')
    if not os.path.isdir(run_dir):
        return results
    for entry in sorted(os.listdir(run_dir)):
        dot_path = os.path.join(run_dir, entry, 'fds_hh_graph.dot')
        if os.path.exists(dot_path):
            results.append((entry, dot_path))
    return results


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # Use command-line args or auto-discover
    if len(sys.argv) > 1:
        test_dirs = []
        for arg in sys.argv[1:]:
            dot_path = os.path.join(base_dir, 'run', arg, 'fds_hh_graph.dot')
            if os.path.exists(dot_path):
                test_dirs.append((arg, dot_path))
            elif os.path.exists(arg):
                test_dirs.append((os.path.basename(os.path.dirname(arg)), arg))
            else:
                print(f"Warning: {arg} not found")
    else:
        test_dirs = find_dot_files(base_dir)

    if not test_dirs:
        print("No dot files found. Run tests first or specify test names as arguments.")
        return

    for test_name, dot_path in test_dirs:
        graph_time, nodes = parse_dot_file(dot_path)
        print_report(test_name, graph_time, nodes)

    # Cross-test comparison
    if len(test_dirs) > 1:
        print(f"\n{'='*90}")
        print(f"  CROSS-TEST OVERHEAD COMPARISON")
        print(f"{'='*90}")
        print(f"  {'Test Case':<30} {'Graph(s)':>8} {'Barrier(ms)':>12} {'Orch(ms)':>10} {'Coll(ms)':>10} {'Overhead%':>10} {'Kernel%':>10}")
        print(f"  {'-'*30} {'-'*8} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

        for test_name, dot_path in test_dirs:
            graph_time, nodes = parse_dot_file(dot_path)
            cats = categorize_nodes(nodes)
            graph_ms = graph_time * 1000.0 if graph_time > 0 else 0.001
            barrier_de = sum(n['dequeue_exec_ms'] for n in cats['barrier_states'])
            orch_de = sum(n['dequeue_exec_ms'] for n in cats['orchestrators'])
            coll_de = sum(n['dequeue_exec_ms'] for n in cats['collectors'])
            total_de = sum(n['dequeue_exec_ms'] for n in nodes)
            kernel_pct = (graph_ms - total_de) / graph_ms * 100
            overhead_pct = total_de / graph_ms * 100
            print(f"  {test_name:<30} {graph_time:>8.3f} {barrier_de:>12.3f} {orch_de:>10.3f} {coll_de:>10.3f} {overhead_pct:>9.1f}% {kernel_pct:>9.1f}%")


if __name__ == '__main__':
    main()
