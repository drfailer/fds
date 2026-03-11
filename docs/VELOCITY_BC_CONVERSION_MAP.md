# VELOCITY_BC to VELOCITY_BC_KERNEL Conversion Map

## Pattern: Following WallBC Example

**Signature:**
```fortran
SUBROUTINE VELOCITY_BC_KERNEL(M, NM, T, APPLY_TO_ESTIMATED_VARIABLES)
  TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
  INTEGER, INTENT(IN) :: NM
  REAL(EB), INTENT(IN) :: T
  LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
```

**Setup:**
- Call `POINT_TO_MESH(NM)` at beginning (sets up module pointers for CC_VELOCITY callees)
- Use explicit `M%...` for mesh-specific arrays
- Keep global arrays/flags as-is

## Variable Substitutions

### Mesh-Specific Arrays (Replace with M%...)

| Original | Replacement | Lines |
|----------|-------------|-------|
| `US, VS, WS` | `M%US, M%VS, M%WS` | 747-749 |
| `RHOS` | `M%RHOS` | 750 |
| `ZZS` | `M%ZZS` | 751 |
| `U, V, W` | `M%U, M%V, M%W` | 753-755 |
| `RHO` | `M%RHO` | 756 |
| `ZZ` | `M%ZZ` | 757 |
| `N_EXTERNAL_WALL_CELLS` | `M%N_EXTERNAL_WALL_CELLS` | 762 |
| `WALL(...)` | `M%WALL(...)` | Throughout |
| `EXTERNAL_WALL(...)` | `M%EXTERNAL_WALL(...)` | Throughout |
| `OMESH(...)` | `M%OMESH(...)` | 767-774, 1307 |
| `BOUNDARY_COORD(...)` | `M%BOUNDARY_COORD(...)` | Throughout |
| `DRAG_UVWMAX` | `M%DRAG_UVWMAX` | 803, 1299 |
| `EDGE_COUNT(NM)` | `M%EDGE_COUNT` | 807 |
| `EDGE(...)` | `M%EDGE(...)` | 809 onwards |
| `CELL(...)` | `M%CELL(...)` | Throughout |
| `DY(...)` | `M%DY(...)` | 853 |
| `DZ(...)` | `M%DZ(...)` | 854, 860 |
| `DX(...)` | `M%DX(...)` | 861, 867-868 |
| `IBAR` | `M%IBAR` | Throughout |
| `JBAR` | `M%JBAR` | Throughout |
| `KBAR` | `M%KBAR` | Throughout |
| `IBP1` | `M%IBP1` | Throughout |
| `JBP1` | `M%JBP1` | Throughout |
| `KBP1` | `M%KBP1` | Throughout |
| `BOUNDARY_PROP1(...)` | `M%BOUNDARY_PROP1(...)` | Throughout |
| `MU(...)` | `M%MU(...)` | 1171, 1274 |
| `TMP(...)` | `M%TMP(...)` | 1274 |
| `ZC(...)` | `M%ZC(...)` | 1224, 1229 |
| `T_USED(...)` | `M%T_USED(...)` | 1508 |

### Global Arrays/Flags (Keep as-is)

| Variable | Note |
|----------|------|
| `VENTS(...)` | Global array |
| `SURFACE(...)` | Global array |
| `U_WIND(...)` | Global array |
| `V_WIND(...)` | Global array |
| `W_WIND(...)` | Global array |
| `PREDICTOR` | Global flag |
| `CORRECTOR` | Global flag |
| `CC_IBM` | Global flag |
| `OPEN_WIND_BOUNDARY` | Global flag |
| `SOLID_PHASE_ONLY` | Global constant |
| `PERIODIC_TEST` | Global variable |
| `T_BEGIN` | Global constant |
| `MU_RSQMW_Z, RSQ_MW_Z` | Global arrays |
| `I_MAX_TEMP` | Global constant |

## Special Cases

### Pointer Assignments (Lines 746-758)
```fortran
! OLD:
IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   ZZP => ZZS
ELSE
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   ZZP => ZZ
ENDIF

! NEW:
IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
   RHOP => M%RHOS
   ZZP => M%ZZS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
   RHOP => M%RHO
   ZZP => M%ZZ
ENDIF
```

### OMESH Access (Lines 767-774, 1307, etc.)
```fortran
! OLD: OMESH(EWC%NOM)%US
! NEW: M%OMESH(EWC%NOM)%US
```

### WALL/CELL/EDGE/BOUNDARY Pointers
All remain as pointer assignments, but point to M%... arrays:
```fortran
! OLD: WC =>WALL(IW)
! NEW: WC => M%WALL(IW)

! OLD: BC => BOUNDARY_COORD(WC%BC_INDEX)
! NEW: BC => M%BOUNDARY_COORD(WC%BC_INDEX)
```

## Wrapper Function

Keep old `VELOCITY_BC` as wrapper:
```fortran
SUBROUTINE VELOCITY_BC(T, NM, APPLY_TO_ESTIMATED_VARIABLES)
  REAL(EB), INTENT(IN) :: T
  INTEGER, INTENT(IN) :: NM
  LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES

  CALL VELOCITY_BC_KERNEL(MESHES(NM), NM, T, APPLY_TO_ESTIMATED_VARIABLES)

END SUBROUTINE VELOCITY_BC
```

## Total Substitutions

Estimated ~150+ substitutions needed:
- WALL: ~30 occurrences
- CELL: ~40 occurrences
- EDGE: ~15 occurrences
- BOUNDARY_COORD: ~10 occurrences
- BOUNDARY_PROP1: ~15 occurrences
- Grid variables (IBAR, JBAR, etc.): ~25 occurrences
- Velocity arrays (U, V, W, etc.): ~20 occurrences
- Other mesh arrays: ~20 occurrences
