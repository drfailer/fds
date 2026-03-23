// UNUSED — Block decomposition disabled. Kept for reference.
#ifndef MESH_BLOCK_STATE_H
#define MESH_BLOCK_STATE_H

#include <hedgehog/hedgehog.h>
#include <unordered_map>
#include "../data/mesh_data.h"
#include "../data/mesh_block_data.h"
#include "../fds_fortran_interface.h"

/// Decomposes a MeshData token into multiple MeshBlockData tokens along K.
/// Each MeshData produces ceil(KBAR / blockSize) blocks, dispatched immediately.
class MeshBlockDecomposeState
    : public hh::AbstractState<1, MeshData, MeshBlockData> {
public:
    /// @param numBlocks Target number of blocks per mesh
    explicit MeshBlockDecomposeState(int numBlocks)
        : hh::AbstractState<1, MeshData, MeshBlockData>(),
          numBlocks_(std::max(1, numBlocks)) {}

    void execute(std::shared_ptr<MeshData> data) override {
        int kbar = fds_get_kbar(data->nm);
        int bs = std::max(1, (kbar + numBlocks_ - 1) / numBlocks_);
        int total = (kbar + bs - 1) / bs;

        for (int b = 0; b < total; ++b) {
            int k1 = b * bs + 1;
            int k2 = std::min((b + 1) * bs, kbar);
            this->addResult(std::make_shared<MeshBlockData>(
                data->nm, k1, k2, data->t, data->dt, data->phase, total, data));
        }
    }

private:
    int numBlocks_;
};

/// Reassembles MeshBlockData tokens back into MeshData.
/// Collects all blocks for a given mesh (identified by nm) and emits
/// the original MeshData when all blocks have arrived.
class MeshBlockReassembleState
    : public hh::AbstractState<1, MeshBlockData, MeshData> {
public:
    MeshBlockReassembleState() = default;

    void execute(std::shared_ptr<MeshBlockData> block) override {
        int nm = block->nm;
        auto &entry = entries_[nm];
        if (entry.count == 0) {
            entry.expected = block->totalBlocks;
            entry.meshData = block->originalMeshData;
        }
        entry.count++;
        if (entry.count == entry.expected) {
            this->addResult(entry.meshData);
            entries_.erase(nm);
        }
    }

private:
    struct Entry {
        int count = 0;
        int expected = 0;
        std::shared_ptr<MeshData> meshData;
    };
    std::unordered_map<int, Entry> entries_;
};

#endif // MESH_BLOCK_STATE_H
