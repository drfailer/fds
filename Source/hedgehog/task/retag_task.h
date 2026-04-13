#ifndef RETAG_TASK_H
#define RETAG_TASK_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"

/// Generic type-conversion task: MeshData<From> -> MeshData<To>.
///
/// Used at subgraph boundaries to bridge between the external pipeline
/// (MeshData<Default>) and typed subgraph inputs/outputs
/// (e.g. MeshData<PredictorPressure>).
///
/// Single-threaded — the conversion is a trivial field copy.
template<MeshState From, MeshState To>
class RetagTask
    : public hh::AbstractTask<1, MeshData<From>, MeshData<To>> {
public:
    explicit RetagTask(std::string name = "Retag")
        : hh::AbstractTask<1, MeshData<From>, MeshData<To>>(
              std::move(name), 1) {}

    void execute(std::shared_ptr<MeshData<From>> md) override {
        this->addResult(md->template retag<To>());
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<From>, MeshData<To>>>
    copy() override {
        return std::make_shared<RetagTask<From, To>>(this->name());
    }
};

#endif // RETAG_TASK_H
