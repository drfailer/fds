#ifndef EXCHANGE_DISPATCH_H
#define EXCHANGE_DISPATCH_H

#include <memory>
#include "../data/mesh_data.h"
#include "pack_data.h"
#include "exchange_bundle.h"

/// Compile-time pair mapping input MeshState K to output MeshState Out.
template<MeshState K, MeshState Out = K>
struct ExchKind {
    static constexpr MeshState input = K;
    static constexpr MeshState output = Out;
};

// ---------------------------------------------------------------------------
// Recursive CRTP dispatch helpers.
//
// Each generates virtual execute() overrides for every MeshState in a
// parameter pack, delegating to Derived's handler via static_cast (CRTP).
//
// Usage:  class MyTask : public MeshDataDispatch<MyTask, HHBase, K1, K2> { ... };
//         MyTask must provide:  template<MeshState K> void handleMeshData(...)
// ---------------------------------------------------------------------------

/// Dispatch layer: generates execute(shared_ptr<MeshData<K>>) for each K.
template<typename Derived, typename HHBase, MeshState... Ks>
struct MeshDataDispatch;

template<typename Derived, typename HHBase>
struct MeshDataDispatch<Derived, HHBase> : HHBase {
    using HHBase::HHBase;
};

template<typename Derived, typename HHBase, MeshState K, MeshState... Rest>
struct MeshDataDispatch<Derived, HHBase, K, Rest...>
    : MeshDataDispatch<Derived, HHBase, Rest...> {
    using MeshDataDispatch<Derived, HHBase, Rest...>::MeshDataDispatch;
    void execute(std::shared_ptr<MeshData<K>> data) override {
        static_cast<Derived*>(this)->template handleMeshData<K>(std::move(data));
    }
};

/// Dispatch layer: generates execute(shared_ptr<PackData<K>>) for each K.
template<typename Derived, typename HHBase, MeshState... Ks>
struct PackDataDispatch;

template<typename Derived, typename HHBase>
struct PackDataDispatch<Derived, HHBase> : HHBase {
    using HHBase::HHBase;
};

template<typename Derived, typename HHBase, MeshState K, MeshState... Rest>
struct PackDataDispatch<Derived, HHBase, K, Rest...>
    : PackDataDispatch<Derived, HHBase, Rest...> {
    using PackDataDispatch<Derived, HHBase, Rest...>::PackDataDispatch;
    void execute(std::shared_ptr<PackData<K>> data) override {
        static_cast<Derived*>(this)->template handlePackData<K>(std::move(data));
    }
};

/// Dispatch layer: generates execute(shared_ptr<ExchangeBundle<K>>) for each K.
template<typename Derived, typename HHBase, MeshState... Ks>
struct BundleDispatch;

template<typename Derived, typename HHBase>
struct BundleDispatch<Derived, HHBase> : HHBase {
    using HHBase::HHBase;
};

template<typename Derived, typename HHBase, MeshState K, MeshState... Rest>
struct BundleDispatch<Derived, HHBase, K, Rest...>
    : BundleDispatch<Derived, HHBase, Rest...> {
    using BundleDispatch<Derived, HHBase, Rest...>::BundleDispatch;
    void execute(std::shared_ptr<ExchangeBundle<K>> data) override {
        static_cast<Derived*>(this)->template handleBundle<K>(std::move(data));
    }
};

#endif // EXCHANGE_DISPATCH_H
