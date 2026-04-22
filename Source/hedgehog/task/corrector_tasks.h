#ifndef CORRECTOR_TASKS_H
#define CORRECTOR_TASKS_H

#include <hedgehog/hedgehog.h>
#include "../data/mesh_data.h"
#include "../fds_fortran_interface.h"

/// Task C1: Mass finite diff + compute viscosity + density (corrector)
class CorrStep1Task : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrStep1Task(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrStep1", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_compute_viscosity(data->nm, 1);  // estimated=true
        fds_mass_finite_differences(data->nm);
        fds_density(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrStep1Task>(this->numberThreads());
    }
};

/// Task C2: Viscosity BC + velocity flux (corrector)
class CorrDivSetupTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrDivSetupTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrDivSetup", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_set_baroclinic_false(data->nm);
        fds_viscosity_bc(data->nm, 1);  // estimated=true
        fds_velocity_flux(data->t, data->dt, data->nm, 1);
        fds_agglomeration(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrDivSetupTask>(this->numberThreads());
    }
};

/// Task C3: Condensation/evaporation
class CorrCondensTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrCondensTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrCondens", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_condensation(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrCondensTask>(this->numberThreads());
    }
};

/// Task C4: Particle mass/energy + move + momentum transfer
class CorrParticleTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrParticleTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrParticle", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_particle_mass_energy(data->t, data->dt, data->nm);
        fds_move_particles(data->t, data->dt, data->nm);
        fds_particle_momentum(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrParticleTask>(this->numberThreads());
    }
};

/// Task C5: Wall BC (corrector)
class CorrWallBCTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrWallBCTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrWallBC", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_wall_bc(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrWallBCTask>(this->numberThreads());
    }
};

/// Task C7: Combustion BC + divergence part 1
class CorrDivPart1Task : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrDivPart1Task(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrDivPart1", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_combustion_bc(data->nm);
        fds_divergence_part_1(data->t, data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrDivPart1Task>(this->numberThreads());
    }
};

/// Task C8: Divergence part 2 (corrector)
class CorrDivPart2Task : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrDivPart2Task(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrDivPart2", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_divergence_part_2(data->dt, data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrDivPart2Task>(this->numberThreads());
    }
};

/// Task C9: Velocity corrector + check divergence
class CorrVelocityTask : public hh::AbstractTask<1, MeshData<>, MeshData<>> {
public:
    explicit CorrVelocityTask(size_t nThreads = 1)
        : hh::AbstractTask<1, MeshData<>, MeshData<>>("CorrVelocity", nThreads) {}

    void execute(std::shared_ptr<MeshData<>> data) override {
        fds_velocity_corrector(data->t, data->dt, data->nm);
        fds_check_divergence(data->nm);
        this->addResult(data);
    }

    std::shared_ptr<hh::AbstractTask<1, MeshData<>, MeshData<>>> copy() override {
        return std::make_shared<CorrVelocityTask>(this->numberThreads());
    }
};

// CorrFinalTask replaced by CorrFinal sub-graph (Pattern B) in velocity_bc_subgraph.h
// DUMP_MESH_OUTPUTS remains in TimestepState (must happen after global UPDATE_CONTROLS)

#endif // CORRECTOR_TASKS_H
