#ifndef EXCHANGE_STRATEGY_H
#define EXCHANGE_STRATEGY_H

#include "../fds_fortran_interface.h"

/// Compile-time trait: determines what data to push/pull for flux exchange (CODE 5).
///
/// Push copies the full 3D slab of FVX, FVY, FVZ, and H (or HS) from
/// MESHES(NM) into a flat buffer.  Pull copies from the buffer into
/// MESHES(NOM)%OMESH(NM) arrays.  This matches the data transferred by
/// MESH_EXCHANGE_FLUX_NEIGHBOR_TS (full slab copy).
///
/// Future exchange codes (species, velocity, etc.) define their own traits with
/// the same interface and different Fortran bindings.
struct FluxExchangeStrategy {
    /// Push: copy full slab from MESHES(NM) source arrays into a flat buffer.
    static void pushToBuffer(int nm, int nom, double *buf, int size) {
        fds_flux_push_slab(nm, nom, buf, size);
    }

    /// Pull: copy full slab from buffer into MESHES(NOM)%OMESH(NM) arrays.
    static void pullFromBuffer(int nm, int nom, const double *buf, int size) {
        fds_flux_pull_slab(nm, nom, buf, size);
    }

    /// Buffer size (doubles) for one (source, dest) pair: 4 * NI * NJ * NK.
    static int bufferSize(int nm, int nom) {
        return fds_flux_slab_size(nm, nom);
    }

    /// Max buffer size across all local mesh pairs (for memory pool pre-alloc).
    static int maxBufferSize() {
        return fds_flux_max_slab_size();
    }

    /// Buffer size using receiver-side data (safe when nm is remote).
    /// Arguments: nom = local dest mesh, nm = source mesh (possibly remote).
    static int bufferSizeRecv(int nom, int nm) {
        return fds_flux_slab_size_recv(nom, nm);
    }

    /// Pull using receiver-side data (safe when nm is remote).
    /// Arguments: nom = local dest mesh, nm = source mesh (possibly remote).
    static void pullFromBufferRecv(int nom, int nm, const double *buf, int size) {
        fds_flux_pull_slab_recv(nom, nm, buf, size);
    }
};

#endif // EXCHANGE_STRATEGY_H
