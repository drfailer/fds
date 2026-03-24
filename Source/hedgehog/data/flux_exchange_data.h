#ifndef FLUX_EXCHANGE_DATA_H
#define FLUX_EXCHANGE_DATA_H

/// Data token for per-neighbor flux exchange.
///
/// Represents a completed flux copy from sourceNM to destNM.
/// For same-process: the copy was already done in FluxPackState.
/// For cross-process (future): carries a buffer for MPI transfer.
struct FluxExchangeData {
    int sourceNM;  ///< Mesh that produced the flux data
    int destNM;    ///< Mesh that needs the flux data

    FluxExchangeData() : sourceNM(0), destNM(0) {}
    FluxExchangeData(int src, int dst) : sourceNM(src), destNM(dst) {}
};

#endif // FLUX_EXCHANGE_DATA_H
