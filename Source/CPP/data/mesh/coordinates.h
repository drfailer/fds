#ifndef COORDINATES_H
#define COORDINATES_H
#include <cstddef>

namespace mesh {

struct Coordinates {
   double DXI;       ///< \f$ \delta \xi = (x_I-x_0)/I \f$
   double DETA;      ///< \f$ \delta \eta = (y_J-y_0)/J \f$
   double DZETA;     ///< \f$ \delta \zeta = (z_K-z_0)/K \f$
   double RDXI;      ///< \f$ 1/ \delta \xi \f$
   double RDETA;     ///< \f$ 1/ \delta \eta \f$
   double RDZETA;    ///< \f$ 1/ \delta \zeta \f$
   double DXMIN;     ///< \f$ \min_i \delta x_i \f$
   double DXMAX;     ///< \f$ \max_i \delta x_i \f$
   double DYMIN;     ///< \f$ \min_j \delta y_j \f$
   double DYMAX;     ///< \f$ \max_j \delta y_j \f$
   double DZMIN;     ///< \f$ \min_k \delta z_k \f$
   double DZMAX;     ///< \f$ \max_k \delta z_k \f$
   double XS;        ///< Lower extent of mesh x coordinate, \f$ x_0 \f$
   double XF;        ///< Upper extent of mesh x coordinate, \f$ x_I \f$
   double YS;        ///< Lower extent of mesh y coordinate, \f$ y_0 \f$
   double YF;        ///< Upper extent of mesh y coordinate, \f$ y_J \f$
   double ZS;        ///< Lower extent of mesh z coordinate, \f$ z_0 \f$
   double ZF;        ///< Upper extent of mesh z coordinate, \f$ z_K \f$
   double RDXINT;    ///< \f$ 500/\delta \xi \f$
   double RDYINT;    ///< \f$ 500/\delta \eta \f$
   double RDZINT;    ///< \f$ 500/\delta \zeta \f$
   double CELL_SIZE; ///< Approximate cell size, \f$ (\delta\xi\,\delta\eta\,\delta\zeta)^{1/3} \f$
   double *R;      ///< Radial coordinate, \f$ r_i \f$, for CYLINDRICAL geometry
   double *RC;     ///< Radial coordinate, cell center, \f$ (r_i+r_{i-1})/2 \f$
   double *RRN;    ///< \f$ 2/(r_i+r_{i-1}) \f$
   double *X;      ///< Position of forward x face of cell (I,J,K), \f$ x_i \f$
   double *Y;      ///< Position of forward y face of cell (I,J,K), \f$ y_j \f$
   double *Z;      ///< Position of forward z face of cell (I,J,K), \f$ z_k \f$
   double *XC;     ///< x coordinate of cell center, \f$ (x_i+x_{i-1})/2 \f$
   double *YC;     ///< y coordinate of cell center, \f$ (y_j+y_{j-1})/2 \f$
   double *ZC;     ///< z coordinate of cell center, \f$ (z_k+z_{k-1})/2 \f$
   double *HX;     ///< Grid stretch factor, \f$ (x_i-x_{i-1})/\delta \xi \f$
   double *HY;     ///< Grid stretch factor, \f$ (y_j-y_{j-1})/\delta \eta \f$
   double *HZ;     ///< Grid stretch factor, \f$ (z_k-z_{k-1})/\delta \zeta \f$
   double *DX;     ///< \f$ \delta x_i = x_i-x_{i-1} \f$
   double *DY;     ///< \f$ \delta y_j = y_j-y_{j-1} \f$
   double *DZ;     ///< \f$ \delta z_k = z_k-z_{k-1} \f$
   double *RDX;    ///< \f$ 1/\delta x_i \f$
   double *RDY;    ///< \f$ 1/\delta y_j \f$
   double *RDZ;    ///< \f$ 1/\delta z_k \f$
   double *DXN;    ///< \f$ (x_i+x_{i+1})/2 \f$
   double *DYN;    ///< \f$ (y_j+y_{j+1})/2 \f$
   double *DZN;    ///< \f$ (z_k+z_{k+1})/2 \f$
   double *RDXN;   ///< \f$ 2/(x_i+x_{i+1}) \f$
   double *RDYN;   ///< \f$ 2/(y_j+y_{j+1}) \f$
   double *RDZN;   ///< \f$ 2/(z_k+z_{k+1}) \f$
   double *CELLSI; ///< Array used to locate the cell index of \f$ x \f$
   double *CELLSJ; ///< Array used to locate the cell index of \f$ y \f$
   double *CELLSK; ///< Array used to locate the cell index of \f$ z \f$
   size_t CELLSI_LO, CELLSI_HI; ///< hold CELLSI array bounds
   size_t CELLSJ_LO, CELLSJ_HI; ///< hold CELLSJ array bounds
   size_t CELLSK_LO, CELLSK_HI; ///< hold CELLSK array bounds
   double *XPLT; ///< 4 byte real array holding \f$ x \f$ mesh coordinates
   double *YPLT; ///< 4 byte real array holding \f$ y \f$ mesh coordinates
   double *ZPLT; ///< 4 byte real array holding \f$ z \f$ mesh coordinates
};

} // end namespace mesh

#endif
