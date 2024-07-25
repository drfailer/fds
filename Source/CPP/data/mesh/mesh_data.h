#ifndef MESH_DATA_H
#define MESH_DATA_H
#include <cstddef>
#include "../edge/edge.h"
#include "adv.h"
#include "ambient.h"
#include "cc/cc_edge.h"
#include "cc/cc_ibm_arrays.h"
#include "cc/cc_inbcf_area.h"
#include "cc/cc_rcface.h"
#include "cc/cc_regface.h"
#include "cc/cc_regfacez.h"
#include "cells.h"
#include "cface.h"
#include "coordinates.h"
#include "density.h"
#include "dif.h"
#include "dimensions.h"
#include "divergence.h"
#include "droplet.h"
#include "lagrangian_particle.h"
#include "lumped_species.h"
#include "momentum_flux_terms.h"
#include "obstruction.h"
#include "../constants.h"
#include "omesh.h"
#include "particle_drag.h"
#include "patch.h"
#include "pressure.h"
#include "rad_file.h"
#include "radiation.h"
#include "slice.h"
#include "velocity_components.h"
#include "vents.h"
#include "viscosity.h"
#include "wall.h"
#include "wind.h"
#include "work.h"
#include "zone_mesh.h"

namespace mesh {

struct Mesh {
  Dimensions dimensions;
  Coordinates coordinates;

  Cells cells;
  edge::Edge *EDGE;

  VelocityComponents velocityComponents;
  Divergence divergence;
  MomentumFluxTerms momentumFluxTerms;
  ParticuleDrag particleDrag;
  Density density;
  Pressure pressure;
  Radiation radiation;
  LumpedSpecies lumpedSpecies;
  Droplet droplet;
  Viscosity viscosity;
  Ambient ambient;
  Wind wind;
  Adv adv;
  Dif dif;
  Work work;

  double *H;    ///< \f$ \tilde{p}_{ijk}/\rho_{ijk} + |\mathbf{u}|^2_{ijk}/2 \f$
  double *HS;   ///< H estimated at next time step
  double *KRES; ///< Resolved kinetic energy, \f$ |\mathbf{u}|^2_{ijk}/2 \f$
  double *TMP;  ///< Gas temperature, \f$ T_{ijk} \f$ (K)
  double *Q;    ///< Heat release rate per unit volume, \f$ \dot{q}_{ijk}''' \f$
  double *RSUM; ///< \f$ R_0 \sum_\alpha Z_{\alpha,ijk}/W_\alpha \f$
  double *CSD2; ///< \f$ C_s \Delta^2 \f$ in Smagorinsky turbulence expression
  double *CHEM_SUBIT;  ///< Number of chemistry sub-iterations
  double *MIX_TIME;    ///< Mixing-controlled combustion reaction time (s)
  double *STRAIN_RATE; ///< Strain rate \f$ |S|_{ijk} \f$ (1/s)
  double *D_Z_MAX;     ///< \f$ \max D_\alpha \f$
  double *Q_DOT_PPP_S; ///< Heat release rate per unit volume in 3D pyrolysis model

  double *REAC_SOURCE_TERM; ///< \f$ \dot{m}_{\alpha,ijk}''' \f$
  double *DEL_RHO_D_DEL_Z;  ///< \f$ (\nabla \cdot \rho D_\alpha \nabla Z_\alpha)_{ijk} \f$
  double *FX;        ///< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ x \f$ face of cell
  double *FY;        ///< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ y \f$ face of cell
  double *FZ;        ///< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ z \f$ face of cell
  double *Q_REAC;    ///< \f$ \dot{q}_{ijk}''' \f$ for a specified reaction
  double *M_DOT_PPP; ///< Mass source term, \f$ \dot{m}_{\alpha,ijk}''' \f$

  double POIS_PTB, POIS_ERR;
  double *SAVE1, *SAVE2, *WORK;
  double *BXS, *BXF, *BYS, *BYF, *BZS, *BZF, *BXST, *BXFT, *BYST, *BYFT, *BZST,
      *BZFT;
  size_t LSAVE, LWORK, LBC, MBC, NBC, LBC2, MBC2, NBC2, ITRN, JTRN, KTRN, IPS;

  double *D_PBAR_DT;   ///< \f$ (\partial \overline{p}_m/\partial t)^n \f$
  double *D_PBAR_DT_S; ///< \f$ (\partial \overline{p}_m/\partial t)^* \f$
  double *U_LEAK;
  double *U_DUCT;

  double *R_PBAR; ///< \f$ 1/\overline{p}_m(z,t) \f$
  ZoneMesh *ZONE_MESH;  ///< Stores Pardiso parameters for specific ZONE on a MESH

  double CFL,DIVMX,DIVMN,VN,RESMAX,PART_UVWMAX=0;
  size_t ICFL, JCFL, KCFL, IMX, JMX, KMX, IMN, JMN, KMN, I_VN, J_VN, K_VN, IRM,
      JRM, KRM, DT_RESTRICT_COUNT = 0, DT_RESTRICT_STORE = 0;
  bool CLIP_RHOMIN = false, CLIP_RHOMAX = false;
  bool BAROCLINIC_TERMS_ATTACHED = false;
  char TRNX_ID[LABEL_LENGTH], TRNY_ID[LABEL_LENGTH], TRNZ_ID[LABEL_LENGTH];

  size_t N_NEIGHBORING_MESHES; ///< Number of meshing abutting the current one
  size_t *NEIGHBORING_MESH;    ///< Array listing the indices of neighboring meshes
  size_t *RGB;                 ///< Color indices of the mesh for Smokeview

  size_t N_OBST=0;          ///< Number of obstructions in the mesh
  Obstruction *OBSTRUCTION; ///< Derived type variable holding obstruction information

  size_t N_VENT=0;          ///< Number of vents in the mesh
  Vents *VENTS;             ///< Derived type variable holding vent information

  bool *CONNECTED_MESH;     ///< T or F if cell is within another mesh

  int NREGFACE_H[MAX_DIM]={0};
  cc::CCRegface *REGFACE_IAXIS_H, *REGFACE_JAXIS_H, *REGFACE_KAXIS_H;

  // CC_IBM mesh arrays
  cc::CCIBMArrays ccIBMArrays;

  size_t CC_NREGFACE_Z[MAX_DIM]={0}, CC_NBBREGFACE_Z[MAX_DIM]={0};
  cc::CCRegfaceZ *CC_REGFACE_IAXIS_Z, *CC_REGFACE_JAXIS_Z, *CC_REGFACE_KAXIS_Z;
  size_t CC_NRCFACE_Z=0, CC_NBBRCFACE_Z=0, CC_NRCFACE_H=0;
  int *RCF_H;
  cc::CCRCface *RC_FACE;
  double *RHO_ZZN;

  // CFACE to be used in conjunction with solid side solvers:
  CFace *CFACE;

  // Edges connected to regular Gas faces to receive wall model shear stress and vorticity.
  size_t CC_NRCEDGE=0, CC_NIBEDGE=0;
  cc::CCEdge *CC_RCEDGE, *CC_IBEDGE;

  // Array with maximum height (Z) of geometry intersections with vertical grid lines in the mesh.
  double *GEOM_ZMAX;

  // Arrays for special cut-cells:
  size_t N_SPCELL=0, N_SPCELL_CF=0;
  int *SPCELL_LIST;

  // Linked face velocity arrays:
  double *EWC_UN_LNK, *UN_LNK, *UN_ULNK;

  // Array with boundary cut-face areas pes cut-cell, for use in cut-cell blocking.
  cc::CCINBCFArea *INBCF_AREA;

  // Arrays for cut-cell blocking:
  size_t N_CC_BLOCKED=0;
  int *JBT_CC_BLOCKED;
  double *XYZ_CC_BLOCKED;

  // ...

  Wall wall;

  size_t HT_3D_SWEEP_DIRECTION = 0;

  size_t N_EXTERNAL_CFACE_CELLS = 0, N_INTWALL_CFACE_CELLS = 0,
         INTERNAL_CFACE_CELLS_LB = 0, N_INTERNAL_CFACE_CELLS = 0,
         N_CFACE_CELLS_DIM = 0;

  double BC_CLOCK;

  double *UVW_SAVE, *U_GHOST, *V_GHOST, *W_GHOST;

  OMesh *OMESH;
  LagrangianParticle *LAGRANGIAN_PARTICLE;
  size_t NLP, NLPDIM, PARTICLE_TAG;

  size_t N_SLCF=0;
  Slice *SLICE;

  size_t N_RADF=0;
  RadFile *RAD_FILE;

  size_t N_PATCH, N_BNDF_POINTS;
  Patch *PATCH;

  double *UIID;
  int RAD_CALL_COUNTER, ANGLE_INC_COUNTER;

  int *INTERPOLATED_MESH;

  char STRING[MESH_STRING_LENGTH];
  size_t N_STRINGS, N_STRINGS_MAX;

  int *K_AGL_SLICE;
  int *LS_KLO_TERRAIN, *LS_KHI_TERRAIN, *K_LS, *LS_SURF_INDEX;
  size_t N_TERRAIN_SLCF=0;

  double *FLUX0_LS, *FLUX1_LS, *PHI_LS, *PHI1_LS, *ROS_BACKU, *TOA, *ROS_HEAD,
      *ROS_FLANK, *WIND_EXP, *SR_X_LS, *SR_Y_LS, *U_LS, *V_LS, *Z_LS, *DZTDX,
      *DZTDY, *MAG_ZT, *PHI_WS, *UMF, *THETA_ELPS, *PHI_S, *PHI_S_X, *PHI_S_Y,
      *PHI_W, *LS_WORK1, *LS_WORK2;
};

} // end namespace mesh

#endif
