#ifndef MESH_OMESH_H
#define MESH_OMESH_H

#include "storage.h"
namespace mesh {

struct OMesh {
  double *MU, *RHO, *RHOS, *TMP, *U, *V, *W, *US, *VS, *WS, *H, *HS, *FVX, *FVY,
      *FVZ, *D, *DS, *KRES, *IL_S, *IL_R, *Q;
  double *ZZ, *ZZS;
  int *IIO_R, *JJO_R, *KKO_R, *IOR_R, *IIO_S, *JJO_S, *KKO_S, *IOR_S;
  int I_MIN_R = -10, I_MAX_R = -10, J_MIN_R = -10, J_MAX_R = -10, K_MIN_R = -10,
      K_MAX_R = -10, NIC_R = 0, I_MIN_S = -10, I_MAX_S = -10, J_MIN_S = -10,
      J_MAX_S = -10, K_MIN_S = -10, K_MAX_S = -10, NIC_S = 0;
  int INTEGER_SEND_BUFFER[7], INTEGER_RECV_BUFFER[7];
  double *REAL_SEND_PKG1, *REAL_SEND_PKG2, *REAL_SEND_PKG3, *REAL_SEND_PKG4,
      *REAL_SEND_PKG5, *REAL_SEND_PKG6, *REAL_SEND_PKG7, *REAL_SEND_PKG8,
      *REAL_RECV_PKG1, *REAL_RECV_PKG2, *REAL_RECV_PKG3, *REAL_RECV_PKG4,
      *REAL_RECV_PKG5, *REAL_RECV_PKG6, *REAL_RECV_PKG7, *REAL_RECV_PKG8;
  int N_EXTERNAL_OBST = 0, N_INTERNAL_OBST = 0;
  Storage PARTICLE_SEND_BUFFER, PARTICLE_RECV_BUFFER;
  Storage WALL_SEND_BUFFER, WALL_RECV_BUFFER;
  Storage THIN_WALL_SEND_BUFFER, THIN_WALL_RECV_BUFFER;

  // CC_IBM data exchange arrays:
  int NICC_S[2] = {0}, NICC_R[2] = {0}, NICF_S[2] = {0}, NICF_R[2] = {0},
      NLKF_S = 0, NLKF_R = 0, NFCC_S[2] = {0}, NFCC_R[2] = {0}, NCC_INT_R = 0,
      NFEP_R[5] = {0}, NFEP_R_G = 0;
  double *REAL_SEND_PKG11, *REAL_SEND_PKG112, *REAL_SEND_PKG12,
      *REAL_SEND_PKG13, *REAL_RECV_PKG11, *REAL_RECV_PKG112, *REAL_RECV_PKG12,
      *REAL_RECV_PKG13;
  int *ICC_UNKZ_CT_S, *ICC_UNKZ_CC_S, *ICC_UNKZ_CT_R, *ICC_UNKZ_CC_R,
      *ICF_UFFB_CF_S, *ICF_UFFB_CF_R;
  int *UNKZ_CT_S, *UNKZ_CC_S, *UNKZ_CT_R, *UNKZ_CC_R;

  // Face variables data (velocities):
  int *IIO_FC_R, *JJO_FC_R, *KKO_FC_R, *AXS_FC_R, *IIO_FC_S, *JJO_FC_S,
      *KKO_FC_S, *AXS_FC_S;
  int *IIO_CC_R, *JJO_CC_R, *KKO_CC_R, *IIO_CC_S, *JJO_CC_S, *KKO_CC_S;
  int *IIO_LF_R, *JJO_LF_R, *KKO_LF_R, *AXS_LF_R, *IIO_LF_S, *JJO_LF_S,
      *KKO_LF_S, *AXS_LF_S;
  int *IFEP_R_1, *IFEP_R_2, *IFEP_R_3, *IFEP_R_4, *IFEP_R_5;
  double *U_LNK, *V_LNK, *W_LNK;

  // Level Set
  double *PHI_LS, *PHI1_LS, *U_LS, *V_LS, *Z_LS, *TOA;
  double *REAL_SEND_PKG14, *REAL_RECV_PKG14;
};

} // end namespace mesh

#endif
