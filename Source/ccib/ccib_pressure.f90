!  +++++++++++++++++++++++ CC_PRESSURE ++++++++++++++++++++++++++

! Pressure solver routines for the
! cut-cell / immersed-boundary method.

MODULE CC_PRESSURE

USE CC_SCALARS_DATA
USE CC_SCALARS, ONLY: GET_CUTFACE_BAROCLINIC_TORQUE, GET_RCFACE_BAROCLINIC_TORQUE
USE CC_DIVERGENCE, ONLY: GET_FN_DIVERGENCE_CUTCELL, GET_VELOC_DIVERGENCE_CUTCELL
USE CC_PRESSURE_KERNELS
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: ADD_CUTCELL_D_PBAR_DT, ADD_LINKEDCELL_D_PBAR_DT, &
          ADD_CUTCELL_PSUM, ADD_LINKEDCELL_PSUM, &
          ADD_INPLACE_NNZ_H_WHLDOM, COPY_UNST_DM_TO_CART, &
          COPY_CC_UNKH_TO_HS, COPY_CC_MUNKH_TO_UNKH, &
          GET_BOUNDFACE_GEOM_INFO_H, GET_CC_IROW, GET_CC_UNKH, &
          GET_CC_MATRIXGRAPH_H, GET_CFACE_OPEN_BC_COEF, &
          GET_CUTCELL_DDDT, GET_CUTCELL_HP, &
          GET_FH_FROM_PRHS_AND_BCS, &
          GET_H_CUTFACES, GET_H_GUARD_CUTCELL, GET_H_MATRIX_CC, &
          GET_PRES_CFACE_BCS, GET_RCFACES_H, &
          NUMBER_UNKH_CUTCELLS, &
          UNSTRUCTURED_POISSON_RESIDUAL, UNSTRUCTURED_POISSON_RESIDUAL_RC

CONTAINS

! ------------------------ COPY_UNST_DM_TO_CART -------------------------------------

SUBROUTINE COPY_UNST_DM_TO_CART(NM)
! Assumes POINT_TO_MESH(NM) has been called.
INTEGER, INTENT(IN) :: NM
INTEGER :: ICC,JCC,I,J,K,NS
REAL(EB) :: VOL

DO ICC=1,MESHES(NM)%N_CUTCELL_MESH+MESHES(NM)%N_GCCUTCELL_MESH
   CC => CUT_CELL(ICC);  I = CC%IJK(IAXIS); J = CC%IJK(JAXIS);  K = CC%IJK(KAXIS)
   IF(I < 0 .OR. I > IBP1) CYCLE
   IF(J < 0 .OR. J > JBP1) CYCLE
   IF(K < 0 .OR. K > KBP1) CYCLE
   IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE ! Cycle in case Cartesian cell inside OBSTS.
   IF (.NOT.ALLOCATED(MESHES(NM)%D_SOURCE)) CYCLE
   VOL=DX(I)*DY(J)*DZ(K); D_SOURCE(I,J,K) =0._EB; M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES)=0._EB
   DO JCC=1,CC%NCELL
      D_SOURCE(I,J,K) = D_SOURCE(I,J,K) + CUT_CELL(ICC)%D_SOURCE(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
      DO NS=1,N_TRACKED_SPECIES
        M_DOT_PPP(I,J,K,NS) = M_DOT_PPP(I,J,K,NS) + CUT_CELL(ICC)%M_DOT_PPP(NS,JCC)*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO
   D_SOURCE(I,J,K) = D_SOURCE(I,J,K)/VOL
   M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES) = M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES)/VOL
ENDDO

RETURN
END SUBROUTINE COPY_UNST_DM_TO_CART


! --------------------------- GET_H_GUARD_CUTCELL -----------------------------------

SUBROUTINE GET_H_GUARD_CUTCELL(IPZ,HP)

! Fill ghost cell H values for cut-cell regions at mesh boundaries.
! Handles: RC faces (cut-cell to regular), cut-faces at interpolated boundaries.
! assumes POINT_TO_MESH(NM) has been called.
INTEGER, INTENT(IN) :: IPZ
REAL(EB), INTENT(INOUT), POINTER, DIMENSION(:,:,:) :: HP

! Local Variables:
INTEGER :: IW, II, JJ, KK, IIG, JJG, KKG, ICC, ICC_GHOST, ICC_INT, NOM, IIO, JJO, KKO, IOR, ICC_EXT
INTEGER :: X1AXIS, ICF_EXT, ISHF(IAXIS:KAXIS)
REAL(EB) :: H_INT, D_INT, D_GHOST, D_EXT, X_FACE, SUM_FLUX, A_INT, A_EXT, H_EXT, H_GHOST
TYPE (WALL_TYPE), POINTER :: WC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (CC_CUTFACE_TYPE), POINTER :: CF



IF (PRES_FLAG==ULMAT_FLAG) THEN
   IF (ONE_UNKH_PER_CUTCELL) THEN
     ! To DO.
  ELSE
      WALL_CELL_LOOP_ULMAT : DO IW=1,N_EXTERNAL_WALL_CELLS
        WC => WALL(IW); BC => BOUNDARY_COORD(WC%BC_INDEX); ICC = CCVAR(BC%II,BC%JJ,BC%KK,CC_IDCC)
         IF (ZONE_MESH(PRESSURE_ZONE(BC%IIG,BC%JJG,BC%KKG))%CONNECTED_ZONE_PARENT/=IPZ .OR. ICC<1) CYCLE WALL_CELL_LOOP_ULMAT
         IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN
            CUT_CELL(ICC)%H(1:CUT_CELL(ICC)%NCELL) = HP(BC%II,BC%JJ,BC%KK)
         ELSE
            CUT_CELL(ICC)%HS(1:CUT_CELL(ICC)%NCELL) = HP(BC%II,BC%JJ,BC%KK)
        ENDIF
      ENDDO WALL_CELL_LOOP_ULMAT
  ENDIF

ELSE
   ! EXTERNAL_WALL_CELL approach for GLMAT with refinement support
   WALL_CELL_LOOP : DO IW=1,N_EXTERNAL_WALL_CELLS
      WC => WALL(IW); BC => BOUNDARY_COORD(WC%BC_INDEX)
      IF (ALL((/INTERPOLATED_BOUNDARY, PERIODIC_BOUNDARY/) /= WC%BOUNDARY_TYPE)) CYCLE
      EWC => EXTERNAL_WALL(IW)
      NOM = EWC%NOM; IF (NOM < 1) CYCLE
      OM => OMESH(NOM)

      II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK
      IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR; X1AXIS = ABS(IOR)
      ICC_GHOST = CCVAR(II,JJ,KK,CC_IDCC)      ! Ghost cell cut-cell index
      ICC_INT   = CCVAR(IIG,JJG,KKG,CC_IDCC)   ! Internal cell cut-cell index
      IF (ICC_GHOST < 1 .AND. ICC_INT < 1) CYCLE ! This is handled by EXTERNAL_WALL_LOOP in pres.f90

      ! Get H_INT (internal cell)
      IF (ICC_INT > 0) THEN
         IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN; H_INT = CUT_CELL(ICC_INT)%H(1)
         ELSE;                H_INT = CUT_CELL(ICC_INT)%HS(1)
         ENDIF
      ELSE
         H_INT = HP(IIG,JJG,KKG)
      ENDIF
      ! Get interface location, D_INT, D_GHOST, and A_INT based on face type
      IF (WC%CUT_FACE_INDEX > 0) THEN
         ! Cut-face at boundary
         CF => CUT_FACE(WC%CUT_FACE_INDEX)
         A_INT = SUM(CF%AREA(1:CF%NFACE))  ! Total cut-face area
      ELSE
         ! RC face (regular Cartesian face)
         SELECT CASE(X1AXIS)
         CASE(1); A_INT = DY(JJG) * DZ(KKG)
         CASE(2); A_INT = DX(IIG) * DZ(KKG)
         CASE(3); A_INT = DX(IIG) * DY(JJG)
         END SELECT
      ENDIF
      ! X_FACE location, D_INT and D_GHOST using cut-cell centroids when available
      SELECT CASE(X1AXIS)
      CASE(1)
         IF (IOR > 0) THEN; X_FACE = X(IIG-1); ELSE; X_FACE = X(IIG); ENDIF
         IF (ICC_INT > 0) THEN; D_INT = ABS(CUT_CELL(ICC_INT)%XYZCEN(IAXIS,1) - X_FACE)
         ELSE;                  D_INT = ABS(XC(IIG) - X_FACE)
         ENDIF
         IF (ICC_GHOST > 0) THEN; D_GHOST = ABS(CUT_CELL(ICC_GHOST)%XYZCEN(IAXIS,1) - X_FACE)
         ELSE;                    D_GHOST = ABS(XC(II) - X_FACE)
         ENDIF
      CASE(2)
         IF (IOR > 0) THEN; X_FACE = Y(JJG-1); ELSE; X_FACE = Y(JJG); ENDIF
         IF (ICC_INT > 0) THEN; D_INT = ABS(CUT_CELL(ICC_INT)%XYZCEN(JAXIS,1) - X_FACE)
         ELSE;                  D_INT = ABS(YC(JJG) - X_FACE)
         ENDIF
         IF (ICC_GHOST > 0) THEN; D_GHOST = ABS(CUT_CELL(ICC_GHOST)%XYZCEN(JAXIS,1) - X_FACE)
         ELSE;                    D_GHOST = ABS(YC(JJ) - X_FACE)
         ENDIF
      CASE(3)
         IF (IOR > 0) THEN; X_FACE = Z(KKG-1); ELSE; X_FACE = Z(KKG); ENDIF
         IF (ICC_INT > 0) THEN; D_INT = ABS(CUT_CELL(ICC_INT)%XYZCEN(KAXIS,1) - X_FACE)
         ELSE;                  D_INT = ABS(ZC(KKG) - X_FACE)
         ENDIF
         IF (ICC_GHOST > 0) THEN; D_GHOST = ABS(CUT_CELL(ICC_GHOST)%XYZCEN(KAXIS,1) - X_FACE)
         ELSE;                    D_GHOST = ABS(ZC(KK) - X_FACE)
         ENDIF
      END SELECT

      ! Flux-matched sum over external cells
      SUM_FLUX = 0._EB
      DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
         DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
            DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
               ICC_EXT = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
               IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN; H_EXT = OM%H(IIO,JJO,KKO)
               ELSE;                H_EXT = OM%HS(IIO,JJO,KKO)
               ENDIF
               SELECT CASE(X1AXIS)
               CASE(1)
                  A_EXT = MESHES(NOM)%DY(JJO) * MESHES(NOM)%DZ(KKO)
                  IF (ICC_EXT > 0) THEN; D_EXT = ABS(MESHES(NOM)%CUT_CELL(ICC_EXT)%XYZCEN(IAXIS,1) - X_FACE)
                  ELSE;                  D_EXT = ABS(MESHES(NOM)%XC(IIO) - X_FACE)
                  ENDIF
               CASE(2)
                  A_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DZ(KKO)
                  IF (ICC_EXT > 0) THEN; D_EXT = ABS(MESHES(NOM)%CUT_CELL(ICC_EXT)%XYZCEN(JAXIS,1) - X_FACE)
                  ELSE;                  D_EXT = ABS(MESHES(NOM)%YC(JJO) - X_FACE)
                  ENDIF
               CASE(3)
                  A_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DY(JJO)
                  IF (ICC_EXT > 0) THEN; D_EXT = ABS(MESHES(NOM)%CUT_CELL(ICC_EXT)%XYZCEN(KAXIS,1) - X_FACE)
                  ELSE;                  D_EXT = ABS(MESHES(NOM)%ZC(KKO) - X_FACE)
                  ENDIF
               END SELECT
               ISHF(IAXIS:KAXIS) = 0; IF(IOR < 0) ISHF(X1AXIS) = -1
               ICF_EXT = MESHES(NOM)%FCVAR(IIO+ISHF(IAXIS),JJO+ISHF(JAXIS),KKO+ISHF(KAXIS),CC_IDCF,X1AXIS)
               IF (ICF_EXT > 0) A_EXT = SUM(MESHES(NOM)%CUT_FACE(ICF_EXT)%AREA(1:MESHES(NOM)%CUT_FACE(ICF_EXT)%NFACE))
               IF (EWC%AREA_RATIO < 0.99_EB) A_EXT = A_INT
               SUM_FLUX = SUM_FLUX + (H_EXT - H_INT) / (D_EXT + D_INT) * A_EXT
            ENDDO
         ENDDO
      ENDDO

      ! Compute ghost cell value
      H_GHOST = H_INT + (D_INT + D_GHOST) / A_INT * SUM_FLUX

      ! Fill ghost cells
      HP(II,JJ,KK) = H_GHOST
      IF (ICC_GHOST > 0) THEN
         IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN; CUT_CELL(ICC_GHOST)%H(1:CUT_CELL(ICC_GHOST)%NCELL)  = H_GHOST
         ELSE;                CUT_CELL(ICC_GHOST)%HS(1:CUT_CELL(ICC_GHOST)%NCELL) = H_GHOST
         ENDIF
      ENDIF
   ENDDO WALL_CELL_LOOP
ENDIF

END SUBROUTINE GET_H_GUARD_CUTCELL





! ---------------------- UNSTRUCTURED_POISSON_RESIDUAL_RC ----------------------

SUBROUTINE UNSTRUCTURED_POISSON_RESIDUAL_RC(I,J,K,HP,RHOP,P,RES,DO_SEPARABLE)

! NOTE: Assumes POINT_TO_MESH(NM) has been called.

INTEGER, INTENT(IN) :: I,J,K
REAL(EB),INTENT(OUT):: RES
LOGICAL, INTENT(IN) :: DO_SEPARABLE
REAL(EB), POINTER, INTENT(IN), DIMENSION(:,:,:) :: HP,RHOP,P

! Local Dummy vars:
REAL(EB):: RHSS, LHSS

RES=0._EB; IF(CCVAR(I,J,K,CC_UNKZ)<=0) RETURN

CALL GET_RHSLHS_POISSON_RC(I,J,K,HP,RHOP,P,RHSS,LHSS,DO_SEPARABLE)
RES  = ABS(RHSS-LHSS)

RETURN
END SUBROUTINE UNSTRUCTURED_POISSON_RESIDUAL_RC

! --------------------------- GET_RHSLHS_POISSON_RC ----------------------------

SUBROUTINE GET_RHSLHS_POISSON_RC(I,J,K,HP,RHOP,P,RHSS,LHSS,DO_SEPARABLE)

! NOTE: Assumes POINT_TO_MESH(NM) has been called.

INTEGER, INTENT(IN) :: I,J,K
REAL(EB),INTENT(OUT):: RHSS,LHSS
LOGICAL, INTENT(IN) :: DO_SEPARABLE
REAL(EB), POINTER, INTENT(IN), DIMENSION(:,:,:) :: HP,RHOP,P

! Local Dummy vars:
INTEGER :: X1AXIS,IRC,ILH,FCT,ICC,JCC
REAL(EB):: PRFCT, AF, VOL, DIV_FN_VOL, RDN, H1, H2, P1, RHO1, KR1, P2, RHO2, KR2, RHOF, FBC, CCM1, CCP1, FN_B

VOL = DX(I)*DY(J)*DZ(K)
DIV_FN_VOL = 0._EB; LHSS = 0._EB
DO_SEPARABLE_IF : IF (DO_SEPARABLE) THEN

   ! X axis:
   X1AXIS = IAXIS; AF = DY(J)*DZ(K)
   DO ILH=-1,0
      FCT  = 2*ILH+1
      RDN  = RDXN(I+ILH)
      IRC  = FCVAR(I+ILH,J,K,CC_IDRC,X1AXIS)
      IF (IRC>0) RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND) - RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
      H1   = HP(I+ILH,J,K); H2 = HP(I+ILH+1,J,K)
      FBC = 1._EB; IF (WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* FVX(I+ILH,J,K) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*(H2-H1)*RDN * AF
   ENDDO
   ! Y axis:
   X1AXIS = JAXIS; AF = DX(I)*DZ(K)
   DO ILH=-1,0
      FCT  = 2*ILH+1
      RDN  = RDYN(J+ILH)
      IRC  = FCVAR(I,J+ILH,K,CC_IDRC,X1AXIS)
      IF (IRC>0) RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND) - RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
      H1   = HP(I,J+ILH,K); H2 = HP(I,J+ILH+1,K)
      FBC = 1._EB; IF (WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* FVY(I,J+ILH,K) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*(H2-H1)*RDN * AF
   ENDDO
   ! Z axis:
   X1AXIS = KAXIS; AF = DX(I)*DY(J)
   DO ILH=-1,0
      FCT  = 2*ILH+1
      RDN  = RDZN(K+ILH)
      IRC  = FCVAR(I,J,K+ILH,CC_IDRC,X1AXIS)
      IF (IRC>0) RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND) - RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
      H1   = HP(I,J,K+ILH); H2 = HP(I,J,K+ILH+1)
      FBC = 1._EB; IF (WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* FVZ(I,J,K+ILH) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*(H2-H1)*RDN * AF
   ENDDO

ELSE DO_SEPARABLE_IF

   PRFCT=0._EB; IF(MESHES(LOWER_MESH_INDEX)%PREDICTOR) PRFCT=1._EB

   ! X axis:
   X1AXIS = IAXIS; AF = DY(J)*DZ(K)
   DO ILH=-1,0
      FCT  = 2*ILH+1; FBC = 1._EB
      RDN  = RDXN(I+ILH); CCM1=0.5_EB; CCP1=0.5_EB

      P1   = P(   I+ILH,J,K); P2   = P(   I+ILH+1,J,K)
      KR1  = KRES(I+ILH,J,K); KR2  = KRES(I+ILH+1,J,K)
      RHO1 = RHOP(I+ILH,J,K); RHO2 = RHOP(I+ILH+1,J,K)

      IF (FCVAR(I+ILH,J,K,CC_UNKF,IAXIS)>0) THEN
         FN_B = F_LINK(FCVAR(I+ILH,J,K,CC_UNKF,IAXIS))
      ELSE
         FN_B = -(P(I+ILH,J,K)*RHOP(I+ILH+1,J,K)+P(I+ILH+1,J,K)*RHOP(I+ILH,J,K))*&
                 (1._EB/RHOP(I+ILH+1,J,K)-1._EB/RHOP(I+ILH,J,K))*RDXN(I+ILH)/(RHOP(I+ILH+1,J,K)+RHOP(I+ILH,J,K))
      ENDIF
      IRC  = FCVAR(I+ILH,J,K,CC_IDRC,X1AXIS)
      IF (IRC>0) THEN
         IF(RC_FACE(IRC)%IWC>0 .AND. &
            ANY(WALL(RC_FACE(IRC)%IWC)%BOUNDARY_TYPE==(/NULL_BOUNDARY,MIRROR_BOUNDARY,SOLID_BOUNDARY/))) THEN
            FN_B = 0._EB
         ELSE
            RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
            CCM1= RDN*(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-X(I+ILH))
            CCP1= RDN*(X(I+ILH)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND) )
            IF(RC_FACE(IRC)%CELL_LIST(1,LOW_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,LOW_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,LOW_IND)
               RHO1 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%CELL_LIST(1,HIGH_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,HIGH_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,HIGH_IND)
               RHO2 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%UNKF>0) THEN
               FN_B = F_LINK(RC_FACE(IRC)%UNKF)
            ELSE
               CALL GET_RCFACE_BAROCLINIC_TORQUE(PRFCT,IRC,CCM1,CCP1,FN_B,HP,RHOP)
            ENDIF
         ENDIF
      ENDIF
      RHOF = CCM1*RHO1 + CCP1*RHO2

      IF ( (I+ILH>0 .AND. I+ILH<IBAR) .AND. &
            WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* ( FVX(I+ILH,J,K) - FBC*FVX_B(I+ILH,J,K) ) * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( 1._EB/RHOF * (P2-P1) + (KR2-KR1) )*RDN * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN + FVX_B(I+ILH,J,K) ) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN + FN_B ) * AF
   ENDDO
   ! Y axis:
   X1AXIS = JAXIS; AF = DX(I)*DZ(K)
   DO ILH=-1,0
      FCT  = 2*ILH+1; FBC = 1._EB
      RDN  = RDYN(J+ILH); CCM1=0.5_EB; CCP1=0.5_EB

      P1   = P(   I,J+ILH,K); P2   = P(   I,J+ILH+1,K)
      KR1  = KRES(I,J+ILH,K); KR2  = KRES(I,J+ILH+1,K)
      RHO1 = RHOP(I,J+ILH,K); RHO2 = RHOP(I,J+ILH+1,K)

      IF (FCVAR(I,J+ILH,K,CC_UNKF,JAXIS)>0) THEN
         FN_B = F_LINK(FCVAR(I,J+ILH,K,CC_UNKF,JAXIS))
      ELSE
         FN_B = -(P(I,J+ILH,K)*RHOP(I,J+ILH+1,K)+P(I,J+ILH+1,K)*RHOP(I,J+ILH,K))*&
                 (1._EB/RHOP(I,J+ILH+1,K)-1._EB/RHOP(I,J+ILH,K))*RDYN(J+ILH)/ (RHOP(I,J+ILH+1,K)+RHOP(I,J+ILH,K))
      ENDIF
      IRC  = FCVAR(I,J+ILH,K,CC_IDRC,X1AXIS)
      IF (IRC>0) THEN
         IF(RC_FACE(IRC)%IWC>0 .AND. &
            ANY(WALL(RC_FACE(IRC)%IWC)%BOUNDARY_TYPE==(/NULL_BOUNDARY,MIRROR_BOUNDARY,SOLID_BOUNDARY/))) THEN
            FN_B = 0._EB
         ELSE
            RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
            CCM1= RDN*(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-Y(J+ILH))
            CCP1= RDN*(Y(J+ILH)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND) )
            IF(RC_FACE(IRC)%CELL_LIST(1,LOW_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,LOW_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,LOW_IND)
               RHO1 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%CELL_LIST(1,HIGH_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,HIGH_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,HIGH_IND)
               RHO2 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%UNKF>0) THEN
               FN_B = F_LINK(RC_FACE(IRC)%UNKF)
            ELSE
               CALL GET_RCFACE_BAROCLINIC_TORQUE(PRFCT,IRC,CCM1,CCP1,FN_B,HP,RHOP)
            ENDIF
         ENDIF
      ENDIF
      RHOF = CCM1*RHO1 + CCP1*RHO2

      IF ( (J+ILH>0 .AND. J+ILH<JBAR) .AND. &
            WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* ( FVY(I,J+ILH,K) - FBC*FVY_B(I,J+ILH,K) ) * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( 1._EB/RHOF * (P2-P1) + (KR2-KR1) )*RDN * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN +FVY_B(I,J+ILH,K) ) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN + FN_B ) * AF
   ENDDO
   ! Z axis:
   X1AXIS = KAXIS; AF = DX(I)*DY(J)
   DO ILH=-1,0
      FCT  = 2*ILH+1; FBC = 1._EB
      RDN  = RDZN(K+ILH); CCM1=0.5_EB; CCP1=0.5_EB

      P1   = P(   I,J,K+ILH); P2   = P(   I,J,K+ILH+1)
      KR1  = KRES(I,J,K+ILH); KR2  = KRES(I,J,K+ILH+1)
      RHO1 = RHOP(I,J,K+ILH); RHO2 = RHOP(I,J,K+ILH+1)

      IF (FCVAR(I,J,K+ILH,CC_UNKF,KAXIS)>0) THEN
         FN_B = F_LINK(FCVAR(I,J,K+ILH,CC_UNKF,KAXIS))
      ELSE
         FN_B = -(P(I,J,K+ILH)*RHOP(I,J,K+ILH+1)+P(I,J,K+ILH+1)*RHOP(I,J,K+ILH))*&
                 (1._EB/RHOP(I,J,K+ILH+1)-1._EB/RHOP(I,J,K+ILH))*RDZN(K+ILH)/ (RHOP(I,J,K+ILH+1)+RHOP(I,J,K+ILH))
      ENDIF
      IRC  = FCVAR(I,J,K+ILH,CC_IDRC,X1AXIS)
      IF (IRC>0) THEN
         IF(RC_FACE(IRC)%IWC>0 .AND. &
            ANY(WALL(RC_FACE(IRC)%IWC)%BOUNDARY_TYPE==(/NULL_BOUNDARY,MIRROR_BOUNDARY,SOLID_BOUNDARY/))) THEN
            FN_B = 0._EB
         ELSE
            RDN = 1._EB/(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND))
            CCM1= RDN*(RC_FACE(IRC)%XCEN(X1AXIS,HIGH_IND)-Z(K+ILH))
            CCP1= RDN*(Z(K+ILH)-RC_FACE(IRC)%XCEN(X1AXIS,LOW_IND) )
            IF(RC_FACE(IRC)%CELL_LIST(1,LOW_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,LOW_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,LOW_IND)
               RHO1 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%CELL_LIST(1,HIGH_IND)==CC_FTYPE_CFGAS) THEN
               ICC  = RC_FACE(IRC)%CELL_LIST(2,HIGH_IND)
               JCC  = RC_FACE(IRC)%CELL_LIST(3,HIGH_IND)
               RHO2 = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ENDIF
            IF(RC_FACE(IRC)%UNKF>0) THEN
               FN_B = F_LINK(RC_FACE(IRC)%UNKF)
            ELSE
               CALL GET_RCFACE_BAROCLINIC_TORQUE(PRFCT,IRC,CCM1,CCP1,FN_B,HP,RHOP)
            ENDIF
         ENDIF
      ENDIF
      RHOF = CCM1*RHO1 + CCP1*RHO2

      IF ( (K+ILH>1 .AND. K+ILH<KBAR) .AND.  &
            WALL(CELL(CELL_INDEX(I,J,K))%WALL_INDEX(FCT*X1AXIS))%BOUNDARY_TYPE==SOLID_BOUNDARY) FBC = 0._EB
      DIV_FN_VOL = DIV_FN_VOL + REAL(FCT,EB)* ( FVZ(I,J,K+ILH) - FBC*FVZ_B(I,J,K+ILH) ) * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( 1._EB/RHOF * (P2-P1) + (KR2-KR1) )*RDN * AF
      !LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN + FVZ_B(I,J,K+ILH) ) * AF
      LHSS       = LHSS       + REAL(FCT,EB)* FBC*( (P2/RHO2-P1/RHO1+KR2-KR1)*RDN + FN_B ) * AF
   ENDDO

ENDIF DO_SEPARABLE_IF

LHSS = LHSS/VOL
RHSS = -(DDDT(I,J,K)*VOL + DIV_FN_VOL)/VOL ! Normalize to cartesian cell volume.

RETURN
END SUBROUTINE GET_RHSLHS_POISSON_RC


! ------------------------UNSTRUCTURED_POISSON_RESIDUAL ------------------------

SUBROUTINE UNSTRUCTURED_POISSON_RESIDUAL(NM,I,J,K,HP,RHOP,P, &
                                        RES,DO_SEPARABLE)

INTEGER, INTENT(IN) :: NM

INTEGER, INTENT(IN) :: I,J,K
REAL(EB),INTENT(OUT):: RES
LOGICAL, INTENT(IN) :: DO_SEPARABLE
REAL(EB), POINTER, INTENT(IN), DIMENSION(:,:,:) :: HP,RHOP,P

! Local Variables:
INTEGER :: ICC,JCC
REAL(EB):: RHSS,LHSS,DIV_FN,DIV_FN_VOL,LHSS_CC,PRFCT

PRFCT=0._EB
IF(MESHES(NM)%PREDICTOR) PRFCT=1._EB
RES=0._EB; LHSS=0._EB; DIV_FN_VOL = 0._EB
! CUT_CELL entry:
ICC=CCVAR(I,J,K,CC_IDCC)
DO_SEPARABLE_IF : IF (DO_SEPARABLE) THEN
   IF(ONE_UNKH_PER_CUTCELL) THEN
      DO JCC=1,CUT_CELL(ICC)%NCELL
        ! Here we add div(F) in the cut-cell and DDDT:
        CALL GET_FN_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC,DIV_FN, &
        SUBSTRACT_BAROCLINIC=.FALSE.)

        ! Compute int(Grad H dot n) ds for cut-cell:
        CALL GET_LHS_CUTCELL(PRFCT,ICC,JCC,LHSS_CC,I,J,K,HP,RHOP,P,DO_SEPARABLE=.TRUE.)
        LHSS       =  LHSS_CC
        ! Add to RHSS:
        RHSS = -(CUT_CELL(ICC)%DDDTVOL(JCC) + DIV_FN*CUT_CELL(ICC)%VOLUME(JCC))
        RES  = RES + ABS(RHSS-LHSS)
      ENDDO
   ELSE
      DO JCC=1,CUT_CELL(ICC)%NCELL
        ! Here we add div(F) in the cut-cell and DDDT:
        CALL GET_FN_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC,DIV_FN, &
        SUBSTRACT_BAROCLINIC=.FALSE.)
        DIV_FN_VOL = DIV_FN_VOL + DIV_FN*CUT_CELL(ICC)%VOLUME(JCC)

        ! Compute int(Grad H dot n) ds for cut-cell:
        CALL GET_LHS_CUTCELL(PRFCT,ICC,JCC,LHSS_CC,I,J,K,HP,RHOP,P,DO_SEPARABLE=.TRUE.)
        LHSS       = LHSS       + LHSS_CC
      ENDDO
      ! Add to RHSS:
      RHSS = -(CUT_CELL(ICC)%DDDTVOL(1) + DIV_FN_VOL)
      RES  = ABS(RHSS-LHSS)
   ENDIF
ELSE DO_SEPARABLE_IF
   IF(SUM(CUT_CELL(ICC)%VOLUME(1:CUT_CELL(ICC)%NCELL))/(DX(I)*DY(J)*DZ(K)) < V_THRESH_INSPRES) RETURN
   IF(ONE_UNKH_PER_CUTCELL) THEN
      DO JCC=1,CUT_CELL(ICC)%NCELL
        ! Here we add div(F) in the cut-cell and DDDT:
        CALL GET_FN_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC,DIV_FN, &
        SUBSTRACT_BAROCLINIC=.TRUE.)

        ! Compute int(Grad H dot n) ds for cut-cell:
        CALL GET_LHS_CUTCELL(PRFCT,ICC,JCC,LHSS_CC,I,J,K,HP,RHOP,P,DO_SEPARABLE=.FALSE.)
        LHSS       =  LHSS_CC
        ! Add to RHSS:
        RHSS = -(CUT_CELL(ICC)%DDDTVOL(JCC) + DIV_FN*CUT_CELL(ICC)%VOLUME(JCC))
        RES  = RES + ABS(RHSS-LHSS)
      ENDDO
   ELSE
      DO JCC=1,CUT_CELL(ICC)%NCELL
        ! Here we add div(F) in the cut-cell and DDDT:
        CALL GET_FN_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC,DIV_FN, &
        SUBSTRACT_BAROCLINIC=.TRUE.)
        DIV_FN_VOL = DIV_FN_VOL + DIV_FN*CUT_CELL(ICC)%VOLUME(JCC)

        ! Compute (1/rho*Grad(p)-Grad(Kres)):
        CALL GET_LHS_CUTCELL(PRFCT,ICC,JCC,LHSS_CC,I,J,K,HP,RHOP,P,DO_SEPARABLE=.FALSE.)
        LHSS       = LHSS       + LHSS_CC
      ENDDO
      ! Add to RHSS:
      RHSS = -(CUT_CELL(ICC)%DDDTVOL(1) + DIV_FN_VOL)
      RES  = ABS(RHSS-LHSS)
   ENDIF

ENDIF DO_SEPARABLE_IF

RES  = RES/(DX(I)*DY(J)*DZ(K)) ! Normalize to cartesian cell volume.

! IF (DO_SEPARABLE) THEN
! !IF(I==38 .AND. J==10 .AND. K==2) THEN
! IF(RES>1.E-6_EB) &
! WRITE(LU_ERR,*) 'Cut-cell I,J,K,RES=',I,J,K,RHSS,LHSS,RES,SUM(CUT_CELL(ICC)%VOLUME(1:CUT_CELL(ICC)%NCELL)) !,':',&
! !MESHES(NM)%ZONE_MESH(0)%F_H(CUT_CELL(ICC)%UNKH(1))
! !ENDIF
! ENDIF

RETURN
END SUBROUTINE UNSTRUCTURED_POISSON_RESIDUAL


! --------------------------------- GET_LHS_CUTCELL --------------------------------


SUBROUTINE GET_LHS_CUTCELL(PRFCT,ICC,JCC,LHSS_CC,I,J,K,HP,RHOP,P,DO_SEPARABLE)

! NOTE: Assumes POINT_TO_MESH(NM) has been called.

REAL(EB),INTENT(IN) :: PRFCT
INTEGER, INTENT(IN) :: ICC,JCC,I,J,K
REAL(EB),INTENT(OUT):: LHSS_CC
LOGICAL, INTENT(IN) :: DO_SEPARABLE
REAL(EB), POINTER, INTENT(IN), DIMENSION(:,:,:) :: HP,RHOP,P


! Local Variables:
INTEGER :: IFC,IFACE,X1AXIS,LOWHIGH,ILH,ICFA,IFC2,IFACE2,ICLO,ICHI,JCLO,JCHI,IRC
REAL(EB):: HC1,HC2,X1,X2,XFC,AF,FN,FCT,P1,P2,KR1,KR2,RHO1,RHO2,RHOF,FN_B,CCM1,CCP1

LHSS_CC = 0._EB

DO_SEPARABLE_IF : IF (DO_SEPARABLE) THEN

   ! Compute int(Grad H dot n) ds for cut-cell:
   IFC_LOOP_1 : DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
      IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
      AF   = 0._EB
      FN   = 0._EB
      SELECT CASE(CUT_CELL(ICC)%FACE_LIST(1,IFACE))
      CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
         LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
         X1AXIS  = CUT_CELL(ICC)%FACE_LIST(3,IFACE)
         ILH     =        LOWHIGH - 1
         FCT     = REAL(2*LOWHIGH - 3, EB)
         HC1  = HP(I,J,K); HC2 = HC1
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            AF   = DY(J)*DZ(K)
            IF(LOWHIGH==HIGH_IND) THEN
               HC1  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               HC2  = HP(I+1,J,K)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = XC(I+1); IRC  = FCVAR(I,J,K,CC_IDRC,X1AXIS)
               IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
            ELSE
               HC1  = HP(I-1,J,K)
               HC2  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               X1   = XC(I-1); IRC  = FCVAR(I-1,J,K,CC_IDRC,X1AXIS)
               IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
            ENDIF
         CASE(JAXIS)
            AF   = DX(I)*DZ(K)
            IF(LOWHIGH==HIGH_IND) THEN
               HC1  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               HC2  = HP(I,J+1,K)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = YC(J+1); IRC  = FCVAR(I,J,K,CC_IDRC,X1AXIS)
               IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
            ELSE
               HC1  = HP(I,J-1,K)
               HC2  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               X1   = YC(J-1); IRC  = FCVAR(I,J-1,K,CC_IDRC,X1AXIS)
               IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
            ENDIF
         CASE(KAXIS)
            AF   = DX(I)*DY(J)
            IF(LOWHIGH==HIGH_IND) THEN
               HC1  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               HC2  = HP(I,J,K+1)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = ZC(K+1); IRC  = FCVAR(I,J,K,CC_IDRC,X1AXIS)
               IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
            ELSE
               HC1  = HP(I,J,K-1)
               HC2  = PRFCT*CUT_CELL(ICC)%H(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%HS(JCC)
               X1   = ZC(K-1); IRC  = FCVAR(I,J,K-1,CC_IDRC,X1AXIS)
               IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
            ENDIF
         END SELECT
         IF(IRC>0) THEN
            IF( RC_FACE(IRC)%IWC>0 .AND. &
                ANY(WALL(RC_FACE(IRC)%IWC)%BOUNDARY_TYPE==(/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/)) ) &
                CYCLE IFC_LOOP_1
         ENDIF
         FN   = FCT*(HC2-HC1)/(X2-X1) ! Grad H dot n
      CASE(CC_FTYPE_CFGAS) ! GASPHASE CUT FACE:
         LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
         FCT     = REAL(2*LOWHIGH - 3, EB)
         IFC2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
         IF(CUT_FACE(IFC2)%IWC>0 .AND. &
            ANY(WALL(CUT_FACE(IFC2)%IWC)%BOUNDARY_TYPE==(/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) &
            CYCLE IFC_LOOP_1
         IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
         AF      = CUT_FACE(IFC2)%AREA(IFACE2)
         ! Low side H:
         ICLO = CUT_FACE(IFC2)%CELL_LIST(2,LOW_IND,IFACE2); JCLO = CUT_FACE(IFC2)%CELL_LIST(3,LOW_IND,IFACE2);
         !HC1  = HP(CUT_CELL(ICLO)%IJK(IAXIS),CUT_CELL(ICLO)%IJK(JAXIS),CUT_CELL(ICLO)%IJK(KAXIS))
         HC1  = PRFCT*CUT_CELL(ICLO)%H(JCLO)+(1._EB-PRFCT)*CUT_CELL(ICLO)%HS(JCLO)
         ! High side H:
         ICHI = CUT_FACE(IFC2)%CELL_LIST(2,HIGH_IND,IFACE2); JCHI = CUT_FACE(IFC2)%CELL_LIST(3,HIGH_IND,IFACE2);
         !HC2  = HP(CUT_CELL(ICHI)%IJK(IAXIS),CUT_CELL(ICHI)%IJK(JAXIS),CUT_CELL(ICHI)%IJK(KAXIS))
         HC2  = PRFCT*CUT_CELL(ICHI)%H(JCHI) + (1._EB-PRFCT)*CUT_CELL(ICHI)%HS(JCHI)

         X1AXIS = CUT_FACE(IFC2)%IJK(KAXIS+1)
         X1   = CUT_FACE(IFC2)%XCENLOW(X1AXIS,IFACE2)
         X2   = CUT_FACE(IFC2)%XCENHIGH(X1AXIS,IFACE2)
         FN   = FCT*(HC2-HC1)/(X2-X1) ! Grad H dot n
      CASE(CC_FTYPE_CFINB) ! INBOUNDARY CUT FACE
         FCT     = 1._EB    ! Normal vector defined into the body.
         IFC2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
         IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
         ICFA    = CUT_FACE(IFC2)%CFACE_INDEX(IFACE2)
         AF      = CUT_FACE(IFC2)%AREA(IFACE2)
         ! FN      = 0._EB
      END SELECT
      LHSS_CC = LHSS_CC + AF*FN
   ENDDO IFC_LOOP_1

ELSE DO_SEPARABLE_IF

   ! Compute (1/rho*Grad(p)+Grad(Kres)):
   CCM1 = 0.5_EB; CCP1 = 0.5_EB
   IFC_LOOP_2 : DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
      IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
      AF    = 0._EB; FN = 0._EB; FN_B = 0._EB
      SELECT CASE(CUT_CELL(ICC)%FACE_LIST(1,IFACE))
      CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
         LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
         X1AXIS  = CUT_CELL(ICC)%FACE_LIST(3,IFACE)
         ILH     =        LOWHIGH - 1
         FCT     = REAL(2*LOWHIGH - 3, EB)
         P1   = P(I,J,K); P2 = P1
         KR1  = KRES(I,J,K); KR2 = KR1
         RHO1 = PRFCT*CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC); RHO2=RHO1
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            AF   = DY(J)*DZ(K); XFC  = X(I-1+ILH)   ! Face location in face normal direction.
            !FN_B = FVX_B(I-1+ILH,J,K)
            IRC  = FCVAR(I-1+ILH,J,K,CC_IDRC,X1AXIS)
            IF(LOWHIGH==HIGH_IND) THEN
               P2   = P(I+1,J,K); KR2 = KRES(I+1,J,K)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = XC(I+1); IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               RHO2 = RHOP(I+1,J,K)
            ELSE
               P1   = P(I-1,J,K); KR1 = KRES(I-1,J,K)
               X1   = XC(I-1); IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               RHO1 = RHOP(I-1,J,K)
            ENDIF
         CASE(JAXIS)
            AF   = DX(I)*DZ(K); XFC  = Y(J-1+ILH)   ! Face location in face normal direction.
            !FN_B = FVY_B(I,J-1+ILH,K)
            IRC  = FCVAR(I,J-1+ILH,K,CC_IDRC,X1AXIS)
            IF(LOWHIGH==HIGH_IND) THEN
               P2   = P(I,J+1,K); KR2 = KRES(I,J+1,K)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = YC(J+1); IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               RHO2 = RHOP(I,J+1,K)
            ELSE
               P1   = P(I,J-1,K); KR1 = KRES(I,J-1,K)
               X1   = YC(J-1); IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               RHO1 = RHOP(I,J-1,K)
            ENDIF
         CASE(KAXIS)
            AF   = DX(I)*DY(J); XFC  = Z(K-1+ILH)   ! Face location in face normal direction.
            !FN_B = FVZ_B(I,J,K-1+ILH)
            IRC  = FCVAR(I,J,K-1+ILH,CC_IDRC,X1AXIS)
            IF(LOWHIGH==HIGH_IND) THEN
               P2   = P(I,J,K+1); KR2 = KRES(I,J,K+1)
               X1   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               X2   = ZC(K+1); IF(IRC>0) X2=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               RHO2 = RHOP(I,J,K+1)
            ELSE
               P1   = P(I,J,K-1); KR1 = KRES(I,J,K-1)
               X1   = ZC(K-1); IF(IRC>0) X1=RC_FACE(IRC)%XCEN(X1AXIS,LOWHIGH)
               X2   = CUT_CELL(ICC)%XYZCEN(X1AXIS,JCC)
               RHO1 = RHOP(I,J,K-1)
            ENDIF
         END SELECT
         IF(IRC>0) THEN
            IF( RC_FACE(IRC)%IWC>0 .AND. &
                ANY(WALL(RC_FACE(IRC)%IWC)%BOUNDARY_TYPE==(/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/)) ) &
                CYCLE IFC_LOOP_2
         ENDIF
         IF(RC_FACE(IRC)%UNKF>0) THEN
            FN_B = F_LINK(RC_FACE(IRC)%UNKF)
         ELSE
            CALL GET_RCFACE_BAROCLINIC_TORQUE(PRFCT,IRC,CCM1,CCP1,FN_B,HP,RHOP)
         ENDIF
         CCM1 = (X2-XFC)/(X2-X1); CCP1 = (XFC-X1)/(X2-X1); RHOF = CCM1*RHO1 + CCP1*RHO2
         !FN   = FCT * ( 1._EB/RHOF * (P2-P1) + (KR2-KR1) )/(X2-X1) ! (1/rho*Grad(p)+Grad(Kres))
         FN   = FCT * ( (P2/RHO2-P1/RHO1+KR2-KR1)/(X2-X1) + FN_B )  ! (Grad(p/rho)-p*Grad(1/rho)+Grad(Kres))
      CASE(CC_FTYPE_CFGAS) ! GASPHASE CUT FACE:
         LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
         FCT     = REAL(2*LOWHIGH - 3, EB)
         IFC2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
         IF(CUT_FACE(IFC2)%IWC>0 .AND. &
            ANY(WALL(CUT_FACE(IFC2)%IWC)%BOUNDARY_TYPE==(/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) &
            CYCLE IFC_LOOP_2
         IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
         AF      = CUT_FACE(IFC2)%AREA(IFACE2)
         ! Low side H:
         ICLO = CUT_FACE(IFC2)%CELL_LIST(2,LOW_IND,IFACE2); JCLO = CUT_FACE(IFC2)%CELL_LIST(3,LOW_IND,IFACE2);
         P1  = P(CUT_CELL(ICLO)%IJK(IAXIS),CUT_CELL(ICLO)%IJK(JAXIS),CUT_CELL(ICLO)%IJK(KAXIS))
         KR1 = KRES(CUT_CELL(ICLO)%IJK(IAXIS),CUT_CELL(ICLO)%IJK(JAXIS),CUT_CELL(ICLO)%IJK(KAXIS))
         RHO1= PRFCT*CUT_CELL(ICLO)%RHO(JCLO)+(1._EB-PRFCT)*CUT_CELL(ICLO)%RHOS(JCLO)
         ! High side H:
         ICHI = CUT_FACE(IFC2)%CELL_LIST(2,HIGH_IND,IFACE2); JCHI = CUT_FACE(IFC2)%CELL_LIST(3,HIGH_IND,IFACE2);
         P2  = P(CUT_CELL(ICHI)%IJK(IAXIS),CUT_CELL(ICHI)%IJK(JAXIS),CUT_CELL(ICHI)%IJK(KAXIS))
         KR2 = KRES(CUT_CELL(ICHI)%IJK(IAXIS),CUT_CELL(ICHI)%IJK(JAXIS),CUT_CELL(ICHI)%IJK(KAXIS))
         RHO2= PRFCT*CUT_CELL(ICHI)%RHO(JCHI)+(1._EB-PRFCT)*CUT_CELL(ICHI)%RHOS(JCHI)

         X1AXIS = CUT_FACE(IFC2)%IJK(KAXIS+1)
         XFC    = CUT_FACE(IFC2)%XYZCEN(X1AXIS,IFACE2)
         X1   = CUT_FACE(IFC2)%XCENLOW( X1AXIS,IFACE2)
         X2   = CUT_FACE(IFC2)%XCENHIGH(X1AXIS,IFACE2)
         RHOF = (X2-XFC)/(X2-X1)*RHO1 + (XFC-X1)/(X2-X1)*RHO2 ! CCM1*RHO1 + CCM2*RHO2
         !FN   = FCT * ( 1._EB/RHOF * (P2-P1) + (KR2-KR1) )/(X2-X1) ! (1/rho*Grad(p)+Grad(Kres))
         !FN   = FCT*((P2/RHO2-P1/RHO1+KR2-KR1)/(X2-X1)+CUT_FACE(IFC2)%FN_B(IFACE2)) ! (Grad(p/rho) -
                                                                                     !  p*Grad(1/rho)+Grad(Kres))
         IF(CUT_FACE(IFC2)%UNKF(IFACE2)>0) FN_B = F_LINK(CUT_FACE(IFC2)%UNKF(IFACE2))
         FN = FCT*( (P2/RHO2-P1/RHO1+KR2-KR1)/(X2-X1) + FN_B )
      CASE(CC_FTYPE_CFINB) ! INBOUNDARY CUT FACE
         FCT     = 1._EB    ! Normal vector defined into the body.
         IFC2    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
         IFACE2  = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
         ICFA    = CUT_FACE(IFC2)%CFACE_INDEX(IFACE2)
         AF      = CUT_FACE(IFC2)%AREA(IFACE2)
         ! FN      = 0._EB
      END SELECT
      LHSS_CC = LHSS_CC + AF*FN
   ENDDO IFC_LOOP_2

ENDIF DO_SEPARABLE_IF

RETURN
END SUBROUTINE GET_LHS_CUTCELL

! ----------------------------- GET_FH_FROM_PRHS_AND_BCS ----------------------------

SUBROUTINE GET_FH_FROM_PRHS_AND_BCS(NM,DT,CYL_FCT,UNKH,NUNKH,IPZ,F_H)

! NOTE : This routine assumes POINT_TO_MESH has been called and F_H has been ZEROED for the MPI process.

INTEGER, INTENT(IN) :: NM,NUNKH,UNKH,IPZ
REAL(EB),INTENT(IN) :: DT,CYL_FCT
REAL(EB),INTENT(OUT):: F_H(NUNKH)

! Local variables:
INTEGER :: I,J,K,IROW,IW,IIG,JJG,KKG,IOR,ICFACE,IFACE,JFACE,ILH,JLH,KLH,IRC
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(CFACE_TYPE), POINTER :: CFA
REAL(EB) :: IDX, AF, VAL, BCV

! Dummy Assignment.
VAL=DT

! First Source on Cartesian cells with CC_UNKH > 0:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,UNKH)<=0 .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
         ! Row number:
         IROW = CCVAR(I,J,K,UNKH) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START) ! Local numeration.
         ! Add to F_H: If CYL_FCT=0. -> Cartesian coordinates volume (RC(I)=1.).
         !             If CYL_FCT=1. -> Cylindrical coords volume.
         F_H(IROW) = F_H(IROW) + PRHS(I,J,K) * ((1._EB-CYL_FCT)*DY(J) + CYL_FCT*RC(I))*DX(I)*DZ(K)
      ENDDO
   ENDDO
ENDDO

! Rebuild F_H for cut-cells using previously computed CFACE boundary conditions.
IF (CC_IBM) CALL GET_CUTCELL_FH(NM,NUNKH,IPZ,F_H) ! Note: CYL_FCT not used for cut-cells.

! Compute FV in boundary and external CFACEs with DIRICHLET external BCs:
CFACE_LOOP_1 : DO ICFACE=1,N_EXTERNAL_CFACE_CELLS

   CFA => CFACE(ICFACE)
   ! Global matrix solve, skip INTERPOLATED boundaries.
   IF (PRES_FLAG/=ULMAT_FLAG .AND. CFA%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE
   ! Here case where SOLID and OPEN or interpolated are mixed on a boundary:
   IF( CFA%BOUNDARY_TYPE==NULL_BOUNDARY .OR. CFA%BOUNDARY_TYPE==SOLID_BOUNDARY) CYCLE

   IFACE= CFA%CUT_FACE_IND1
   JFACE= CFA%CUT_FACE_IND2
   WC  => WALL(CUT_FACE(IFACE)%IWC)
   EWC => EXTERNAL_WALL(CUT_FACE(IFACE)%IWC)

   ! DIRICHLET boundaries:
   IF_CFACE_DIRICHLET: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN

      ! Gasphase cell indexes:
      BC => BOUNDARY_COORD(WC%BC_INDEX); IF(ZONE_SOLVE(PRESSURE_ZONE(BC%IIG,BC%JJG,BC%KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR
      ! Define centroid to centroid distance, normal to WC:
      IDX=1._EB/(CUT_FACE(IFACE)%XCENHIGH(ABS(IOR),JFACE)-CUT_FACE(IFACE)%XCENLOW(ABS(IOR),JFACE))

      VAL = -2._EB*IDX * CFA%AREA * CFA%PRES_BXN

      ! Row number:
      IROW = CCVAR(IIG,JJG,KKG,UNKH) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START) ! Local numeration.
      IF (IROW <= 0 .AND. CC_IBM) CALL GET_CC_IROW(MESHES(NM),IIG,JJG,KKG,IPZ,IROW)
      IF (IROW <= 0) CYCLE

      ! Add to F_H:
      F_H(IROW) = F_H(IROW) + VAL

   ENDIF IF_CFACE_DIRICHLET

ENDDO CFACE_LOOP_1

! Finally add External Wall cell BCs:
WALL_CELL_LOOP_1: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC => WALL(IW)
   EWC => EXTERNAL_WALL(IW)
   ! Drop if this is a cut-face or NULL Boundary. Dealt with external CFACE.
   IF (WC%CUT_FACE_INDEX>0 .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
   ! Gasphase cell indexes:
   BC => BOUNDARY_COORD(WC%BC_INDEX); IF(ZONE_SOLVE(PRESSURE_ZONE(BC%IIG,BC%JJG,BC%KKG))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR

   ! NEUMANN boundaries:
   IF_NEUMANN: IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN
      ! Define cell size, normal to WC:
      SELECT CASE (IOR)
      CASE(-1) ! -IAXIS oriented, high face of IIG cell.
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG)
         VAL = -BXF(JJG,KKG)*AF
      CASE( 1) ! +IAXIS oriented, low face of IIG cell.
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG)
         VAL =  BXS(JJG,KKG)*AF
      CASE(-2) ! -JAXIS oriented, high face of JJG cell.
         AF  =  DX(IIG)*DZ(KKG)
         VAL = -BYF(IIG,KKG)*AF
      CASE( 2) ! +JAXIS oriented, low face of JJG cell.
         AF  =  DX(IIG)*DZ(KKG)
         VAL =  BYS(IIG,KKG)*AF
      CASE(-3) ! -KAXIS oriented, high face of KKG cell.
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
         VAL = -BZF(IIG,JJG)*AF
      CASE( 3) ! +KAXIS oriented, low face of KKG cell.
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
         VAL =  BZS(IIG,JJG)*AF
      END SELECT

      ! Row number:
      IROW = CCVAR(IIG,JJG,KKG,UNKH) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START) ! Local numeration.
      IF (IROW <= 0 .AND. CC_IBM) THEN
         CALL GET_CC_IROW(MESHES(NM),IIG,JJG,KKG,IPZ,IROW)
         IF (IROW <= 0) CYCLE
      ENDIF

      ! Add to F_H:
      F_H(IROW) = F_H(IROW) + VAL

   ENDIF IF_NEUMANN

   ! DIRICHLET boundaries:
   IF_DIRICHLET: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN
      ! Global matrix solve, skip INTERPOLATED boundaries.
      IF (PRES_FLAG/=ULMAT_FLAG .AND. WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE
      ! Here case where SOLID and OPEN or interpolated are mixed on a boundary:
      IF( WC%BOUNDARY_TYPE==NULL_BOUNDARY  .OR. &
          WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==MIRROR_BOUNDARY ) CYCLE
      ! Define cell size, normal to WC:
      ILH    = 0; JLH = 0; KLH = 0
      SELECT CASE (IOR)
      CASE(-1) ! -IAXIS oriented, high face of IIG cell.
         IDX = RDXN(IIG+ILH); BCV = BXF(JJG,KKG)
         AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG  )) * DZ(KKG)
      CASE( 1) ! +IAXIS oriented, low face of IIG cell.
         ILH = -1; IDX = RDXN(IIG+ILH); BCV = BXS(JJG,KKG)
         AF  = ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*R(IIG-1)) * DZ(KKG)
      CASE(-2) ! -JAXIS oriented, high face of JJG cell.
         IDX = RDYN(JJG+JLH); BCV = BYF(IIG,KKG)
         AF  = DX(IIG)*DZ(KKG)
      CASE( 2) ! +JAXIS oriented, low face of JJG cell.
         JLH = -1; IDX = RDYN(JJG+JLH); BCV = BYS(IIG,KKG)
         AF  = DX(IIG)*DZ(KKG)
      CASE(-3) ! -KAXIS oriented, high face of KKG cell.
         IDX = RDZN(KKG+KLH); BCV = BZF(IIG,JJG)
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
      CASE( 3) ! +KAXIS oriented, low face of KKG cell.
         KLH = -1; IDX = RDZN(KKG+KLH); BCV = BZS(IIG,JJG)
         AF  =  ((1._EB-CYL_FCT)*DY(JJG) + CYL_FCT*RC(IIG  ))* DX(IIG)
      END SELECT
      ! Address case of RC face in the boundary:
      IF (CC_IBM) THEN
         IRC = FCVAR(IIG+ILH,JJG+JLH,KKG+KLH,CC_IDRC,ABS(BC%IOR))
         IF(IRC > 0) IDX = 1._EB / ( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
      ENDIF
      ! Row number:
      IROW = CCVAR(IIG,JJG,KKG,UNKH) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START) ! Local numeration.
      IF (IROW <= 0 .AND. CC_IBM) CALL GET_CC_IROW(MESHES(NM),IIG,JJG,KKG,IPZ,IROW)
      IF (IROW <= 0) CYCLE
      ! Add to F_H:
      F_H(IROW) = F_H(IROW) + (-2._EB*IDX*AF*BCV)

   ENDIF IF_DIRICHLET

ENDDO WALL_CELL_LOOP_1


END SUBROUTINE GET_FH_FROM_PRHS_AND_BCS


! -------------------------------- GET_PRES_CFACE_BCS -------------------------------

SUBROUTINE GET_PRES_CFACE_BCS(NM,T,DT)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

! NOTE : This routine assumes POINT_TO_MESH has been called.

INTEGER, INTENT(IN) :: NM
REAL(EB),INTENT(IN) :: T,DT

! Local Variables:
INTEGER :: ICF,I,J,K,IOR,IFACE,JFACE,NFACE
REAL(EB):: IDX,HN,TSI,VEL_EDDY,TIME_RAMP_FACTOR,P_EXTERNAL,H0,DX_1
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE (CFACE_TYPE), POINTER :: CFA
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: CFA_B1
TYPE(MESH_TYPE), POINTER :: M2
TYPE (CC_CUTFACE_TYPE), POINTER :: CF2
INTEGER :: IIO,JJO,KKO,NOM,IW,ICFO,JCFO
REAL(EB):: H_OTHER,DA_OTHER,DX_OTHER,DY_OTHER,DZ_OTHER,A_FACE

! Dummy assignment.
I=NM; IDX=DT

IF (MESHES(NM)%PREDICTOR) THEN
   UU => U
   VV => V
   WW => W
   HP => H
ELSE
   UU => US
   VV => VS
   WW => WS
   HP => HS
ENDIF

! Apply pressure boundary conditions at external cells.
! If Neumann, CFACE(ICF)%PRES_BXN contains dH/dx.
! If Dirichlet, CFACE(ICF)%PRES_BXN contains H.

! External mesh CFACEs:
CFACE_LOOP_1 : DO ICF=1,N_EXTERNAL_CFACE_CELLS+N_INTWALL_CFACE_CELLS+N_INTERNAL_CFACE_CELLS

   CFA => CFACE(ICF)
   BC  => BOUNDARY_COORD(CFA%BC_INDEX)
   I   = BC%II
   J   = BC%JJ
   K   = BC%KK
   IOR = BC%IOR

   IFACE=CFA%CUT_FACE_IND1
   JFACE=CFA%CUT_FACE_IND2
   NFACE=CUT_FACE(IFACE)%NFACE

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (CFA%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      HN    = 1._EB
      IF (ICF<=N_EXTERNAL_CFACE_CELLS) THEN
         SELECT CASE(IOR)
         CASE( 1); HN = HX(0)
         CASE(-1); HN = HX(IBP1)
         CASE( 2); HN = HY(0)
         CASE(-2); HN = HY(JBP1)
         CASE( 3); HN = HZ(0)
         CASE(-3); HN = HZ(KBP1)
         END SELECT
      ENDIF
      ! dH/Dn where n is pointing out of the gas region.
      CFA%PRES_BXN = -HN*(CUT_FACE(IFACE)%FN(JFACE) + CFA%DUNDT)
   ENDIF IF_NEUMANN

   ! Interpolated and OPEN BCs only in External wall cells:
   EXT_CFACE_IF : IF (ICF<=N_EXTERNAL_CFACE_CELLS) THEN

      IW = CUT_FACE(IFACE)%IWC
      WC=>WALL(IW)
      EWC=>EXTERNAL_WALL(IW)

      ! Interpolated boundary -- set boundary value of H to be average of neighboring cells from previous time step

      INTERPOLATED_ONLY:  IF(PRES_FLAG==ULMAT_FLAG .AND. CFA%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

         NOM     = EWC%NOM
         M2      =>MESHES(NOM)

         H_OTHER = 0._EB; DA_OTHER = 0._EB; DX_OTHER = 0._EB; DY_OTHER = 0._EB; DZ_OTHER = 0._EB
         DX_1 = CUT_FACE(IFACE)%XCENHIGH(ABS(IOR),JFACE)-CUT_FACE(IFACE)%XCENLOW(ABS(IOR),JFACE)
         SELECT CASE(IOR)
            CASE( 1)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF(M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF(M2%FCVAR(IIO,JJO,KKO,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO,JJO,KKO,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 => M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE   = A_FACE + CF2%AREA(JCFO)
                              DX_OTHER = DX_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DY(JJO)*M2%DZ(KKO)
                           DX_OTHER = DX_OTHER + M2%DX(EWC%IIO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER  = H_OTHER/DA_OTHER
               DX_OTHER = DX_OTHER/DA_OTHER
               CFA%PRES_BXN = (DX_OTHER*HP(1,J,K) + DX_1*H_OTHER)/(DX_1+DX_OTHER) + WALL_WORK1(IW)
               BXS(J,K) = CFA%PRES_BXN
            CASE(-1)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF (M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF (M2%FCVAR(IIO-1,JJO,KKO,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO-1,JJO,KKO,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 =>M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE   = A_FACE + CF2%AREA(JCFO)
                              DX_OTHER = DX_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DY(JJO)*M2%DZ(KKO)
                           DX_OTHER = DX_OTHER + M2%DX(EWC%IIO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER  = H_OTHER/DA_OTHER
               DX_OTHER = DX_OTHER/DA_OTHER
               CFA%PRES_BXN = (DX_OTHER*HP(IBAR,J,K) + DX_1*H_OTHER)/(DX_1+DX_OTHER) + WALL_WORK1(IW)
               BXF(J,K) = CFA%PRES_BXN
            CASE( 2)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF (M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF (M2%FCVAR(IIO,JJO,KKO,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO,JJO,KKO,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 =>M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE   = A_FACE + CF2%AREA(JCFO)
                              DY_OTHER = DY_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DX(IIO)*M2%DZ(KKO)
                           DY_OTHER = DY_OTHER + M2%DY(EWC%JJO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER = H_OTHER/DA_OTHER
               DY_OTHER = DY_OTHER/DA_OTHER
               CFA%PRES_BXN = (DY_OTHER*HP(I,1,K) + DX_1*H_OTHER)/(DX_1+DY_OTHER) + WALL_WORK1(IW)
               BYS(I,K) = CFA%PRES_BXN
            CASE(-2)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF (M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF (M2%FCVAR(IIO,JJO-1,KKO,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO,JJO-1,KKO,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 =>M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE   = A_FACE + CF2%AREA(JCFO)
                              DY_OTHER = DY_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DX(IIO)*M2%DZ(KKO)
                           DY_OTHER = DY_OTHER + M2%DY(EWC%JJO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER  = H_OTHER/DA_OTHER
               DY_OTHER = DY_OTHER/DA_OTHER
               CFA%PRES_BXN = (DY_OTHER*HP(I,JBAR,K) + DX_1*H_OTHER)/(DX_1+DY_OTHER) + WALL_WORK1(IW)
               BYF(I,K) = CFA%PRES_BXN
            CASE( 3)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF (M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF (M2%FCVAR(IIO,JJO,KKO,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO,JJO,KKO,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 =>M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE   = A_FACE + CF2%AREA(JCFO)
                              DZ_OTHER = DZ_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DX(IIO)*M2%DY(JJO)
                           DZ_OTHER = DZ_OTHER + M2%DZ(EWC%KKO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER  = H_OTHER/DA_OTHER
               DZ_OTHER = DZ_OTHER/DA_OTHER
               CFA%PRES_BXN = (DZ_OTHER*HP(I,J,1) + DX_1*H_OTHER)/(DX_1+DZ_OTHER) + WALL_WORK1(IW)
               BZS(I,J) = CFA%PRES_BXN
            CASE(-3)

               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                        IF (M2%CELL(M2%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
                        IF (M2%FCVAR(IIO,JJO,KKO-1,CC_FGSC,ABS(IOR))==CC_SOLID) CYCLE

                        ICFO = M2%FCVAR(IIO,JJO,KKO-1,CC_IDCF,ABS(IOR))
                        IF (ICFO>0) THEN
                           CF2 =>M2%CUT_FACE(ICFO); A_FACE = 0._EB
                           DO JCFO=1,CF2%NFACE
                              A_FACE = A_FACE + CF2%AREA(JCFO)
                              DZ_OTHER = DZ_OTHER + &
                              (CF2%XCENHIGH(ABS(IOR),JCFO)-CF2%XCENLOW(ABS(IOR),JCFO))*CF2%AREA(JCFO)
                           ENDDO
                        ELSE
                           A_FACE   = M2%DX(IIO)*M2%DY(JJO)
                           DZ_OTHER = DZ_OTHER + M2%DZ(EWC%KKO_MIN)*A_FACE
                        ENDIF
                        DA_OTHER = DA_OTHER + A_FACE
                        IF (MESHES(NM)%PREDICTOR) H_OTHER = H_OTHER + OMESH(NOM)%H(IIO,JJO,KKO)*A_FACE
                        IF (MESHES(NM)%CORRECTOR) H_OTHER = H_OTHER + OMESH(NOM)%HS(IIO,JJO,KKO)*A_FACE
                     ENDDO
                  ENDDO
               ENDDO
               H_OTHER  = H_OTHER/DA_OTHER
               DZ_OTHER = DZ_OTHER/DA_OTHER
               CFA%PRES_BXN = (DZ_OTHER*HP(I,J,KBAR) + DX_1*H_OTHER)/(DX_1+DZ_OTHER) + WALL_WORK1(IW)
               BZF(I,J) = CFA%PRES_BXN
         END SELECT

      ENDIF INTERPOLATED_ONLY

      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN_IF: IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN

         VT => VENTS(WC%VENT_INDEX)
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         IF (ABS(B1%T_IGN-T_BEGIN)<=TWENTY_EPSILON_EB .AND. VT%PRESSURE_RAMP_INDEX >=1) THEN
            TSI = T
         ELSE
            TSI = T - T_BEGIN
         ENDIF
         TIME_RAMP_FACTOR = EVALUATE_RAMP(TSI,VT%PRESSURE_RAMP_INDEX)
         P_EXTERNAL = TIME_RAMP_FACTOR*VT%DYNAMIC_PRESSURE

         ! Synthetic eddy method for OPEN inflow boundaries
         VEL_EDDY = 0._EB
         IF (VT%N_EDDY>0) THEN
            SELECT CASE(ABS(VT%IOR))
               CASE(1); VEL_EDDY = VT%U_EDDY(J,K)
               CASE(2); VEL_EDDY = VT%V_EDDY(I,K)
               CASE(3); VEL_EDDY = VT%W_EDDY(I,J)
            END SELECT
         ENDIF

         ! Wind inflow boundary conditions

         H0 = 0.5_EB*(U0**2+V0**2+W0**2)
         IF (OPEN_WIND_BOUNDARY) &
         H0 = 0.5_EB*((U_WIND(K)+VEL_EDDY)**2 + (V_WIND(K)+VEL_EDDY)**2 + (W_WIND(K)+VEL_EDDY)**2)

         CFA_B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
         SELECT CASE(IOR)
            CASE( 1)
               IF (UU(0,J,K)<0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(1,J,K)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BXS(J,K) = CFA%PRES_BXN
            CASE(-1)
               IF (UU(IBAR,J,K)>0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(IBAR,J,K)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BXF(J,K) = CFA%PRES_BXN
            CASE( 2)
               IF (VV(I,0,K)<0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(I,1,K)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BYS(I,K) = CFA%PRES_BXN
            CASE(-2)
               IF (VV(I,JBAR,K)>0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(I,JBAR,K)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BYF(I,K) = CFA%PRES_BXN
            CASE( 3)
               IF (WW(I,J,0)<0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(I,J,1)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BZS(I,J) = CFA%PRES_BXN
            CASE(-3)
               IF (WW(I,J,KBAR)>0._EB) THEN
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + KRES(I,J,KBAR)
               ELSE
                  CFA%PRES_BXN = P_EXTERNAL/CFA_B1%RHO_F + H0
               ENDIF
               BZF(I,J) = CFA%PRES_BXN
         END SELECT

      ENDIF OPEN_IF

   ENDIF EXT_CFACE_IF

ENDDO CFACE_LOOP_1

END SUBROUTINE GET_PRES_CFACE_BCS

! --------------------------------- GET_CUTCELL_DDDT --------------------------------

SUBROUTINE GET_CUTCELL_DDDT(M,DT,NM)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER,  INTENT(IN) :: NM

! Local Variables:
REAL(EB) :: RDT, PRFCT, VOL, DIVVOL, DPCC, DIV_JCC
INTEGER  :: I,J,K,IPZ,NCELL,ICC,JCC
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P

RDT = 1._EB/DT

SELECT CASE(M%PREDICTOR)
   CASE(.TRUE.)
      D_PBAR_DT_P => M%D_PBAR_DT_S
      PRFCT = 1._EB
   CASE(.FALSE.)
      D_PBAR_DT_P => M%D_PBAR_DT
      PRFCT = 0._EB
END SELECT

PRED_CORR_IF : IF (M%PREDICTOR) THEN

   ICC_LOOP_1 : DO ICC=1,M%N_CUTCELL_MESH
      I      = M%CUT_CELL(ICC)%IJK(IAXIS)
      J      = M%CUT_CELL(ICC)%IJK(JAXIS)
      K      = M%CUT_CELL(ICC)%IJK(KAXIS)
      IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE ICC_LOOP_1
      IPZ    = M%PRESSURE_ZONE(I,J,K)
      NCELL  = M%CUT_CELL(ICC)%NCELL
      DIVVOL = 0._EB; DPCC   = 0._EB; VOL    = 0._EB
      IF (ONE_UNKH_PER_CUTCELL) THEN ! DDDTVOL(JCC) defined per cut-cell.
         DO JCC=1,NCELL
            CALL GET_VELOC_DIVERGENCE_CUTCELL(M, &
            ICC,JCC,0._EB,DIV_JCC) ! Velocity divg of U,V,W -> PRFCT=0._EB
            M%CUT_CELL(ICC)%DVOL_PR(JCC) = DIV_JCC*M%CUT_CELL(ICC)%VOLUME(JCC)
            DIVVOL = DIVVOL + M%CUT_CELL(ICC)%DVOL_PR(JCC)
            ! Thermodynamic divergence * vol:
            DPCC= ( (1._EB-PRFCT)*M%CUT_CELL(ICC)%D(JCC) + PRFCT*M%CUT_CELL(ICC)%DS(JCC) )*M%CUT_CELL(ICC)%VOLUME(JCC)
            ! Add Pressure derivative to divergence:
            M%CUT_CELL(ICC)%DDDTVOL(JCC)  = (DPCC-M%CUT_CELL(ICC)%DVOL_PR(JCC))*RDT
         ENDDO
      ELSE
         DO JCC=1,NCELL
            VOL  = VOL + M%CUT_CELL(ICC)%VOLUME(JCC)
            CALL GET_VELOC_DIVERGENCE_CUTCELL(M, &
            ICC,JCC,0._EB,DIV_JCC) ! Velocity divg of U,V,W -> PRFCT=0._EB
            DIVVOL = DIVVOL + DIV_JCC*M%CUT_CELL(ICC)%VOLUME(JCC)
            ! Thermodynamic divergence * vol:
            DPCC= DPCC + ( (1._EB-PRFCT)*M%CUT_CELL(ICC)%D(JCC) + PRFCT*M%CUT_CELL(ICC)%DS(JCC) )*M%CUT_CELL(ICC)%VOLUME(JCC)
         ENDDO
         ! Define average DDDT for CUT_CELL(ICC):
         M%CUT_CELL(ICC)%DDDTVOL(1)  = (DPCC-DIVVOL)*RDT
         M%CUT_CELL(ICC)%DVOL_PR(1)  = DIVVOL
      ENDIF

   ENDDO ICC_LOOP_1

ELSEIF (M%CORRECTOR) THEN PRED_CORR_IF

   ICC_LOOP_2 : DO ICC=1,M%N_CUTCELL_MESH
      I      = M%CUT_CELL(ICC)%IJK(IAXIS)
      J      = M%CUT_CELL(ICC)%IJK(JAXIS)
      K      = M%CUT_CELL(ICC)%IJK(KAXIS)
      IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE ICC_LOOP_2
      IPZ    = M%PRESSURE_ZONE(I,J,K)
      NCELL  = M%CUT_CELL(ICC)%NCELL
      DIVVOL = 0._EB; DPCC   = 0._EB; VOL    = 0._EB
      IF (ONE_UNKH_PER_CUTCELL) THEN ! DDDTVOL(JCC) defined per cut-cell.
         DO JCC=1,NCELL
            CALL GET_VELOC_DIVERGENCE_CUTCELL(M, &
            ICC,JCC,1._EB,DIV_JCC) ! Velocity divg of US,VS,WS -> PRFCT=1._EB
            DPCC= ( (1._EB-PRFCT)*M%CUT_CELL(ICC)%D(JCC) + PRFCT*M%CUT_CELL(ICC)%DS(JCC) )*M%CUT_CELL(ICC)%VOLUME(JCC)
            M%CUT_CELL(ICC)%DDDTVOL(JCC)  =  &
            (2._EB*DPCC-(DIV_JCC*M%CUT_CELL(ICC)%VOLUME(JCC)+M%CUT_CELL(ICC)%DVOL_PR(JCC)))*RDT
         ENDDO
      ELSE
         DO JCC=1,NCELL
            VOL  = VOL + M%CUT_CELL(ICC)%VOLUME(JCC)
            CALL GET_VELOC_DIVERGENCE_CUTCELL(M, &
            ICC,JCC,1._EB,DIV_JCC) ! Velocity divg of US,VS,WS -> PRFCT=1._EB
            DIVVOL = DIVVOL + DIV_JCC*M%CUT_CELL(ICC)%VOLUME(JCC)
            ! Thermodynamic divergence * vol:
            DPCC= DPCC + ( (1._EB-PRFCT)*M%CUT_CELL(ICC)%D(JCC) + PRFCT*M%CUT_CELL(ICC)%DS(JCC) )*M%CUT_CELL(ICC)%VOLUME(JCC)
         ENDDO
         ! Define average DDDT for CUT_CELL(ICC):
         M%CUT_CELL(ICC)%DDDTVOL(1)  = (2._EB*DPCC-(DIVVOL+M%CUT_CELL(ICC)%DVOL_PR(1)))*RDT
      ENDIF
   ENDDO ICC_LOOP_2

ENDIF PRED_CORR_IF

RETURN
END SUBROUTINE GET_CUTCELL_DDDT

! --------------------- GET_BOUNDFACE_GEOM_INFO_H --------------------------------
SUBROUTINE GET_BOUNDFACE_GEOM_INFO_H

! Work deferred.

RETURN
END SUBROUTINE GET_BOUNDFACE_GEOM_INFO_H

! ----------------------------- GET_CUTCELL_HP ------------------------------------


SUBROUTINE GET_CUTCELL_HP(NM,IPZ,HP)

INTEGER, INTENT(IN) :: NM,IPZ
REAL(EB), INTENT(INOUT), POINTER, DIMENSION(:,:,:) :: HP

! Local Variables:
INTEGER :: I,J,K,IROW,ICC

IF (MESHES(NM)%PREDICTOR) THEN

   ! Note does not take into account ONE_UNKH_PER_CUTCELL:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => MESHES(NM)%CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      IF(CELL(CELL_INDEX(I,J,K))%SOLID .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IROW = CC%UNKH(1) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START)
      ! Assign to cut-cell H:
      CC%H(1:CC%NCELL) = -ZONE_SOLVE(IPZ)%X_H(IROW)
      ! Assign to HP:
      HP(I,J,K) = -ZONE_SOLVE(IPZ)%X_H(IROW)
   ENDDO

ELSE

   ! Note does not take into account ONE_UNKH_PER_CUTCELL:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => MESHES(NM)%CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      IF(CELL(CELL_INDEX(I,J,K))%SOLID .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      IROW     = CC%UNKH(1) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START)
      ! Assign to cut-cell HS:
      CC%HS(1:CC%NCELL) = -ZONE_SOLVE(IPZ)%X_H(IROW)
      ! Assign to HP:
      HP(I,J,K) = -ZONE_SOLVE(IPZ)%X_H(IROW)
   ENDDO

ENDIF

RETURN
END SUBROUTINE GET_CUTCELL_HP

! ----------------------------- GET_CUTCELL_FH ------------------------------------

SUBROUTINE GET_CUTCELL_FH(NM,NUNKH,IPZ,F_H)

! NOTE : Assumes POINT_TO_MESH(NM) has been called.

INTEGER, INTENT(IN)     :: NM,NUNKH,IPZ
REAL(EB), INTENT(INOUT) :: F_H(1:NUNKH)

! Local Variables:
INTEGER :: IROW,ICC,JCC,I,J,K
REAL(EB):: DIV_FN, DIV_FN_VOL

CUTCELL_LOOP_A : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   CC   => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
   IF(CELL(CELL_INDEX(I,J,K))%SOLID .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
   IROW =  CC%UNKH(1) - ZONE_SOLVE(IPZ)%UNKH_IND(NM_START)
   ! Here we add div(F) in the cut-cell and DDDT:
   DIV_FN_VOL = 0._EB
   DO JCC=1,CC%NCELL
      CALL GET_FN_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC,DIV_FN, &
        SUBSTRACT_BAROCLINIC=.FALSE.)
      DIV_FN_VOL = DIV_FN_VOL + DIV_FN*CC%VOLUME(JCC)
   ENDDO
   ! Add to F_H:
   F_H(IROW) = -(CC%DDDTVOL(1) + DIV_FN_VOL)
ENDDO CUTCELL_LOOP_A

RETURN
END SUBROUTINE GET_CUTCELL_FH

! ---------------------------- GET_H_MATRIX_CC ------------------------------------

SUBROUTINE GET_H_MATRIX_CC(NM,NM1,IPZ)

! This routine assumes the calling subroutine has called POINT_TO_MESH for NM.

INTEGER, INTENT(IN) :: NM,NM1,IPZ

! Local Variables:
INTEGER :: X1AXIS,IFACE,ICF,I,J,K,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),ILOC,JLOC,JCOL,IROW, &
           LOCROW,LOCROW_1,LOCROW_2,IW,NOM,IIO,JJO,KKO,ICC_EXT,IND_INT,IND_EXT,ICFO,ISHF(IAXIS:KAXIS), &
           ICC,II,JJ,KK,IIG,JJG,KKG
REAL(EB) :: AF,IDX,BIJ,KFACE(2,2),X_FACE,D_INT,D_EXT,AF_EXT
TYPE(CC_CUTFACE_TYPE), POINTER :: CF
TYPE(CC_RCFACE_TYPE), POINTER :: RCF
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! X direction bounds:
ILO_FACE = 0                ! Low mesh boundary face index.
IHI_FACE = IBAR             ! High mesh boundary face index.
ILO_CELL = ILO_FACE + 1     ! First internal cell index. See notes.
IHI_CELL = IHI_FACE         ! Last internal cell index.
! Y direction bounds:
JLO_FACE = 0                ! Low mesh boundary face index.
JHI_FACE = JBAR             ! High mesh boundary face index.
JLO_CELL = JLO_FACE + 1     ! First internal cell index. See notes.
JHI_CELL = JHI_FACE         ! Last internal cell index.
! Z direction bounds:
KLO_FACE = 0                ! Low mesh boundary face index.
KHI_FACE = KBAR             ! High mesh boundary face index.
KLO_CELL = KLO_FACE + 1     ! First internal cell index. See notes.
KHI_CELL = KHI_FACE         ! Last internal cell index.

! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
   RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE));
   I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
   IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
   ! Unknowns on related cells:
   IND(LOW_IND)  = RCF%UNKH(LOW_IND)
   IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
   IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
   IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
   ! Row ind(1),ind(2):
   LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
   SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF  = DY(J)*DZ(K)
         IDX = RDXN(I)
         IF ( I == ILO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
         IF ( I == IHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      CASE(JAXIS)
         AF  = DX(I)*DZ(K)
         IDX = RDYN(J)
         IF ( J == JLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
         IF ( J == JHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
      CASE(KAXIS)
         AF  = DX(I)*DY(J)
         IDX = RDZN(K)
         IF ( K == KLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
         IF ( K == KHI_FACE)  LOCROW_2 =  LOW_IND ! Only low side unknown row.
   ENDSELECT
   ! Skip all boundary RC faces (handled by EXTERNAL_WALL loop below)
   IF (LOCROW_1 == LOCROW_2) CYCLE
   IDX = 1._EB / ( RCF%XCEN(X1AXIS,HIGH_IND) - RCF%XCEN(X1AXIS,LOW_IND) )

   ! Now add to Adiff corresponding coeff:
   BIJ   = IDX*AF
   !    Cols 1,2: ind(LOW_IND) ind(HIGH_IND), Rows 1,2: ind_loc(LOW_IND) ind_loc(HIGH_IND)
   KFACE(1,1) = BIJ; KFACE(2,1) =-BIJ; KFACE(1,2) =-BIJ; KFACE(2,2) = BIJ
   DO ILOC=LOCROW_1,LOCROW_2   ! Local row number in Kface
      IROW=IND_LOC(ILOC)       ! Process Local Unknown number.
      DO JLOC=LOW_IND,HIGH_IND ! Local col number in Kface
         ! Find column position on-the-fly
         DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
            IF (IND(JLOC) == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
               ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) + KFACE(ILOC,JLOC)
               EXIT
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Now Gasphase CUT_FACES:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   CF =>  MESHES(NM)%CUT_FACE(ICF); IF ( CF%STATUS /= CC_GASPHASE ) CYCLE
   I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE

   ! Handle external cut-faces (at mesh boundaries) - same-level and refined
   IF( CF%IWC > 0 ) THEN
      WC=>MESHES(NM)%WALL(CF%IWC)
      IF (.NOT.( WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) ) CYCLE
      EWC => EXTERNAL_WALL(CF%IWC)
      NOM = EWC%NOM

      ! Determine which side is internal
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF (I == ILO_FACE) THEN; LOCROW = HIGH_IND; X_FACE = XS
            ELSE;                    LOCROW = LOW_IND;  X_FACE = XF
            ENDIF
         CASE(JAXIS)
            IF (J == JLO_FACE) THEN; LOCROW = HIGH_IND; X_FACE = YS
            ELSE;                    LOCROW = LOW_IND;  X_FACE = YF
            ENDIF
         CASE(KAXIS)
            IF (K == KLO_FACE) THEN; LOCROW = HIGH_IND; X_FACE = ZS
            ELSE;                    LOCROW = LOW_IND;  X_FACE = ZF
            ENDIF
      END SELECT

      DO IFACE=1,CF%NFACE
         ! Internal unknown
         IND_INT = CF%UNKH(LOCROW,IFACE)
         IROW = IND_INT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
         AF = CF%AREA(IFACE)
         ! Distance from internal cell centroid to interface
         IF (LOCROW == HIGH_IND) THEN
            D_INT = ABS(CF%XCENHIGH(X1AXIS,IFACE) - X_FACE)
         ELSE
            D_INT = ABS(CF%XCENLOW(X1AXIS,IFACE) - X_FACE)
         ENDIF

         ! Loop over external cells (works for both same-level and refined)
         DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
            DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
               DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
                  ICC_EXT = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
                  IND_EXT = OMESH(NOM)%MUNKH(IIO,JJO,KKO)
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS)
                        AF_EXT = MESHES(NOM)%DY(JJO) * MESHES(NOM)%DZ(KKO)
                        D_EXT = ABS(MESHES(NOM)%XC(IIO) - X_FACE)
                  CASE(JAXIS)
                        AF_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DZ(KKO)
                        D_EXT = ABS(MESHES(NOM)%YC(JJO) - X_FACE)
                  CASE(KAXIS)
                        AF_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DY(JJO)
                        D_EXT = ABS(MESHES(NOM)%ZC(KKO) - X_FACE)
                  END SELECT
                  IF (ICC_EXT > 0) THEN
                     ! Note: IND_EXT comes from OMESH(NOM)%MUNKH which was populated
                     ! via COPY_CC_MUNKH_TO_UNKH from communicated HS data
                     D_EXT = ABS(MESHES(NOM)%CUT_CELL(ICC_EXT)%XYZCEN(X1AXIS,1) - X_FACE)
                  ENDIF
                  ! Check for cut-face in NOM
                  ISHF(IAXIS:KAXIS) = 0; IF(BOUNDARY_COORD(WC%BC_INDEX)%IOR < 0) ISHF(X1AXIS) = -1
                  ICFO = MESHES(NOM)%FCVAR(IIO+ISHF(IAXIS),JJO+ISHF(JAXIS),KKO+ISHF(KAXIS),CC_IDCF,X1AXIS)
                  IF (ICFO > 0) AF_EXT = SUM(MESHES(NOM)%CUT_FACE(ICFO)%AREA(1:MESHES(NOM)%CUT_FACE(ICFO)%NFACE))
                  IF (EWC%AREA_RATIO < 0.99_EB) AF_EXT = AF
                  BIJ = AF_EXT / (D_INT + D_EXT)
                  ! Add to diagonal (IND_INT)
                  DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
                     IF (IND_INT == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
                        ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) + BIJ
                        EXIT
                     ENDIF
                  ENDDO
                  ! Add to off-diagonal (IND_EXT)
                  DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
                     IF (IND_EXT == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
                        ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) - BIJ
                        EXIT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CYCLE  ! Skip normal processing - external cut-face handled
   ENDIF

   ! Normal processing (internal cut-faces only)
   SELECT CASE(X1AXIS)
      CASE(IAXIS); IDX = RDXN(I)
      CASE(JAXIS); IDX = RDYN(J)
      CASE(KAXIS); IDX = RDZN(K)
   ENDSELECT
   DO IFACE=1,CF%NFACE
      ! Unknowns on related cells:
      IND(LOW_IND)     = CF%UNKH(LOW_IND,IFACE)
      IND(HIGH_IND)    = CF%UNKH(HIGH_IND,IFACE)
      IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
      AF = CF%AREA(IFACE)
      IDX= 1._EB/ ( CF%XCENHIGH(X1AXIS,IFACE) - CF%XCENLOW(X1AXIS, IFACE) )

      ! Now add to Adiff corresponding coeff:
      BIJ   = IDX*AF
      !    Cols 1,2: ind(LOW_IND) ind(HIGH_IND), Rows 1,2: ind_loc(LOW_IND) ind_loc(HIGH_IND)
      KFACE(1,1) = BIJ; KFACE(2,1) =-BIJ; KFACE(1,2) =-BIJ; KFACE(2,2) = BIJ
      DO ILOC=LOW_IND,HIGH_IND ! Local row number in Kface (both rows for internal faces)
         IROW=IND_LOC(ILOC)
         DO JLOC=LOW_IND,HIGH_IND ! Local col number in Kface
            ! Find column position on-the-fly
            DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
               IF (IND(JLOC) == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
                  ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) + KFACE(ILOC,JLOC)
                  EXIT
               ENDIF
           ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

! Handle RC faces at mesh boundaries (via EXTERNAL_WALL loop)
! This handles RC faces connecting cut-cells to external meshes (same-level and refined)
! RC face = regular Cartesian face where at least one side is a cut-cell
DO IW = 1, N_EXTERNAL_WALL_CELLS
   WC => WALL(IW)
   IF (.NOT.(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   IF (WC%CUT_FACE_INDEX > 0) CYCLE  ! Skip cut-faces (handled above)

   EWC => EXTERNAL_WALL(IW)
   NOM = EWC%NOM; IF (NOM < 1) CYCLE

   ! Get cell indices from boundary coord
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   II = BC%II; JJ = BC%JJ; KK = BC%KK         ! Ghost cell
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG   ! Internal cell
   X1AXIS = ABS(BC%IOR)

   IF (ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT /= IPZ) CYCLE

   ! Check if this is an RC face: at least one side must be a cut-cell
   ! Skip if neither side is a cut-cell (handled by pres.f90)
   ICC = CCVAR(IIG,JJG,KKG,CC_IDCC)
   IF (ICC < 1 .AND. CCVAR(II,JJ,KK,CC_IDCC) < 1) CYCLE

   ! Get internal unknown and distance
   IF (ICC > 0) THEN
      IND_INT = CUT_CELL(ICC)%UNKH(1)
   ELSE
      IND_INT = CCVAR(IIG,JJG,KKG,CC_UNKH)
   ENDIF
   IROW = IND_INT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)

   ! Get interface location and internal cell distance
   SELECT CASE(X1AXIS)
      CASE(IAXIS)
         IF (BC%IOR > 0) THEN; X_FACE = XS; ELSE; X_FACE = XF; ENDIF
         IF (ICC > 0) THEN; D_INT = ABS(CUT_CELL(ICC)%XYZCEN(IAXIS,1) - X_FACE)
         ELSE;              D_INT = ABS(XC(IIG) - X_FACE)
         ENDIF
         AF = DY(JJG) * DZ(KKG)
      CASE(JAXIS)
         IF (BC%IOR > 0) THEN; X_FACE = YS; ELSE; X_FACE = YF; ENDIF
         IF (ICC > 0) THEN; D_INT = ABS(CUT_CELL(ICC)%XYZCEN(JAXIS,1) - X_FACE)
         ELSE;              D_INT = ABS(YC(JJG) - X_FACE)
         ENDIF
         AF = DX(IIG) * DZ(KKG)
      CASE(KAXIS)
         IF (BC%IOR > 0) THEN; X_FACE = ZS; ELSE; X_FACE = ZF; ENDIF
         IF (ICC > 0) THEN; D_INT = ABS(CUT_CELL(ICC)%XYZCEN(KAXIS,1) - X_FACE)
         ELSE;              D_INT = ABS(ZC(KKG) - X_FACE)
         ENDIF
         AF = DX(IIG) * DY(JJG)
   END SELECT

   ! Loop over external cells (refinement)
   DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
      DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
         DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
            ICC_EXT = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
            IND_EXT = OMESH(NOM)%MUNKH(IIO,JJO,KKO)
            SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  AF_EXT = MESHES(NOM)%DY(JJO) * MESHES(NOM)%DZ(KKO)
                  D_EXT = ABS(MESHES(NOM)%XC(IIO) - X_FACE)
               CASE(JAXIS)
                  AF_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DZ(KKO)
                  D_EXT = ABS(MESHES(NOM)%YC(JJO) - X_FACE)
               CASE(KAXIS)
                  AF_EXT = MESHES(NOM)%DX(IIO) * MESHES(NOM)%DY(JJO)
                  D_EXT = ABS(MESHES(NOM)%ZC(KKO) - X_FACE)
            END SELECT
            IF (ICC_EXT > 0) THEN
               ! Note: IND_EXT comes from OMESH(NOM)%MUNKH which was populated
               ! via COPY_CC_MUNKH_TO_UNKH from communicated HS data
               D_EXT = ABS(MESHES(NOM)%CUT_CELL(ICC_EXT)%XYZCEN(X1AXIS,1) - X_FACE)
            ENDIF
            ! Check for cut-face in NOM
            ISHF(IAXIS:KAXIS) = 0; IF (BC%IOR < 0) ISHF(X1AXIS) = -1
            ICFO = MESHES(NOM)%FCVAR(IIO+ISHF(IAXIS),JJO+ISHF(JAXIS),KKO+ISHF(KAXIS),CC_IDCF,X1AXIS)
            IF (ICFO > 0) AF_EXT = SUM(MESHES(NOM)%CUT_FACE(ICFO)%AREA(1:MESHES(NOM)%CUT_FACE(ICFO)%NFACE))
            IF (EWC%AREA_RATIO < 0.99_EB) AF_EXT = AF
            BIJ = AF_EXT / (D_INT + D_EXT)
            ! Add to diagonal (IND_INT)
            DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
               IF (IND_INT == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
                  ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) + BIJ
                  EXIT
               ENDIF
            ENDDO
            ! Add to off-diagonal (IND_EXT)
            DO JCOL = 1, ZONE_SOLVE(IPZ)%ROW_H(IROW)%NNZ
               IF (IND_EXT == ZONE_SOLVE(IPZ)%ROW_H(IROW)%JD(JCOL)) THEN
                  ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) = ZONE_SOLVE(IPZ)%ROW_H(IROW)%D(JCOL) - BIJ
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
RETURN
END SUBROUTINE GET_H_MATRIX_CC


! -------------------------- GET_CC_MATRIXGRAPH_H ---------------------------------

SUBROUTINE GET_CC_MATRIXGRAPH_H(NM,NM1,IPZ,LOOP_FLAG)

INTEGER, INTENT(IN) :: NM,NM1,IPZ
LOGICAL, INTENT(IN) :: LOOP_FLAG

! Local Variables:
INTEGER :: X1AXIS,IFACE,ICF,I,J,K,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND)
INTEGER :: LOCROW_1,LOCROW_2,LOCROW,IIND,NII,ILOC
INTEGER :: IW,IIG,JJG,KKG,II,JJ,KK,NOM,IIO,JJO,KKO,ICC_INT,ICC_EXT,IND_INT,IND_EXT

! X direction bounds:
ILO_FACE = 0                    ! Low mesh boundary face index.
IHI_FACE = MESHES(NM)%IBAR      ! High mesh boundary face index.
ILO_CELL = ILO_FACE + 1     ! First internal cell index. See notes.
IHI_CELL = IHI_FACE ! Last internal cell index.

! Y direction bounds:
JLO_FACE = 0                    ! Low mesh boundary face index.
JHI_FACE = MESHES(NM)%JBAR      ! High mesh boundary face index.
JLO_CELL = JLO_FACE + 1     ! First internal cell index. See notes.
JHI_CELL = JHI_FACE ! Last internal cell index.

! Z direction bounds:
KLO_FACE = 0                    ! Low mesh boundary face index.
KHI_FACE = MESHES(NM)%KBAR      ! High mesh boundary face index.
KLO_CELL = KLO_FACE + 1     ! First internal cell index. See notes.
KHI_CELL = KHI_FACE ! Last internal cell index.

LOOP_FLAG_COND : IF ( LOOP_FLAG ) THEN ! MESH_LOOP_1 in calling routine.
   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE));
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      ! Unknowns on related cells:
      IND(LOW_IND)  = RCF%UNKH(LOW_IND)
      IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
      IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! Row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == ILO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == JLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == KLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KHI_FACE)  LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      ! Add to global matrix arrays:
      CALL ADD_INPLACE_NNZ_H_WHLDOM(LOCROW_1,LOCROW_2,IND,IND_LOC,IPZ)
   ENDDO

   DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      CF => MESHES(NM)%CUT_FACE(ICF); IF (CF%STATUS/=CC_GASPHASE) CYCLE
      ! Drop if cut-face on a wall-cell, and type different than INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY.
      IF(CF%IWC>0) THEN
         WC=>MESHES(NM)%WALL(CF%IWC)
         IF (.NOT.(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      ENDIF
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == ILO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == JLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == KLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KHI_FACE)  LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT

      ! External cut-face at mesh boundary - need to add entries for ALL external cells
      ! This handles both same-level and refinement cases
      IF (CF%IWC > 0 .AND. LOCROW_1 == LOCROW_2) THEN
         EWC => EXTERNAL_WALL(CF%IWC)
         NOM = EWC%NOM
         IF (NOM > 0) THEN
            DO IFACE=1,CF%NFACE
               IND_INT = CF%UNKH(LOCROW_1,IFACE)
               IND_LOC(LOW_IND) = IND_INT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
               ! Loop over ALL external cells
               DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
                  DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
                     DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
                        ! External unknown from communicated MUNKH array
                        ! (MUNKH is populated from HS which contains cut-cell UNKH)
                        IND_EXT = OMESH(NOM)%MUNKH(IIO,JJO,KKO)
                        IF (IND_EXT < 1) CYCLE
                        IND_LOC(HIGH_IND) = IND_EXT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
                        IND(LOW_IND)  = IND_INT
                        IND(HIGH_IND) = IND_EXT
                        CALL ADD_INPLACE_NNZ_H_WHLDOM(LOW_IND,LOW_IND,IND,IND_LOC,IPZ)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         CYCLE  ! Skip normal processing for this external cut-face
      ENDIF

      ! Internal cut-faces (normal processing)
      DO IFACE=1,CF%NFACE
         !% Unknowns on related cells:
         IND(LOW_IND)     = CF%UNKH(LOW_IND,IFACE)
         IND(HIGH_IND)    = CF%UNKH(HIGH_IND,IFACE)
         IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! Row indexes refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
         ! Add to global matrix arrays:
         CALL ADD_INPLACE_NNZ_H_WHLDOM(LOCROW_1,LOCROW_2,IND,IND_LOC,IPZ)
      ENDDO
   ENDDO

   ! Handle refinement interfaces via EXTERNAL_WALL loop (RC faces only).
   ! For refinement, one coarse cell connects to multiple fine cells in the neighboring mesh.
   ! We need to add graph entries for ALL external cells, not just one.
   ! Note: Cut-faces are handled by the CUT_FACE loop above.
   DO IW = 1, N_EXTERNAL_WALL_CELLS
      WC => WALL(IW)
      IF (.NOT.(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      IF (WC%CUT_FACE_INDEX > 0) CYCLE  ! Skip cut-faces (handled above)
      EWC => EXTERNAL_WALL(IW)
      NOM = EWC%NOM; IF (NOM < 1) CYCLE

      BC => BOUNDARY_COORD(WC%BC_INDEX)
      IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG   ! Internal cell
      II  = BC%II;  JJ  = BC%JJ;  KK  = BC%KK    ! Ghost cell

      IF (ZONE_SOLVE(PRESSURE_ZONE(IIG,JJG,KKG))%CONNECTED_ZONE_PARENT /= IPZ) CYCLE

      ! Skip if both cells are regular gas (handled by pres.f90)
      ICC_INT = CCVAR(IIG,JJG,KKG,CC_IDCC)
      ICC_EXT = CCVAR(II,JJ,KK,CC_IDCC)
      IF (ICC_INT < 1 .AND. ICC_EXT < 1) CYCLE

      ! Get internal unknown
      IF (ICC_INT > 0) THEN
         IND_INT = CUT_CELL(ICC_INT)%UNKH(1)
      ELSE
         IND_INT = CCVAR(IIG,JJG,KKG,CC_UNKH)
      ENDIF
      IND_LOC(LOW_IND) = IND_INT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)

      ! Loop over ALL cells in the neighboring mesh that share this boundary
      DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
         DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
            DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
               ! External unknown from communicated MUNKH array
               ! (MUNKH is populated from HS which contains cut-cell UNKH)
               IND_EXT = OMESH(NOM)%MUNKH(IIO,JJO,KKO)
               IF (IND_EXT < 1) CYCLE

               IND_LOC(HIGH_IND) = IND_EXT - ZONE_SOLVE(IPZ)%UNKH_IND(NM1)

               ! Add graph entries: internal row -> external column
               IND(LOW_IND)  = IND_INT
               IND(HIGH_IND) = IND_EXT
               CALL ADD_INPLACE_NNZ_H_WHLDOM(LOW_IND,LOW_IND,IND,IND_LOC,IPZ)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

ELSE ! MESH_LOOP_2 in calling routine.

   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   DO IFACE=1,MESHES(NM)%CC_NRCFACE_H
      RCF => RC_FACE(MESHES(NM)%RCF_H(IFACE));
      I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      ! Unknowns on related cells:
      IND(LOW_IND)  = RCF%UNKH(LOW_IND)
      IND(HIGH_IND) = RCF%UNKH(HIGH_IND)
      IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! Row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == ILO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == JLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == KLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KHI_FACE)  LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      RCF%JDH(1:2,1:2) = 0
      ! Add to global matrix arrays:
      DO LOCROW=LOCROW_1,LOCROW_2
         DO IIND=LOW_IND,HIGH_IND
            NII = ZONE_SOLVE(IPZ)%ROW_H(IND_LOC(LOCROW))%NNZ
            DO ILOC=1,NII
               IF ( IND(IIND) == ZONE_SOLVE(IPZ)%ROW_H(IND_LOC(LOCROW))%JD(ILOC) ) THEN
                   RCF%JDH(LOCROW,IIND) = ILOC
                   EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Now Gasphase CUT_FACES:
   DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      CF => MESHES(NM)%CUT_FACE(ICF); IF (CF%STATUS/=CC_GASPHASE) CYCLE
      ! Drop if cut-face on a wall-cell, and type different than INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY.
      IF(CF%IWC>0) THEN
         WC=>MESHES(NM)%WALL(CF%IWC)
         IF (.NOT.(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      ENDIF
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND; LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF ( I == ILO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( I == IHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(JAXIS)
            IF ( J == JLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( J == JHI_FACE ) LOCROW_2 =  LOW_IND ! Only low side unknown row.
         CASE(KAXIS)
            IF ( K == KLO_FACE ) LOCROW_1 = HIGH_IND ! Only high side unknown row.
            IF ( K == KHI_FACE)  LOCROW_2 =  LOW_IND ! Only low side unknown row.
      ENDSELECT
      CF%JDH(:,:,:) = 0
      DO IFACE=1,CF%NFACE
         !% Unknowns on related cells:
         IND(LOW_IND)     = CF%UNKH(LOW_IND,IFACE)
         IND(HIGH_IND)    = CF%UNKH(HIGH_IND,IFACE)
         IND_LOC(LOW_IND) = IND(LOW_IND) - ZONE_SOLVE(IPZ)%UNKH_IND(NM1) ! Row indexes refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- ZONE_SOLVE(IPZ)%UNKH_IND(NM1)
         ! Add to global matrix arrays:
         DO LOCROW=LOCROW_1,LOCROW_2
            DO IIND=LOW_IND,HIGH_IND
               NII = ZONE_SOLVE(IPZ)%ROW_H(IND_LOC(LOCROW))%NNZ
               DO ILOC=1,NII
                  IF ( IND(IIND) == ZONE_SOLVE(IPZ)%ROW_H(IND_LOC(LOCROW))%JD(ILOC) ) THEN
                        CF%JDH(LOCROW,IIND,IFACE) = ILOC
                        EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDIF LOOP_FLAG_COND

RETURN
END SUBROUTINE GET_CC_MATRIXGRAPH_H

! ------------------------ ADD_INPLACE_NNZ_H_WHLDOM -------------------------------

SUBROUTINE ADD_INPLACE_NNZ_H_WHLDOM(LOCROW_1,LOCROW_2,IND,IND_LOC,IPZ)

INTEGER, INTENT(IN) :: LOCROW_1,LOCROW_2,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),IPZ

! Local Variables:
INTEGER LOCROW, IIND, ILOC
LOGICAL INLIST

LOCROW_LOOP : DO LOCROW=LOCROW_1,LOCROW_2
   DO IIND=LOW_IND,HIGH_IND
      ! Populate variable per-row storage ROW_H only:
      ILOC = IND_LOC(LOCROW)
      IF (ILOC>=1 .AND. ILOC<=SIZE(ZONE_SOLVE(IPZ)%ROW_H)) THEN
         IF (.NOT. ALLOCATED(ZONE_SOLVE(IPZ)%ROW_H(ILOC)%JD)) THEN
            ALLOCATE(ZONE_SOLVE(IPZ)%ROW_H(ILOC)%JD(7))
            ALLOCATE(ZONE_SOLVE(IPZ)%ROW_H(ILOC)%D(7))
            ZONE_SOLVE(IPZ)%ROW_H(ILOC)%JD(1) = IND(IIND)
            ZONE_SOLVE(IPZ)%ROW_H(ILOC)%D(1)  = 0._EB
            ZONE_SOLVE(IPZ)%ROW_H(ILOC)%NNZ   = 1
         ELSE
            INLIST = ANY(ZONE_SOLVE(IPZ)%ROW_H(ILOC)%JD(1:ZONE_SOLVE(IPZ)%ROW_H(ILOC)%NNZ) == IND(IIND))
            IF (.NOT. INLIST) THEN
               CALL INSERT_SORTED_INT_REAL(            &
                     ZONE_SOLVE(IPZ)%ROW_H(ILOC)%JD,   &
                     ZONE_SOLVE(IPZ)%ROW_H(ILOC)%D,    &
                     ZONE_SOLVE(IPZ)%ROW_H(ILOC)%NNZ,  &
                     IND(IIND))
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO LOCROW_LOOP

RETURN
END SUBROUTINE ADD_INPLACE_NNZ_H_WHLDOM

! Helper: insert NEWCOL into sorted integer array JD and mirror REAL array D, grow allocs
SUBROUTINE INSERT_SORTED_INT_REAL(JD, D, NNZ, NEWCOL)
USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)
INTEGER, ALLOCATABLE, INTENT(INOUT) :: JD(:)
REAL(EB), ALLOCATABLE, INTENT(INOUT) :: D(:)
INTEGER, INTENT(INOUT) :: NNZ
INTEGER, INTENT(IN) :: NEWCOL

INTEGER :: POS, CAP
INTEGER :: I
INTEGER, ALLOCATABLE :: JD_TMP(:)
REAL(EB), ALLOCATABLE :: D_TMP(:)

! Determine current capacity from allocated size
CAP = SIZE(JD)

! Grow capacity geometrically if needed
IF (NNZ+1 > CAP) THEN
   CAP = MAX(2*MAX(CAP,1), NNZ+1)
   ALLOCATE(JD_TMP(CAP))
   ALLOCATE(D_TMP(CAP))
   IF (NNZ > 0) THEN
      JD_TMP(1:NNZ) = JD(1:NNZ)
      D_TMP(1:NNZ)  = D(1:NNZ)
   ENDIF
   IF (ALLOCATED(JD)) DEALLOCATE(JD)
   IF (ALLOCATED(D))  DEALLOCATE(D)
   CALL MOVE_ALLOC(JD_TMP, JD)
   CALL MOVE_ALLOC(D_TMP,  D)
ENDIF

! Find insertion position (ascending order)
POS = 1
DO WHILE (POS <= NNZ)
   IF(JD(POS) < NEWCOL) THEN
      POS = POS + 1
   ELSE
      EXIT
   ENDIF
ENDDO

! Shift in-place to make room
IF (NNZ >= POS) THEN
   DO I = NNZ, POS, -1
      JD(I+1) = JD(I)
      D(I+1)  = D(I)
   ENDDO
ENDIF

! Insert new entry
JD(POS) = NEWCOL
D(POS)  = 0._EB
NNZ = NNZ + 1

END SUBROUTINE INSERT_SORTED_INT_REAL


! --------------------------- GET_H_CUTFACES ------------------------------------

SUBROUTINE GET_H_CUTFACES(ONE_NM)

INTEGER, OPTIONAL, INTENT(IN) :: ONE_NM

! Local variables:
INTEGER :: NM,NM_LO,NM_HI
INTEGER :: NCELL,ICC,JCC,IFC,IFACE,LOWHIGH,ICF1,ICF2,IRC
INTEGER :: IW,II,JJ,KK,IIF,JJF,KKF,IOR,LOWHIGH_TEST,X1AXIS

IF (PRESENT(ONE_NM)) THEN
   NM_LO = ONE_NM
   NM_HI = ONE_NM
ELSE
   NM_LO = LOWER_MESH_INDEX
   NM_HI = UPPER_MESH_INDEX
ENDIF

! Mesh loop:
MAIN_MESH_LOOP : DO NM=NM_LO,NM_HI

   CALL POINT_TO_MESH(NM)

   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      NCELL = CUT_CELL(ICC)%NCELL
      DO JCC=1,NCELL
         DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
            IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
            ! If face type in face_list is not CC_FTYPE_CFGAS, drop:
            IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_CFGAS ) CYCLE

            ! Which face?
            LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
            ICF1    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
            ICF2    = CUT_CELL(ICC)%FACE_LIST(5,IFACE)

            IF ( LOWHIGH == LOW_IND) THEN ! Cut-face on low side of cut-cell:
               CUT_FACE(ICF1)%UNKH(HIGH_IND,ICF2) = CUT_CELL(ICC)%UNKH(JCC)
            ELSE ! HIGH
               CUT_FACE(ICF1)%UNKH(LOW_IND,ICF2) = CUT_CELL(ICC)%UNKH(JCC)
            ENDIF

         ENDDO
      ENDDO
   ENDDO

   ! Now Apply external wall cell loop for guard-cell cut cells:
   ULMAT_IF : IF (PRES_FLAG/=ULMAT_FLAG) THEN
      GUARD_CUT_CELL_LOOP :  DO IW=1,N_EXTERNAL_WALL_CELLS
         WC=>WALL(IW)
         IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE GUARD_CUT_CELL_LOOP

         BC => BOUNDARY_COORD(WC%BC_INDEX)
         II  = BC%II
         JJ  = BC%JJ
         KK  = BC%KK
         IOR = BC%IOR

         ! Drop if face is not of type CC_CUTCFE:
         X1AXIS=ABS(IOR)
         SELECT CASE(IOR)
         CASE( IAXIS)
            IIF=II  ; JJF=JJ  ; KKF=KK
            LOWHIGH_TEST=HIGH_IND ! Face on high side of Guard-Cell
         CASE(-IAXIS)
            IIF=II-1; JJF=JJ  ; KKF=KK
            LOWHIGH_TEST=LOW_IND
         CASE( JAXIS)
            IIF=II  ; JJF=JJ  ; KKF=KK
            LOWHIGH_TEST=HIGH_IND
         CASE(-JAXIS)
            IIF=II  ; JJF=JJ-1; KKF=KK
            LOWHIGH_TEST=LOW_IND
         CASE( KAXIS)
            IIF=II  ; JJF=JJ  ; KKF=KK
            LOWHIGH_TEST=HIGH_IND
         CASE(-KAXIS)
            IIF=II  ; JJF=JJ  ; KKF=KK-1
            LOWHIGH_TEST=LOW_IND
         END SELECT

         ! Copy CCVAR(II,JJ,KK,CC_CGSC) to guard cell:
         ICC = MESHES(NM)%CCVAR(II,JJ,KK,CC_IDCC)

         IF (FCVAR(IIF,JJF,KKF,CC_FGSC,X1AXIS) == CC_CUTCFE) THEN

         DO JCC=1,CUT_CELL(ICC)%NCELL
            ! Loop faces and test:
            DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
               IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
               ! Which face ?
               LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
               IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_CFGAS) CYCLE ! Must Be gasphase cut-face
               IF ( LOWHIGH                              /= LOWHIGH_TEST) CYCLE ! In same side as EWC from guard-cell
               IF ( CUT_CELL(ICC)%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE ! Normal to same axis as EWC
               ICF1    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
               ICF2    = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
               IF ( LOWHIGH == LOW_IND) THEN ! Cut-face on low side of cut-cell:
                  CUT_FACE(ICF1)%UNKH(HIGH_IND,ICF2) = CUT_CELL(ICC)%UNKH(JCC)
               ELSE ! HIGH
                  CUT_FACE(ICF1)%UNKH(LOW_IND,ICF2) = CUT_CELL(ICC)%UNKH(JCC)
               ENDIF
            ENDDO
         ENDDO

         ELSEIF (FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS) > 0) THEN ! RC_FACE
            IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
            IF(ICC>0) THEN
               DO JCC=1,CUT_CELL(ICC)%NCELL
                  ! Loop faces and test:
                  DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
                     IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
                     ! Which face ?
                     LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
                     IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE ! Must Be gasphase rc-face
                     IF ( LOWHIGH                          /= LOWHIGH_TEST) CYCLE ! In same side as EWC from guard-cell
                     IF ( CUT_CELL(ICC)%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE ! Normal to same axis as EWC
                     MESHES(NM)%RC_FACE(IRC)%UNKH(3-LOWHIGH_TEST)  = CUT_CELL(ICC)%UNKH(JCC)
                  ENDDO
               ENDDO
            ELSE
               MESHES(NM)%RC_FACE(IRC)%UNKH(3-LOWHIGH_TEST)  = CCVAR(II,JJ,KK,CC_UNKH)
            ENDIF
         ENDIF

      ENDDO GUARD_CUT_CELL_LOOP
   ENDIF ULMAT_IF

ENDDO MAIN_MESH_LOOP


RETURN
END SUBROUTINE GET_H_CUTFACES


! ---------------------------- GET_RCFACES_H ------------------------------------

SUBROUTINE GET_RCFACES_H(NM)

INTEGER, INTENT(IN) :: NM

! Local Variables:
INTEGER :: IRC,IIFC,X1AXIS,X2AXIS,X3AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: IJKFACE
INTEGER :: ICC,JCC,IJK(MAX_DIM),IFC,IFACE,LOWHIGH,IADD,JADD,KADD
INTEGER :: XIAXIS,XJAXIS,XKAXIS,INDXI1(MAX_DIM),INCELL,JNCELL,KNCELL,INFACE,JNFACE,KNFACE
INTEGER :: ISTR, IEND, JSTR, JEND, KSTR, KEND
LOGICAL :: INLIST

INTEGER :: IW,II,JJ,KK,IIF,JJF,KKF,IOR,LOWHIGH_TEST,CELL_UNKH
TYPE (WALL_TYPE), POINTER :: WC
TYPE (MESH_TYPE), POINTER :: M

! Test for Pressure Solver:
SELECT CASE(PRES_FLAG)
   CASE DEFAULT; RETURN ! no need to build unstructured matrix
   CASE (UGLMAT_FLAG,ULMAT_FLAG)
END SELECT

M => MESHES(NM)

! Mesh sizes:
NXB=M%IBAR; NYB=M%JBAR; NZB=M%KBAR

! X direction bounds:
ILO_FACE = 0                    ! Low mesh boundary face index.
IHI_FACE = M%IBAR               ! High mesh boundary face index.
ILO_CELL = ILO_FACE + 1         ! First internal cell index. See notes.
IHI_CELL = IHI_FACE             ! Last internal cell index.
ISTR     = ILO_FACE - NGUARD    ! Allocation start x arrays.
IEND     = IHI_FACE + NGUARD    ! Allocation end x arrays.

! Y direction bounds:
JLO_FACE = 0                    ! Low mesh boundary face index.
JHI_FACE = M%JBAR               ! High mesh boundary face index.
JLO_CELL = JLO_FACE + 1         ! First internal cell index. See notes.
JHI_CELL = JHI_FACE             ! Last internal cell index.
JSTR     = JLO_FACE - NGUARD    ! Allocation start y arrays.
JEND     = JHI_FACE + NGUARD    ! Allocation end y arrays.

! Z direction bounds:
KLO_FACE = 0                    ! Low mesh boundary face index.
KHI_FACE = M%KBAR               ! High mesh boundary face index.
KLO_CELL = KLO_FACE + 1         ! First internal cell index. See notes.
KHI_CELL = KHI_FACE             ! Last internal cell index.
KSTR     = KLO_FACE - NGUARD    ! Allocation start z arrays.
KEND     = KHI_FACE + NGUARD    ! Allocation end z arrays.

! First count for allocation:
ALLOCATE( IJKFACE(ILO_FACE:IHI_FACE,JLO_FACE:JHI_FACE,KLO_FACE:KHI_FACE,IAXIS:KAXIS) )
IJKFACE(:,:,:,:) = 0
DO ICC=1,M%N_CUTCELL_MESH
   CC => M%CUT_CELL(ICC); IJK(IAXIS:KAXIS) = CC%IJK(IAXIS:KAXIS)
   IF(CELL(CELL_INDEX(IJK(IAXIS),IJK(JAXIS),IJK(KAXIS)))%SOLID) CYCLE
   DO JCC=1,CC%NCELL
      ! Loop faces and test:
      DO IFC=1,CC%CCELEM(1,JCC)
         IFACE = CC%CCELEM(IFC+1,JCC)
         ! If face type in face_list is not CC_FTYPE_RCGAS, drop:
         IF(CC%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE
         ! Which face?
         LOWHIGH = CC%FACE_LIST(2,IFACE)
         X1AXIS  = CC%FACE_LIST(3,IFACE)
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            X2AXIS = JAXIS; X3AXIS = KAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
         CASE(JAXIS)
            X2AXIS = KAXIS; X3AXIS = IAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
         CASE(KAXIS)
            X2AXIS = IAXIS; X3AXIS = JAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
         END SELECT

         IF (LOWHIGH == LOW_IND) THEN
            ! Face indexes:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
            INFACE = INDXI1(XIAXIS); JNFACE = INDXI1(XJAXIS); KNFACE = INDXI1(XKAXIS)

            ! Location of next Cartesian cell:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
            INCELL = INDXI1(XIAXIS); JNCELL = INDXI1(XJAXIS); KNCELL = INDXI1(XKAXIS)

            IF (PRES_FLAG==UGLMAT_FLAG) CELL_UNKH = M%CCVAR(INCELL,JNCELL,KNCELL,CC_UNKH)
            IF (PRES_FLAG== ULMAT_FLAG) CELL_UNKH = M%MUNKH(INCELL,JNCELL,KNCELL)
            IF ( CELL_UNKH > 0 ) THEN
               IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) = 1
            ELSEIF ( M%CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE ) THEN ! Cut-cell.
               IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) = 1
            ENDIF
         ELSE ! HIGH_IND
            ! Face indexes:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS), IJK(X2AXIS), IJK(X3AXIS) /)
            INFACE = INDXI1(XIAXIS); JNFACE = INDXI1(XJAXIS); KNFACE = INDXI1(XKAXIS)

            ! Location of next Cartesian cell:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)+1, IJK(X2AXIS), IJK(X3AXIS) /)
            INCELL = INDXI1(XIAXIS); JNCELL = INDXI1(XJAXIS); KNCELL = INDXI1(XKAXIS)

            IF (PRES_FLAG==UGLMAT_FLAG) CELL_UNKH = M%CCVAR(INCELL,JNCELL,KNCELL,CC_UNKH)
            IF (PRES_FLAG== ULMAT_FLAG) CELL_UNKH = M%MUNKH(INCELL,JNCELL,KNCELL)
            IF ( CELL_UNKH > 0 ) THEN
               IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) = 1
            ELSEIF ( M%CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE ) THEN ! Cut-cell.
               IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) = 1
            ENDIF
         ENDIF

      ENDDO
   ENDDO
ENDDO

! Check for RCF_H on the boundary of the domain, where the cut-cell is in the cut-cell region.
! Now Apply external wall cell loop for guard-cell cut cells:
GRD_CC_LOOP_1 :  DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
   WC=>M%WALL(IW); BC => M%BOUNDARY_COORD(WC%BC_INDEX); II = BC%II; JJ = BC%JJ; KK = BC%KK; IOR = BC%IOR
   ! Which face:
   X1AXIS=ABS(IOR)
   IADD = 0; JADD = 0; KADD = 0
   SELECT CASE(IOR)
   CASE( IAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; IADD=-1
   CASE(-IAXIS); IIF=II-1; JJF=JJ  ; KKF=KK  ; IADD=-1
   CASE( JAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; JADD=-1
   CASE(-JAXIS); IIF=II  ; JJF=JJ-1; KKF=KK  ; JADD=-1
   CASE( KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; KADD=-1
   CASE(-KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK-1; KADD=-1
   END SELECT
   IF(ALL(CCVAR(IIF+IADD:IIF,JJF+JADD:JJF,KKF+KADD:KKF,CC_CGSC)==CC_SOLID)) CYCLE GRD_CC_LOOP_1
   IF( PRES_FLAG==UGLMAT_FLAG .AND. &
      (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY) ) THEN
      ! Drop if FACE is not type CC_GASPHASE
      IF (M%FCVAR(IIF,JJF,KKF,CC_FGSC,X1AXIS) /= CC_GASPHASE) CYCLE GRD_CC_LOOP_1
      ! Is this an actual RCF_H laying on the mesh boundary, where the cut-cell is in the guard-cell region?
      IF(.NOT.((M%CCVAR(II,JJ,KK,CC_CGSC)==CC_CUTCFE).AND.(M%CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC)==CC_GASPHASE))) &
      CYCLE GRD_CC_LOOP_1
      IJKFACE(IIF,JJF,KKF,X1AXIS) = 1

   ELSE ! All other types of BCs (SOLID_BOUNDARY, NULL_BOUNDARY, OPEN_BOUNDARY) will not be added to RCF_H.
      IJKFACE(IIF,JJF,KKF,X1AXIS) = 0

   ENDIF
ENDDO GRD_CC_LOOP_1

IRC = SUM(IJKFACE(:,:,:,:))
IF (IRC == 0) THEN
   DEALLOCATE(IJKFACE)
   RETURN
ELSE
   ! Compute xc, yc, zc:
   ! Populate position and cell size arrays: Uniform grid implementation.
   ! X direction:
   ALLOCATE(XCELL(ISTR:IEND));  XCELL = 1._EB/GEOMEPS ! Initialize huge.
   XCELL(ILO_CELL-1:IHI_CELL+1) = M%XC(ILO_CELL-1:IHI_CELL+1)

   ! Y direction:
   ALLOCATE(YCELL(JSTR:JEND));  YCELL = 1._EB/GEOMEPS ! Initialize huge.
   YCELL(JLO_CELL-1:JHI_CELL+1) = M%YC(JLO_CELL-1:JHI_CELL+1)

   ! Z direction:
   ALLOCATE(ZCELL(KSTR:KEND));  ZCELL = 1._EB/GEOMEPS ! Initialize huge.
   ZCELL(KLO_CELL-1:KHI_CELL+1) = M%ZC(KLO_CELL-1:KHI_CELL+1)
ENDIF

M%CC_NRCFACE_H = IRC ! Same number of regular - cut cell faces as scalars.
IF (ALLOCATED(M%RCF_H)) DEALLOCATE(M%RCF_H); ALLOCATE( M%RCF_H(IRC) )
IRC = 0
DO ICC=1,M%N_CUTCELL_MESH
   CC => M%CUT_CELL(ICC); IJK(IAXIS:KAXIS) = CC%IJK(IAXIS:KAXIS)
   IF(CELL(CELL_INDEX(IJK(IAXIS),IJK(JAXIS),IJK(KAXIS)))%SOLID) CYCLE
   DO JCC=1,CC%NCELL
      ! Loop faces and test:
      DO IFC=1,CC%CCELEM(1,JCC)
         IFACE = CC%CCELEM(IFC+1,JCC)
         ! If face type in face_list is not CC_FTYPE_RCGAS, drop:
         IF(CC%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE
         ! Which face?
         LOWHIGH = CC%FACE_LIST(2,IFACE)
         X1AXIS  = CC%FACE_LIST(3,IFACE)
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            X2AXIS = JAXIS; X3AXIS = KAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
         CASE(JAXIS)
            X2AXIS = KAXIS; X3AXIS = IAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
         CASE(KAXIS)
            X2AXIS = IAXIS; X3AXIS = JAXIS
            ! location in I,J,K od x2,x2,x3 axes:
            XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
         END SELECT

         IF_LOW_HIGH_H : IF (LOWHIGH == LOW_IND) THEN

            ! Face indexes:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
            INFACE = INDXI1(XIAXIS); JNFACE = INDXI1(XJAXIS); KNFACE = INDXI1(XKAXIS)

            IF(IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) /= 1) CYCLE ! This is to cycle external WALL CELLs of types other
                                                                ! than INTERPOLATED_BOUNDARY or PERIODIC_BOUNDARY.

            ! Location of next Cartesian cell:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
            INCELL = INDXI1(XIAXIS); JNCELL = INDXI1(XJAXIS); KNCELL = INDXI1(XKAXIS)

            IF (PRES_FLAG==UGLMAT_FLAG) CELL_UNKH = M%CCVAR(INCELL,JNCELL,KNCELL,CC_UNKH)
            IF (PRES_FLAG== ULMAT_FLAG) CELL_UNKH = M%MUNKH(INCELL,JNCELL,KNCELL)
            IF ( CELL_UNKH > 0 ) THEN
               ! Add face to RCF_H list:
               IRC = IRC + 1
               M%RCF_H(IRC) = FCVAR(INFACE,JNFACE,KNFACE,CC_IDRC,X1AXIS)
               ! Cell at i-1, i.e. regular GASPHASE:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(LOW_IND) = CELL_UNKH
               ! Cell at i+1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(HIGH_IND) = CC%UNKH(JCC)

            ELSEIF ( M%CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE ) THEN ! Cut-cell.
               ! Test that Cut-cell to Cut-cell reg face hasn't been added before:
               INLIST = .FALSE.
               DO IIFC=1,IRC
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(IAXIS)   /= INFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(JAXIS)   /= JNFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(KAXIS)   /= KNFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(KAXIS+1) /= X1AXIS ) CYCLE
                  INLIST = .TRUE.
                  EXIT
               ENDDO
               IF (INLIST) THEN
                  ! Cell at i+1, i.e. cut-cell:
                  M%RC_FACE(M%RCF_H(IIFC))%UNKH(HIGH_IND) = CC%UNKH(JCC)
                  CYCLE
               ENDIF
               ! Add face to RCF_H list:
               IRC = IRC + 1
               M%RCF_H(IRC) = FCVAR(INFACE,JNFACE,KNFACE,CC_IDRC,X1AXIS)
               ! Cell at i+1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(HIGH_IND) = CC%UNKH(JCC)

            ENDIF

         ELSE ! IF_LOW_HIGH_H : HIGH_IND

            ! Face indexes:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS), IJK(X2AXIS), IJK(X3AXIS) /)
            INFACE = INDXI1(XIAXIS); JNFACE = INDXI1(XJAXIS); KNFACE = INDXI1(XKAXIS)

            IF(IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS) /= 1) CYCLE ! This is to cycle external WALL CELLs of types other
                                                                ! than INTERPOLATED_BOUNDARY of PERIODIC_BOUNDARY.

            ! Location of next Cartesian cell:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)+1, IJK(X2AXIS), IJK(X3AXIS) /)
            INCELL = INDXI1(XIAXIS); JNCELL = INDXI1(XJAXIS); KNCELL = INDXI1(XKAXIS)

            IF (PRES_FLAG==UGLMAT_FLAG) CELL_UNKH = M%CCVAR(INCELL,JNCELL,KNCELL,CC_UNKH)
            IF (PRES_FLAG== ULMAT_FLAG) CELL_UNKH = M%MUNKH(INCELL,JNCELL,KNCELL)
            IF ( CELL_UNKH > 0 ) THEN
               ! Add face to RCF_H list:
               IRC = IRC + 1
               M%RCF_H(IRC) = FCVAR(INFACE,JNFACE,KNFACE,CC_IDRC,X1AXIS)
               ! Cell at i-1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(LOW_IND) = CC%UNKH(JCC)
               ! Cell at i+1, i.e. regular GASPHASE:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(HIGH_IND) = CELL_UNKH

            ELSEIF ( M%CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE ) THEN ! Cut-cell.
               ! Test that Cut-cell to Cut-cell reg face hasn't been added before:
               INLIST = .FALSE.
               DO IIFC=1,IRC
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(IAXIS)   /= INFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(JAXIS)   /= JNFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(KAXIS)   /= KNFACE ) CYCLE
                  IF ( M%RC_FACE(M%RCF_H(IIFC))%IJK(KAXIS+1) /= X1AXIS ) CYCLE
                  INLIST = .TRUE.
                  EXIT
               ENDDO
               IF (INLIST) THEN
                  ! Cell at i-1, i.e. cut-cell:
                  M%RC_FACE(M%RCF_H(IIFC))%UNKH(LOW_IND) = CC%UNKH(JCC)
                  CYCLE
               ENDIF
               ! Add face to RCF_H list:
               IRC = IRC + 1
               M%RCF_H(IRC) = FCVAR(INFACE,JNFACE,KNFACE,CC_IDRC,X1AXIS)
               ! Cell at i-1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(LOW_IND) = CC%UNKH(JCC)

            ENDIF
         ENDIF IF_LOW_HIGH_H

      ENDDO
   ENDDO
ENDDO

IF (PRES_FLAG==UGLMAT_FLAG) THEN
   GRD_CC_LOOP_2 :  DO IW=1,M%N_EXTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF(.NOT.(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY.OR.WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE GRD_CC_LOOP_2
      BC => M%BOUNDARY_COORD(WC%BC_INDEX); II = BC%II; JJ = BC%JJ; KK = BC%KK; IOR = BC%IOR

      ! Which face:
      X1AXIS=ABS(IOR)
      SELECT CASE(IOR)
      CASE( IAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND ! Face on high side of Guard-Cell
      CASE(-IAXIS)
         IIF=II-1; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST= LOW_IND
      CASE( JAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND
      CASE(-JAXIS)
         IIF=II  ; JJF=JJ-1; KKF=KK  ; LOWHIGH_TEST= LOW_IND
      CASE( KAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND
      CASE(-KAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK-1; LOWHIGH_TEST= LOW_IND
      END SELECT

      ! Drop if FACE is not type CC_GASPHASE
      IF (M%FCVAR(IIF,JJF,KKF,CC_FGSC,X1AXIS) /= CC_GASPHASE) CYCLE GRD_CC_LOOP_2

      ! Is this an actual RCF_H laying on the mesh boundary, where the cut-cell is in the guard-cell region?
      IF(.NOT.((M%CCVAR(II,JJ,KK,CC_CGSC)==CC_CUTCFE).AND.(M%CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_UNKH)>0))) &
      CYCLE GRD_CC_LOOP_2

      CC => M%CUT_CELL(M%CCVAR(II,JJ,KK,CC_IDCC))
      DO JCC=1,CC%NCELL
         ! Loop faces and test:
         DO IFC=1,CC%CCELEM(1,JCC)
            IFACE = CC%CCELEM(IFC+1,JCC)
            ! Which face ?
            LOWHIGH = CC%FACE_LIST(2,IFACE)
            IF ( CC%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE  ! Must Be gasphase RCFACE
            IF ( LOWHIGH                              /= LOWHIGH_TEST) CYCLE ! In same side as EWC from guard-cell
            IF ( CC%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE           ! Normal to same axis as EWC

             ! If so, we need to add it to RCF_H list:
            IF (LOWHIGH == LOW_IND) THEN ! Face on low side of guard cut-cell
               IRC = IRC + 1
               M%RCF_H(IRC) = M%FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               ! Cell at i-1, i.e. regular GASPHASE:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(LOW_IND)  = M%CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_UNKH)
               ! Cell at i+1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(HIGH_IND) = CC%UNKH(JCC)

            ELSEIF(LOWHIGH == HIGH_IND) THEN ! Face on high side of guard cut-cell
               IRC = IRC + 1
               M%RCF_H(IRC) = M%FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               ! Cell at i-1, i.e. cut-cell:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(LOW_IND)  = CC%UNKH(JCC)
               ! Cell at i+1, i.e. regular GASPHASE:
               M%RC_FACE(M%RCF_H(IRC))%UNKH(HIGH_IND) = M%CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_UNKH)

            ENDIF
            ! At this point the face has been found, cycle:
            CYCLE GRD_CC_LOOP_2
         ENDDO
      ENDDO

   ENDDO GRD_CC_LOOP_2
ENDIF

DEALLOCATE(XCELL,YCELL,ZCELL)
DEALLOCATE(IJKFACE)

RETURN
END SUBROUTINE GET_RCFACES_H


! ------------------ NUMBER_UNKH_CUTCELLS ---------------------------

SUBROUTINE NUMBER_UNKH_CUTCELLS(FLAG12,NM,IPZ,NUNKH_LC)

LOGICAL, INTENT(IN) :: FLAG12
INTEGER, INTENT(IN) :: NM,IPZ
INTEGER, INTENT(INOUT) :: NUNKH_LC(LOWER_MESH_INDEX:UPPER_MESH_INDEX)

! Local Variables:
INTEGER :: ICC, JCC, I, J, K

FLAG12_COND : IF (FLAG12) THEN
   ! Initialize Cut-cell unknown numbers as undefined.
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CUT_CELL(ICC)%UNKH(:) = CC_UNDEFINED
   ENDDO
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ ) CYCLE
      NUNKH_LC(NM) = NUNKH_LC(NM) + 1
      CUT_CELL(ICC)%UNKH(1:CC%NCELL) = NUNKH_LC(NM)
   ENDDO
ELSE
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IF(CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      IF(ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ ) CYCLE
      DO JCC=1,CC%NCELL
         CUT_CELL(ICC)%UNKH(JCC) = CUT_CELL(ICC)%UNKH(JCC) + ZONE_SOLVE(IPZ)%UNKH_IND(NM)
      ENDDO
   ENDDO
ENDIF FLAG12_COND

RETURN
END SUBROUTINE NUMBER_UNKH_CUTCELLS

! ------------------- COPY_CC_MUNKH_TO_UNKH ----------------------------

SUBROUTINE COPY_CC_MUNKH_TO_UNKH

! Local Variables:
INTEGER :: NOM,ICC,IW,IIO,JJO,KKO,II,JJ,KK
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
! Loop over external wall cells:
EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   NOM = EWC%NOM
   OM => OMESH(NOM)

   BC => BOUNDARY_COORD(WC%BC_INDEX)
   II = BC%II; JJ = BC%JJ; KK = BC%KK
   ! Copy MUNKH to current mesh guard cell storage (for same-level CF%UNKH access)
   ICC = CCVAR(II,JJ,KK,CC_IDCC)
   IF (ICC > 0) THEN ! Cut-cells on this guard-cell Cartesian cell.
      CUT_CELL(ICC)%UNKH(1) = OM%MUNKH(EWC%IIO_MIN,EWC%JJO_MIN,EWC%KKO_MIN)
   ELSE
      CCVAR(II,JJ,KK,CC_UNKH) = OM%MUNKH(EWC%IIO_MIN,EWC%JJO_MIN,EWC%KKO_MIN)
   ENDIF

   ! Loop over all cells in mesh NOM that correspond to this boundary face.
   ! This populates MESHES(NOM)%CUT_CELL%UNKH for access during matrix assembly.
   ! (Supports grid refinement where multiple external cells map to one guard cell)
   DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
      DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
         DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
            ! Check if this cell in mesh NOM is a cut-cell:
            ICC = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
            IF (ICC > 0) THEN
               ! Copy to ghost cut-cell storage for mesh NOM:
               MESHES(NOM)%CUT_CELL(ICC)%UNKH(1) = OM%MUNKH(IIO,JJO,KKO)
            ELSE
               ! Copy to ghost regular cell storage for mesh NOM:
               MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_UNKH) = OM%MUNKH(IIO,JJO,KKO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDDO EXTERNAL_WALL_LOOP

! Loop over external wall cells:
EXTERNAL_WALL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS
   WC=>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP2
   NOM = EWC%NOM; OM => OMESH(NOM)
   DO KKO = EWC%KKO_MIN, EWC%KKO_MAX
      DO JJO = EWC%JJO_MIN, EWC%JJO_MAX
         DO IIO = EWC%IIO_MIN, EWC%IIO_MAX
            OM%HS(IIO,JJO,KKO) = 0._EB ! (VAR_CC == UNKH)
         ENDDO
      ENDDO
   ENDDO
ENDDO EXTERNAL_WALL_LOOP2

RETURN
END SUBROUTINE COPY_CC_MUNKH_TO_UNKH

! ------------------- COPY_CC_UNKH_TO_HS ----------------------------

SUBROUTINE COPY_CC_UNKH_TO_HS(NM)

INTEGER, INTENT(IN) :: NM

! Local Variables:
INTEGER :: I,J,K,ICC

DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   I = MESHES(NM)%CUT_CELL(ICC)%IJK(IAXIS)
   J = MESHES(NM)%CUT_CELL(ICC)%IJK(JAXIS)
   K = MESHES(NM)%CUT_CELL(ICC)%IJK(KAXIS)
   HS(I,J,K)= REAL(MESHES(NM)%CUT_CELL(ICC)%UNKH(1),EB)
ENDDO

RETURN
END SUBROUTINE COPY_CC_UNKH_TO_HS

END MODULE CC_PRESSURE
