!  +++++++++++++++++++++++ CC_VELOCITY_KERNELS ++++++++++++++++++++++++++

! Thread-safe computation kernels for CC velocity routines.
! Each routine takes TYPE(MESH_TYPE) as explicit argument
! instead of relying on MESH_POINTERS.

MODULE CC_VELOCITY_KERNELS

USE PRECISION_PARAMETERS
USE TYPES
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE COMPLEX_GEOMETRY, ONLY: CC_FTYPE_RCGAS, CC_FTYPE_CFGAS, &
                            CC_FTYPE_CFINB, CC_FGSC, CC_SOLID, &
                            CC_GASPHASE, CC_INBOUNDARY, CC_CGSC, &
                            CC_IDCF, CC_IDRC, CC_UNKZ, NM_START, &
                            CC_VELOCITY_FLUX_TIME_INDEX, &
                            CC_COMPUTE_VISCOSITY_TIME_INDEX, &
                            T_CC_USED
USE CC_SCALARS_DATA, ONLY: TIME_CC_IBM, &
                           UNKZ_IND, RHO_0_CV
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CUTFACE_VELOCITIES, CC_CUTCELL_VELOCITY, CC_COMPUTE_KRES, &
          CC_COMPUTE_VISCOSITY, CC_STORE_FACE_FV, CC_VELOCITY_FLUX, &
          CC_PROJECT_VELOCITY_KERNEL

CONTAINS


! ------------------------------ CUTFACE_VELOCITIES --------------------------------

SUBROUTINE CUTFACE_VELOCITIES(M,UU,VV,WW,CUTFACES)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), POINTER, DIMENSION(:,:,:), INTENT(INOUT) :: UU,VV,WW
LOGICAL, INTENT(IN) :: CUTFACES
TYPE(CC_CUTFACE_TYPE), POINTER :: CF
TYPE(CFACE_TYPE), POINTER :: CFA
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
INTEGER :: ICF,ICFA,JCF,I,J,K,X1AXIS
REAL(EB):: AREA,VELN(3),PREDFCT


CUTFACES_IF : IF (CUTFACES) THEN ! USE CUT_FACE(ICF)%VEL_CF
   DO ICF=1,M%N_CUTFACE_MESH+M%N_GCCUTFACE_MESH
      CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
      I      = CF%IJK(IAXIS); IF(I<0 .OR. I>M%IBP1) CYCLE
      J      = CF%IJK(JAXIS); IF(J<0 .OR. J>M%JBP1) CYCLE
      K      = CF%IJK(KAXIS); IF(K<0 .OR. K>M%KBP1) CYCLE
      X1AXIS = CF%IJK(KAXIS+1)
      SELECT CASE(X1AXIS)
      CASE(IAXIS); UU(I,J,K) = CF%VEL_CF
      CASE(JAXIS); VV(I,J,K) = CF%VEL_CF
      CASE(KAXIS); WW(I,J,K) = CF%VEL_CF
      END SELECT
   ENDDO

   PREDFCT = 0._EB; IF(PREDICTOR) PREDFCT = 1._EB
   ! CFACEs, set velocity in underlaying solid cartesian faces to be used in VELOCITY_FLUX:
   DO ICF=1,M%N_CUTFACE_MESH
      CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_INBOUNDARY) CYCLE
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS)
      ! Area Average velocity for boundary CFACEs:
      AREA = 0._EB; VELN(IAXIS:KAXIS) = 0._EB
      DO JCF=1,CF%NFACE
         ICFA=CF%CFACE_INDEX(JCF); IF(ICFA<1) CYCLE
         CFA => M%CFACE(ICFA)
         BC => M%BOUNDARY_COORD(CFA%BC_INDEX)
         B1 => M%BOUNDARY_PROP1(CFA%B1_INDEX)
         AREA = AREA+CFA%AREA
         VELN(IAXIS:KAXIS) = VELN(IAXIS:KAXIS) - &
            (PREDFCT*B1%U_NORMAL + &
            (1._EB-PREDFCT)*B1%U_NORMAL_S) * &
            CFA%AREA*BC%NVEC(IAXIS:KAXIS)
      ENDDO
      VELN(IAXIS:KAXIS) = VELN(IAXIS:KAXIS)/(AREA+TWENTY_EPSILON_EB)
      ! Distribute into Solid cartesian faces when SOLID cell is present behind:
      IF(M%FCVAR(I-1,J,K,CC_FGSC,IAXIS)==CC_SOLID .AND. &
         M%CCVAR(I-1,J,K,CC_CGSC)==CC_SOLID) &
         UU(I-1,J,K) = VELN(IAXIS)
      IF(M%FCVAR(I  ,J,K,CC_FGSC,IAXIS)==CC_SOLID .AND. &
         M%CCVAR(I+1,J,K,CC_CGSC)==CC_SOLID) &
         UU(I  ,J,K) = VELN(IAXIS)
      IF(M%FCVAR(I,J-1,K,CC_FGSC,JAXIS)==CC_SOLID .AND. &
         M%CCVAR(I,J-1,K,CC_CGSC)==CC_SOLID) &
         VV(I,J-1,K) = VELN(JAXIS)
      IF(M%FCVAR(I,J  ,K,CC_FGSC,JAXIS)==CC_SOLID .AND. &
         M%CCVAR(I,J+1,K,CC_CGSC)==CC_SOLID) &
         VV(I,J  ,K) = VELN(JAXIS)
      IF(M%FCVAR(I,J,K-1,CC_FGSC,KAXIS)==CC_SOLID .AND. &
         M%CCVAR(I,J,K-1,CC_CGSC)==CC_SOLID) &
         WW(I,J,K-1) = VELN(KAXIS)
      IF(M%FCVAR(I,J,K  ,CC_FGSC,KAXIS)==CC_SOLID .AND. &
         M%CCVAR(I,J,K+1,CC_CGSC)==CC_SOLID) &
         WW(I,J,K  ) = VELN(KAXIS)
   ENDDO

ELSE CUTFACES_IF ! USE CUT_FACE(ICF)%VEL_CRT
   DO ICF=1,M%N_CUTFACE_MESH+M%N_GCCUTFACE_MESH
      CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
      I      = CF%IJK(IAXIS); IF(I<0 .OR. I>M%IBP1) CYCLE
      J      = CF%IJK(JAXIS); IF(J<0 .OR. J>M%JBP1) CYCLE
      K      = CF%IJK(KAXIS); IF(K<0 .OR. K>M%KBP1) CYCLE
      X1AXIS = CF%IJK(KAXIS+1)
      SELECT CASE(X1AXIS)
      CASE(IAXIS); UU(I,J,K) = CF%VEL_CRT
      CASE(JAXIS); VV(I,J,K) = CF%VEL_CRT
      CASE(KAXIS); WW(I,J,K) = CF%VEL_CRT
      END SELECT
   ENDDO

   DO ICF=1,M%N_CUTFACE_MESH
      CF => M%CUT_FACE(ICF)
      IF(CF%STATUS/=CC_INBOUNDARY) CYCLE
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS)
      WHERE(M%FCVAR(I-1:I,J,K,CC_FGSC,IAXIS)==CC_SOLID) &
         UU(I-1:I,J,K) = 0._EB
      WHERE(M%FCVAR(I,J-1:J,K,CC_FGSC,JAXIS)==CC_SOLID) &
         VV(I,J-1:J,K) = 0._EB
      WHERE(M%FCVAR(I,J,K-1:K,CC_FGSC,KAXIS)==CC_SOLID) &
         WW(I,J,K-1:K) = 0._EB
   ENDDO

ENDIF CUTFACES_IF

RETURN
END SUBROUTINE CUTFACE_VELOCITIES


! ------------------------------- CC_CUTCELL_VELOCITY -------------------------------

SUBROUTINE CC_CUTCELL_VELOCITY(M,PRFCT,ICC,JCC, &
                               UVWAV,ATOTV,RETURN_INTEGRALS)

! Routine computes area averaged velocity vector for a cut-cell ICC, JCC.
! Areas used for each component are the projected areas in the component direction.

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER,  INTENT(IN) :: ICC,JCC
REAL(EB), INTENT(IN) :: PRFCT
REAL(EB), INTENT(OUT):: UVWAV(MAX_DIM),ATOTV(MAX_DIM)
LOGICAL,  INTENT(IN) :: RETURN_INTEGRALS

! Local Variables:
TYPE(CC_CUTCELL_TYPE), POINTER :: CC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER :: I,J,K,IFC,IFACE,LOWHIGH,X1AXIS,ILH,IFC2,IFACE2,ICFA
REAL(EB):: AUI,AF,VELN

ATOTV(:) = 0._EB; UVWAV(:) = 0._EB
CC => M%CUT_CELL(ICC)
I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)

IFC_LOOP : DO IFC=1,CC%CCELEM(1,JCC)
   IFACE = CC%CCELEM(IFC+1,JCC)
   SELECT CASE(CC%FACE_LIST(1,IFACE))
   CASE(CC_FTYPE_RCGAS) ! REGULAR GASPHASE
      LOWHIGH = CC%FACE_LIST(2,IFACE)
      X1AXIS  = CC%FACE_LIST(3,IFACE)
      ILH     = LOWHIGH - 1
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF   = M%DY(J)*M%DZ(K)
         VELN = PRFCT*M%US(I-1+ILH,J,K) + &
                (1._EB-PRFCT)*M%U(I-1+ILH,J,K)
      CASE(JAXIS)
         AF   = M%DX(I)*M%DZ(K)
         VELN = PRFCT*M%VS(I,J-1+ILH,K) + &
                (1._EB-PRFCT)*M%V(I,J-1+ILH,K)
      CASE(KAXIS)
         AF   = M%DX(I)*M%DY(J)
         VELN = PRFCT*M%WS(I,J,K-1+ILH) + &
                (1._EB-PRFCT)*M%W(I,J,K-1+ILH)
      END SELECT
      ATOTV(X1AXIS) = ATOTV(X1AXIS) + AF
      UVWAV(X1AXIS) = UVWAV(X1AXIS) + VELN * AF
   CASE(CC_FTYPE_CFGAS) ! GASPHASE CUT FACE:
      LOWHIGH = CC%FACE_LIST(2,IFACE)
      IFC2    = CC%FACE_LIST(4,IFACE)
      IFACE2  = CC%FACE_LIST(5,IFACE)
      X1AXIS  = M%CUT_FACE(IFC2)%IJK(KAXIS+1)
      AF      = M%CUT_FACE(IFC2)%AREA(IFACE2)
      VELN    = PRFCT*M%CUT_FACE(IFC2)%VELS(IFACE2) + &
                (1._EB-PRFCT)*M%CUT_FACE(IFC2)%VEL(IFACE2)
      ATOTV(X1AXIS) = ATOTV(X1AXIS) + AF
      UVWAV(X1AXIS) = UVWAV(X1AXIS) + VELN * AF
   CASE(CC_FTYPE_CFINB) ! INBOUNDARY CUT FACE
      IFC2    = CC%FACE_LIST(4,IFACE)
      IFACE2  = CC%FACE_LIST(5,IFACE)
      ICFA    = M%CUT_FACE(IFC2)%CFACE_INDEX(IFACE2)
      AF      = M%CUT_FACE(IFC2)%AREA(IFACE2)
      VELN = (1._EB-PRFCT)*M%CUT_FACE(IFC2)%VEL(IFACE2) &
           + PRFCT*M%CUT_FACE(IFC2)%VELS(IFACE2)
      ! - to use velocity into gasphase, projected area.
      BC => M%BOUNDARY_COORD(M%CFACE(ICFA)%BC_INDEX)
      DO X1AXIS=IAXIS,KAXIS
         AUI  = ABS(BC%NVEC(X1AXIS)*AF)
         ATOTV(X1AXIS) = ATOTV(X1AXIS) + AUI
         UVWAV(X1AXIS) = UVWAV(X1AXIS) - &
                         VELN*BC%NVEC(X1AXIS)*AUI
      ENDDO
   END SELECT
ENDDO IFC_LOOP
IF(.NOT.RETURN_INTEGRALS) &
   WHERE(ATOTV>TWENTY_EPSILON_EB) UVWAV = UVWAV / ATOTV

END SUBROUTINE CC_CUTCELL_VELOCITY


! ------------------------------- CC_COMPUTE_KRES -----------------------------

SUBROUTINE CC_COMPUTE_KRES(M,APPLY_TO_ESTIMATED_VARIABLES)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES

! Local Variables:
TYPE(CC_CUTCELL_TYPE), POINTER :: CC
INTEGER :: ICC,I,J,K,JCC
REAL(EB):: T_NOW,PRFCT,UVWAV(MAX_DIM),ATOTV(MAX_DIM), &
           UVWAV_AUX(MAX_DIM),ATOTV_AUX(MAX_DIM)

T_NOW = CURRENT_TIME()

PRFCT = 0._EB; IF(APPLY_TO_ESTIMATED_VARIABLES) PRFCT = 1._EB
CUTCELL_DO : DO ICC=1,M%N_CUTCELL_MESH
   CC => M%CUT_CELL(ICC)
   I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
   IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
   IF (ONE_UNKH_PER_CUTCELL) THEN
      DO JCC=1,CC%NCELL
         CALL CC_CUTCELL_VELOCITY(M,PRFCT,ICC,JCC, &
            UVWAV,ATOTV,RETURN_INTEGRALS=.FALSE.)
         CC%KRES(JCC) = 0.5_EB * &
            (UVWAV(IAXIS)**2._EB + &
             UVWAV(JAXIS)**2._EB + &
             UVWAV(KAXIS)**2._EB)
      ENDDO

   ELSE
      ATOTV(:) = 0._EB; UVWAV(:) = 0._EB
      DO JCC=1,CC%NCELL
         CALL CC_CUTCELL_VELOCITY(M,PRFCT,ICC,JCC, &
            UVWAV_AUX,ATOTV_AUX,RETURN_INTEGRALS=.TRUE.)
         ! Accumulate individual cut-cell averages for this ICC.
         UVWAV = UVWAV + UVWAV_AUX
         ATOTV = ATOTV + ATOTV_AUX
      ENDDO
      WHERE(ATOTV>TWENTY_EPSILON_EB) UVWAV = UVWAV / ATOTV
      CC%KRES(1:CC%NCELL) = 0.5_EB * &
         (UVWAV(IAXIS)**2._EB + &
          UVWAV(JAXIS)**2._EB + &
          UVWAV(KAXIS)**2._EB)
   ENDIF
   ! Note we use an average KRES per cartesian cell.
   M%KRES(I,J,K) = DOT_PRODUCT( &
      CC%KRES(1:CC%NCELL), &
      CC%VOLUME(1:CC%NCELL)) / &
      SUM(CC%VOLUME(1:CC%NCELL))
ENDDO CUTCELL_DO

T_USED(14) = T_USED(14) + CURRENT_TIME() - T_NOW
RETURN
END SUBROUTINE CC_COMPUTE_KRES


! ----------------------------- CC_COMPUTE_VISCOSITY -------------------------

SUBROUTINE CC_COMPUTE_VISCOSITY(M,DT)

USE TURB_KERNELS, ONLY: WALE_VISCOSITY

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN):: DT

! Local Variables:
INTEGER :: I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: NU_EDDY,DELTA,A_IJ(3,3), &
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ

REAL(EB) :: TNOW

! Dummy assignments:
TNOW = DT
I    = 0

TNOW = CURRENT_TIME()

IF (PREDICTOR) THEN
   RHOP => M%RHO
   UU   => M%U
   VV   => M%V
   WW   => M%W
   ZZP  => M%ZZ
ELSE
   RHOP => M%RHOS
   UU   => M%US
   VV   => M%VS
   WW   => M%WS
   ZZP  => M%ZZS
ENDIF

! No need to compute WALE model turbulent viscosity on cut-cell region.
LES_IF : IF (SIM_MODE/=DNS_MODE) THEN

   ! WALE model on cells belonging to cut-cell region:
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
             IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
             IF (M%CCVAR(I,J,K,CC_CGSC)==CC_SOLID) THEN
                M%MU(I,J,K) = M%MU_DNS(I,J,K); CYCLE
             ENDIF
             IF (M%CCVAR(I,J,K,CC_IDCF)<1) CYCLE

             DELTA = M%LES_FILTER_WIDTH(I,J,K)
             ! compute velocity gradient tensor
             DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
             DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
             DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
             DUDY = 0.25_EB*M%RDY(J) * &
                (UU(I,J+1,K)-UU(I,J-1,K) + &
                 UU(I-1,J+1,K)-UU(I-1,J-1,K))
             DUDZ = 0.25_EB*M%RDZ(K) * &
                (UU(I,J,K+1)-UU(I,J,K-1) + &
                 UU(I-1,J,K+1)-UU(I-1,J,K-1))
             DVDX = 0.25_EB*M%RDX(I) * &
                (VV(I+1,J,K)-VV(I-1,J,K) + &
                 VV(I+1,J-1,K)-VV(I-1,J-1,K))
             DVDZ = 0.25_EB*M%RDZ(K) * &
                (VV(I,J,K+1)-VV(I,J,K-1) + &
                 VV(I,J-1,K+1)-VV(I,J-1,K-1))
             DWDX = 0.25_EB*M%RDX(I) * &
                (WW(I+1,J,K)-WW(I-1,J,K) + &
                 WW(I+1,J,K-1)-WW(I-1,J,K-1))
             DWDY = 0.25_EB*M%RDY(J) * &
                (WW(I,J+1,K)-WW(I,J-1,K) + &
                 WW(I,J+1,K-1)-WW(I,J-1,K-1))
             A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
             A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
             A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ

             CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)

             M%MU(I,J,K) = M%MU_DNS(I,J,K) + &
                            RHOP(I,J,K)*NU_EDDY

         ENDDO
      ENDDO
   ENDDO

ENDIF LES_IF

T_USED(14) = T_USED(14) + CURRENT_TIME() - TNOW
IF (TIME_CC_IBM) &
   T_CC_USED(CC_COMPUTE_VISCOSITY_TIME_INDEX) = &
   T_CC_USED(CC_COMPUTE_VISCOSITY_TIME_INDEX) + &
   CURRENT_TIME() - TNOW


RETURN
END SUBROUTINE CC_COMPUTE_VISCOSITY


! ----------------------- CC_STORE_FACE_FV --------------------------------

SUBROUTINE CC_STORE_FACE_FV(M,SUBSTRACT_BAROCLINIC)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
LOGICAL, INTENT(IN) :: SUBSTRACT_BAROCLINIC

INTEGER :: ICF, I, J, K, X1AXIS

IF (.NOT.SUBSTRACT_BAROCLINIC) THEN

   CF_LOOP_1 : DO ICF=1,M%N_CUTFACE_MESH
      IF (M%CUT_FACE(ICF)%STATUS /= CC_GASPHASE) &
         CYCLE CF_LOOP_1
      M%CUT_FACE(ICF)%FV = 0._EB
      IF(M%CUT_FACE(ICF)%IWC>0) THEN
         IF(M%WALL(M%CUT_FACE(ICF)%IWC)%BOUNDARY_TYPE &
            ==MIRROR_BOUNDARY) CYCLE CF_LOOP_1
      ENDIF
      I      = M%CUT_FACE(ICF)%IJK(IAXIS)
      J      = M%CUT_FACE(ICF)%IJK(JAXIS)
      K      = M%CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = M%CUT_FACE(ICF)%IJK(KAXIS+1)
      SELECT CASE(X1AXIS)
            CASE(IAXIS)
               M%CUT_FACE(ICF)%FV = M%FVX(I,J,K)
            CASE(JAXIS)
               M%CUT_FACE(ICF)%FV = M%FVY(I,J,K)
            CASE(KAXIS)
               M%CUT_FACE(ICF)%FV = M%FVZ(I,J,K)
         END SELECT
   ENDDO CF_LOOP_1

ELSE

   CF_LOOP_2 : DO ICF=1,M%N_CUTFACE_MESH
      IF (M%CUT_FACE(ICF)%STATUS /= CC_GASPHASE) &
         CYCLE CF_LOOP_2
      M%CUT_FACE(ICF)%FV = 0._EB
      IF(M%CUT_FACE(ICF)%IWC>0) THEN
         IF(M%WALL(M%CUT_FACE(ICF)%IWC)%BOUNDARY_TYPE &
            ==MIRROR_BOUNDARY) CYCLE CF_LOOP_2
      ENDIF
      I      = M%CUT_FACE(ICF)%IJK(IAXIS)
      J      = M%CUT_FACE(ICF)%IJK(JAXIS)
      K      = M%CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = M%CUT_FACE(ICF)%IJK(KAXIS+1)
      SELECT CASE(X1AXIS)
            CASE(IAXIS)
               M%CUT_FACE(ICF)%FV = &
                  M%FVX(I,J,K) - M%FVX_B(I,J,K)
            CASE(JAXIS)
               M%CUT_FACE(ICF)%FV = &
                  M%FVY(I,J,K) - M%FVY_B(I,J,K)
            CASE(KAXIS)
               M%CUT_FACE(ICF)%FV = &
                  M%FVZ(I,J,K) - M%FVZ_B(I,J,K)
         END SELECT
   ENDDO CF_LOOP_2

ENDIF

RETURN
END SUBROUTINE CC_STORE_FACE_FV


! --------------------------- CC_VELOCITY_FLUX ------------------------------

SUBROUTINE CC_VELOCITY_FLUX(M,DT, &
   APPLY_TO_ESTIMATED_VARIABLES,RHOP,CORRECT_GRAV, &
   GX,GY,GZ)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES,CORRECT_GRAV
REAL(EB),INTENT(IN) :: DT
REAL(EB),INTENT(IN), OPTIONAL :: &
   GX(0:IBAR_MAX),GY(0:IBAR_MAX),GZ(0:IBAR_MAX)
REAL(EB), INTENT(IN), POINTER, DIMENSION(:,:,:) :: RHOP

! Local Variables:
TYPE(CC_CUTCELL_TYPE), POINTER :: CC
TYPE(CC_CUTFACE_TYPE), POINTER :: CF
TYPE(CC_RCFACE_TYPE),  POINTER :: RCF
INTEGER  :: ICF,IRC,I,J,K,X1AXIS,JCF,ICC,JCC,ISIDE
REAL(EB) :: TNOW, TNOW2, RRHO, X1F, IDX, CCM1, CCP1, &
            RHOV(-1:0), RHO0V(-1:0), FVCC, PRFCT

IF ( FREEZE_VELOCITY .OR. SOLID_PHASE_ONLY ) RETURN
IF (PERIODIC_TEST == 103 .OR. PERIODIC_TEST == 11 .OR. &
    PERIODIC_TEST==7) RETURN

! Dummy for now:
TNOW = DT
TNOW = CURRENT_TIME()

CORRECT_GRAV_IF : IF(CORRECT_GRAV) THEN

PRFCT = 0._EB
IF(APPLY_TO_ESTIMATED_VARIABLES) PRFCT=1._EB
WHERE(M%FCVAR(0:M%IBAR,1:M%JBAR,1:M%KBAR,CC_FGSC,IAXIS) &
   ==CC_SOLID) &
   M%FVX(0:M%IBAR,1:M%JBAR,1:M%KBAR) = 0._EB
WHERE(M%FCVAR(1:M%IBAR,0:M%JBAR,1:M%KBAR,CC_FGSC,JAXIS) &
   ==CC_SOLID) &
   M%FVY(1:M%IBAR,0:M%JBAR,1:M%KBAR) = 0._EB
WHERE(M%FCVAR(1:M%IBAR,1:M%JBAR,0:M%KBAR,CC_FGSC,KAXIS) &
   ==CC_SOLID) &
   M%FVZ(1:M%IBAR,1:M%JBAR,0:M%KBAR) = 0._EB
! Correct the gravity terms for RC and cut-faces:
CUTFACE_LOOP_0 : DO ICF=1,M%N_CUTFACE_MESH
   CF => M%CUT_FACE(ICF)
   IF(CF%STATUS /= CC_GASPHASE) CYCLE CUTFACE_LOOP_0
   I = CF%IJK(IAXIS); J = CF%IJK(JAXIS)
   K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   ! Define cut-face gravity terms sum:
   DO JCF=1,CF%NFACE
      X1F= CF%XYZCEN(X1AXIS,JCF)
      IDX= 1._EB / &
         (CF%XCENHIGH(X1AXIS,JCF)-CF%XCENLOW(X1AXIS,JCF))
      CCM1= IDX*(CF%XCENHIGH(X1AXIS,JCF)-X1F)
      CCP1= IDX*(X1F-CF%XCENLOW(X1AXIS, JCF))
      RHOV(-1:0) = 0._EB; RHO0V(-1:0) = 0._EB
      DO ISIDE=-1,0
         SELECT CASE(CF%CELL_LIST(1,ISIDE+2,JCF))
         CASE(CC_FTYPE_CFGAS)
            ICC = CF%CELL_LIST(2,ISIDE+2,JCF)
            CC=>M%CUT_CELL(ICC)
            JCC = CF%CELL_LIST(3,ISIDE+2,JCF)
            RHOV(ISIDE) = PRFCT*CC%RHOS(JCC) + &
                           (1._EB-PRFCT)*CC%RHO(JCC)
            RHO0V(ISIDE)= CC%RHO_0(JCC)
         END SELECT
      ENDDO
      CF%FN(JCF) = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
                   (CCM1*RHOV(-1) + CCP1*RHOV(0))
   ENDDO
   FVCC = DOT_PRODUCT(CF%FN(1:CF%NFACE), &
      CF%AREA(1:CF%NFACE)) / SUM(CF%AREA(1:CF%NFACE))
   ! Apply unstructured - cartesian correction:
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I+1,J,K))
      M%FVX(I,J,K) = M%FVX(I,J,K) + &
         GX(I)*(FVCC - RRHO*M%RHO_0(K))
   CASE(JAXIS)
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J+1,K))
      M%FVY(I,J,K) = M%FVY(I,J,K) + &
         GY(I)*(FVCC - RRHO*M%RHO_0(K))
   CASE(KAXIS)
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J,K+1))
      M%FVZ(I,J,K) = M%FVZ(I,J,K) + &
         GZ(I)*(FVCC - RRHO*0.5_EB * &
         (M%RHO_0(K)+M%RHO_0(K+1)))
   END SELECT
ENDDO CUTFACE_LOOP_0

RC_FACE_LOOP : DO IRC=1,M%CC_NRCFACE_Z
   RCF    =>M%RC_FACE(IRC)
   I      = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS)
   K      = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
   IDX = 1._EB / &
      (RCF%XCEN(X1AXIS,HIGH_IND)-RCF%XCEN(X1AXIS,LOW_IND))
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      CCM1= IDX*(RCF%XCEN(X1AXIS,HIGH_IND)-M%X(I))
      CCP1= IDX*(M%X(I)-RCF%XCEN(X1AXIS,LOW_IND))
      RHOV(-1:0) = RHOP(I:I+1,J,K)
      RHO0V(-1:0) = M%RHO_0(K)
      IF(M%CCVAR(I  ,J,K,CC_UNKZ)>0) &
         RHO0V(-1) = RHO_0_CV( &
         M%CCVAR(I  ,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
      IF(M%CCVAR(I+1,J,K,CC_UNKZ)>0) &
         RHO0V( 0) = RHO_0_CV( &
         M%CCVAR(I+1,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
      DO ISIDE=-1,0
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            CC=>M%CUT_CELL(ICC)
            RHOV(ISIDE) = PRFCT *CC%RHOS(JCC) + &
               (1._EB-PRFCT) *CC%RHO(JCC)
            RHO0V(ISIDE)= CC%RHO_0(JCC)
         END SELECT
      ENDDO
      FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
             (CCM1*RHOV(-1) + CCP1*RHOV(0))
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I+1,J,K))
      M%FVX(I,J,K) = M%FVX(I,J,K) + &
         GX(I)*(FVCC - RRHO*M%RHO_0(K))
   CASE(JAXIS)
      CCM1= IDX*(RCF%XCEN(X1AXIS,HIGH_IND)-M%Y(J))
      CCP1= IDX*(M%Y(J)-RCF%XCEN(X1AXIS,LOW_IND))
      RHOV(-1:0) = RHOP(I,J:J+1,K)
      RHO0V(-1:0) = M%RHO_0(K)
      IF(M%CCVAR(I,J  ,K,CC_UNKZ)>0) &
         RHO0V(-1) = RHO_0_CV( &
         M%CCVAR(I,J  ,K,CC_UNKZ)-UNKZ_IND(NM_START))
      IF(M%CCVAR(I,J+1,K,CC_UNKZ)>0) &
         RHO0V( 0) = RHO_0_CV( &
         M%CCVAR(I,J+1,K,CC_UNKZ)-UNKZ_IND(NM_START))
      DO ISIDE=-1,0
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            CC=>M%CUT_CELL(ICC)
            RHOV(ISIDE) = PRFCT *CC%RHOS(JCC) + &
               (1._EB-PRFCT) *CC%RHO(JCC)
            RHO0V(ISIDE)= CC%RHO_0(JCC)
         END SELECT
      ENDDO
      FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
             (CCM1*RHOV(-1) + CCP1*RHOV(0))
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J+1,K))
      M%FVY(I,J,K) = M%FVY(I,J,K) + &
         GY(I)*(FVCC - RRHO*M%RHO_0(K))
   CASE(KAXIS)
      CCM1= IDX*(RCF%XCEN(X1AXIS,HIGH_IND)-M%Z(K))
      CCP1= IDX*(M%Z(K)-RCF%XCEN(X1AXIS,LOW_IND))
      RHOV(-1:0) = RHOP(I,J,K:K+1)
      RHO0V(-1:0) = M%RHO_0(K:K+1)
      IF(M%CCVAR(I,J,K  ,CC_UNKZ)>0) &
         RHO0V(-1) = RHO_0_CV( &
         M%CCVAR(I,J,K  ,CC_UNKZ)-UNKZ_IND(NM_START))
      IF(M%CCVAR(I,J,K+1,CC_UNKZ)>0) &
         RHO0V( 0) = RHO_0_CV( &
         M%CCVAR(I,J,K+1,CC_UNKZ)-UNKZ_IND(NM_START))
      DO ISIDE=-1,0
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            CC=>M%CUT_CELL(ICC)
            RHOV(ISIDE) = PRFCT *CC%RHOS(JCC) + &
               (1._EB-PRFCT) *CC%RHO(JCC)
            RHO0V(ISIDE)= CC%RHO_0(JCC)
         END SELECT
      ENDDO
      FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
             (CCM1*RHOV(-1) + CCP1*RHOV(0))
      RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J,K+1))
      M%FVZ(I,J,K) = M%FVZ(I,J,K) + &
         GZ(I)*(FVCC - RRHO*0.5_EB * &
         (M%RHO_0(K)+M%RHO_0(K+1)))
   END SELECT
ENDDO RC_FACE_LOOP

! Faces with two regular gas cells, one belonging to a linked CV:
CCM1=0.5_EB; CCP1=0.5_EB
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         IF(ANY(M%CCVAR(I:I+1,J,K,CC_CGSC) /= &
            CC_GASPHASE)) CYCLE
         IF(ALL(M%CCVAR(I:I+1,J,K,CC_UNKZ) <= 0)) CYCLE
         RHOV(-1:0) = RHOP(I:I+1,J,K)
         RHO0V(-1:0) = M%RHO_0(K)
         IF(M%CCVAR(I  ,J,K,CC_UNKZ)>0) &
            RHO0V(-1) = RHO_0_CV( &
            M%CCVAR(I  ,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
         IF(M%CCVAR(I+1,J,K,CC_UNKZ)>0) &
            RHO0V( 0) = RHO_0_CV( &
            M%CCVAR(I+1,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
         FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
                (CCM1*RHOV(-1) + CCP1*RHOV(0))
         RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I+1,J,K))
         M%FVX(I,J,K) = M%FVX(I,J,K) + &
            GX(I)*(FVCC - RRHO*M%RHO_0(K))
      ENDDO
   ENDDO
ENDDO

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         IF(ANY(M%CCVAR(I,J:J+1,K,CC_CGSC) /= &
            CC_GASPHASE)) CYCLE
         IF(ALL(M%CCVAR(I,J:J+1,K,CC_UNKZ) <= 0)) CYCLE
         RHOV(-1:0) = RHOP(I,J:J+1,K)
         RHO0V(-1:0) = M%RHO_0(K)
         IF(M%CCVAR(I,J  ,K,CC_UNKZ)>0) &
            RHO0V(-1) = RHO_0_CV( &
            M%CCVAR(I,J  ,K,CC_UNKZ)-UNKZ_IND(NM_START))
         IF(M%CCVAR(I,J+1,K,CC_UNKZ)>0) &
            RHO0V( 0) = RHO_0_CV( &
            M%CCVAR(I,J+1,K,CC_UNKZ)-UNKZ_IND(NM_START))
         FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
                (CCM1*RHOV(-1) + CCP1*RHOV(0))
         RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J+1,K))
         M%FVY(I,J,K) = M%FVY(I,J,K) + &
            GY(I)*(FVCC - RRHO*M%RHO_0(K))
      ENDDO
   ENDDO
ENDDO

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF(ANY(M%CCVAR(I,J,K:K+1,CC_CGSC) /= &
            CC_GASPHASE)) CYCLE
         IF(ALL(M%CCVAR(I,J,K:K+1,CC_UNKZ) <= 0)) CYCLE
         RHOV(-1:0) = RHOP(I,J,K:K+1)
         RHO0V(-1:0) = M%RHO_0(K:K+1)
         IF(M%CCVAR(I,J,K  ,CC_UNKZ)>0) &
            RHO0V(-1) = RHO_0_CV( &
            M%CCVAR(I,J,K  ,CC_UNKZ)-UNKZ_IND(NM_START))
         IF(M%CCVAR(I,J,K+1,CC_UNKZ)>0) &
            RHO0V( 0) = RHO_0_CV( &
            M%CCVAR(I,J,K+1,CC_UNKZ)-UNKZ_IND(NM_START))
         FVCC = (CCM1*RHO0V(-1) + CCP1*RHO0V(0)) / &
                (CCM1*RHOV(-1) + CCP1*RHOV(0))
         RRHO = 2._EB / (RHOP(I,J,K)+RHOP(I,J,K+1))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + &
            GZ(I)*(FVCC - RRHO*0.5_EB * &
            (M%RHO_0(K)+M%RHO_0(K+1)))
      ENDDO
   ENDDO
ENDDO

ELSE CORRECT_GRAV_IF

! For now add the value of FV into CUT_FACE(ICF)%FN
CALL CC_STORE_FACE_FV(M,SUBSTRACT_BAROCLINIC=.FALSE.)

CUTFACE_LOOP : DO ICF=1,M%N_CUTFACE_MESH
   IF (M%CUT_FACE(ICF)%STATUS /= CC_GASPHASE) &
      CYCLE CUTFACE_LOOP
   M%CUT_FACE(ICF)%FN(1:M%CUT_FACE(ICF)%NFACE) = &
      M%CUT_FACE(ICF)%FV
ENDDO CUTFACE_LOOP

ENDIF CORRECT_GRAV_IF

T_USED(14) = T_USED(14) + CURRENT_TIME() - TNOW
IF (TIME_CC_IBM) &
   T_CC_USED(CC_VELOCITY_FLUX_TIME_INDEX) = &
   T_CC_USED(CC_VELOCITY_FLUX_TIME_INDEX) + &
   CURRENT_TIME() - TNOW2
RETURN
END SUBROUTINE CC_VELOCITY_FLUX


! -------------------------------- CC_PROJECT_VELOCITY_KERNEL -----------------------------------
!> Thread-safe version of CC_PROJECT_VELOCITY.
!> Projects velocities onto cut-cell faces for CC_IBM geometry.
!> Uses local pointer aliases (Pattern 3) to shadow module-level MESH_POINTERS names.
!>
!> Three modes:
!>   STORE_FLG=.TRUE.: Corrector store — save U/V/W into M%U_STORE_CC etc.
!>   STORE_FLG=.FALSE., PREDICTOR_FLAG=.TRUE.: Predictor — update US/VS/WS cut-face velocities
!>   STORE_FLG=.FALSE., PREDICTOR_FLAG=.FALSE.: Corrector — update U/V/W using stored values

RECURSIVE SUBROUTINE CC_PROJECT_VELOCITY_KERNEL(M, NM, DT, STORE_FLG, PREDICTOR_FLAG)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
LOGICAL, INTENT(IN) :: STORE_FLG, PREDICTOR_FLAG

! Local aliases shadowing module-level MESH_POINTERS names
REAL(EB), POINTER, DIMENSION(:,:,:) :: U, V, W, US, VS, WS, H, HS
REAL(EB), POINTER, DIMENSION(:,:,:) :: FVX, FVY, FVZ
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ
INTEGER :: IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, N_EXTERNAL_WALL_CELLS
INTEGER, POINTER, DIMENSION(:,:,:,:,:) :: FCVAR
TYPE(CC_CUTFACE_TYPE), POINTER, DIMENSION(:) :: CUT_FACE
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_RCFACE_TYPE), POINTER, DIMENSION(:) :: RC_FACE
TYPE(WALL_TYPE), POINTER, DIMENSION(:) :: WALL
TYPE(BOUNDARY_COORD_TYPE), POINTER, DIMENSION(:) :: BOUNDARY_COORD

! Local variables
TYPE(CC_CUTFACE_TYPE), POINTER :: CF
TYPE(CC_RCFACE_TYPE), POINTER :: RCF
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER :: I,J,K,ICF,JCF,X1AXIS,IFACE,IOR,IW,IRC
REAL(EB) :: IDX,H_HI,H_LO,FCTH

! Set aliases from mesh
U => M%U; V => M%V; W => M%W
US => M%US; VS => M%VS; WS => M%WS
H => M%H; HS => M%HS
FVX => M%FVX; FVY => M%FVY; FVZ => M%FVZ
DX => M%DX; DY => M%DY; DZ => M%DZ
IBAR = M%IBAR; JBAR = M%JBAR; KBAR = M%KBAR
IBP1 = M%IBP1; JBP1 = M%JBP1; KBP1 = M%KBP1
N_EXTERNAL_WALL_CELLS = M%N_EXTERNAL_WALL_CELLS
FCVAR => M%FCVAR
CUT_FACE => M%CUT_FACE
CUT_CELL => M%CUT_CELL
RC_FACE => M%RC_FACE
WALL => M%WALL
BOUNDARY_COORD => M%BOUNDARY_COORD

STORE_IF : IF (STORE_FLG) THEN

   IF (ALLOCATED(M%U_STORE_CC)) DEALLOCATE(M%U_STORE_CC)
   IF (ALLOCATED(M%V_STORE_CC)) DEALLOCATE(M%V_STORE_CC)
   IF (ALLOCATED(M%W_STORE_CC)) DEALLOCATE(M%W_STORE_CC)
   ALLOCATE(M%U_STORE_CC(0:IBP1,0:JBP1,0:KBP1))
   ALLOCATE(M%V_STORE_CC(0:IBP1,0:JBP1,0:KBP1))
   ALLOCATE(M%W_STORE_CC(0:IBP1,0:JBP1,0:KBP1))

   M%U_STORE_CC = U
   M%V_STORE_CC = V
   M%W_STORE_CC = W

ELSE STORE_IF

   PRED_CORR_IF : IF (PREDICTOR_FLAG) THEN

      ! Update INBOUNDARY faces:
      DO ICF=1,M%N_CUTFACE_MESH
         CF => CUT_FACE(ICF); IF(CF%STATUS /= CC_INBOUNDARY) CYCLE
         CF%VELS(1:CF%NFACE) = CF%VEL(1:CF%NFACE) - DT*CF%FN(1:CF%NFACE)
      ENDDO

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,IAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(IAXIS,JCF)-CF%XCENLOW(IAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%H( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%H( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(IAXIS,JCF)-CF%XCENLOW(IAXIS,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H(I+1,J,K)-H(I,J,K)) )
                     ENDDO
                  ENDIF
                  US(I,J,K) = DOT_PRODUCT(CF%VELS(1:CF%NFACE), &
                               CF%AREA(1:CF%NFACE)) / (DY(J)*DZ(K))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,JAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(JAXIS,JCF)-CF%XCENLOW(JAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%H( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%H( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(JAXIS,JCF)-CF%XCENLOW(JAXIS,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H(I,J+1,K)-H(I,J,K)) )
                     ENDDO
                  ENDIF
                  VS(I,J,K) = DOT_PRODUCT(CF%VELS(1:CF%NFACE), &
                               CF%AREA(1:CF%NFACE)) / (DX(I)*DZ(K))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,KAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(KAXIS,JCF)-CF%XCENLOW(KAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%H( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%H( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX  = 1._EB/(CF%XCENHIGH(KAXIS,JCF)-CF%XCENLOW(KAXIS,JCF))
                        CF%VELS(JCF) = CF%VEL(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H(I,J,K+1)-H(I,J,K)) )
                     ENDDO
                  ENDIF
                  WS(I,J,K) = DOT_PRODUCT(CF%VELS(1:CF%NFACE), &
                               CF%AREA(1:CF%NFACE)) / (DY(J)*DX(I))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Regular faces connecting gasphase-gasphase or gasphase-cut-cells:
      DO IFACE=1,M%CC_NRCFACE_H
         RCF => RC_FACE(M%RCF_H(IFACE))
         FCTH = 1._EB; IF(RCF%IWC>0 .AND. &
            ANY(WALL(RCF%IWC)%BOUNDARY_TYPE== &
            (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
         I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS)
         K   = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
         IDX = 1._EB / ( RCF%XCEN(X1AXIS,HIGH_IND) - &
                          RCF%XCEN(X1AXIS,LOW_IND) )
         SELECT CASE(X1AXIS)
            CASE(IAXIS)
               US(I,J,K) = U(I,J,K) - DT*( FVX(I,J,K) + &
                            FCTH*IDX*(H(I+1,J,K)-H(I,J,K)) )
            CASE(JAXIS)
               VS(I,J,K) = V(I,J,K) - DT*( FVY(I,J,K) + &
                            FCTH*IDX*(H(I,J+1,K)-H(I,J,K)) )
            CASE(KAXIS)
               WS(I,J,K) = W(I,J,K) - DT*( FVZ(I,J,K) + &
                            FCTH*IDX*(H(I,J,K+1)-H(I,J,K)) )
         END SELECT
      ENDDO

      ! RC faces in OPEN Boundaries:
      WALL_CELL_LOOP_1 : DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF(.NOT.(WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. &
            (PRES_FLAG==ULMAT_FLAG .AND. &
             WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY)) ) &
            CYCLE WALL_CELL_LOOP_1
         BC  => BOUNDARY_COORD(WC%BC_INDEX)
         I = BC%IIG; J = BC%JJG; K = BC%KKG; IOR = BC%IOR
         SELECT CASE (IOR)
         CASE( 1); I = BC%IIG-1
         CASE( 2); J = BC%JJG-1
         CASE( 3); K = BC%KKG-1
         END SELECT
         IRC = FCVAR(I,J,K,CC_IDRC,ABS(IOR))
         IF(IRC < 1) CYCLE WALL_CELL_LOOP_1
         IDX = 1._EB/( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - &
                        RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
         SELECT CASE (ABS(IOR))
         CASE(1); US(I,J,K)= U(I,J,K) - DT*( FVX(I,J,K) + &
                              IDX*(H(I+1,J,K)-H(I,J,K)) )
         CASE(2); VS(I,J,K)= V(I,J,K) - DT*( FVY(I,J,K) + &
                              IDX*(H(I,J+1,K)-H(I,J,K)) )
         CASE(3); WS(I,J,K)= W(I,J,K) - DT*( FVZ(I,J,K) + &
                              IDX*(H(I,J,K+1)-H(I,J,K)) )
         END SELECT
      ENDDO WALL_CELL_LOOP_1

      WHERE(FCVAR(0:IBAR,1:JBAR,1:KBAR,CC_FGSC,IAXIS)==CC_SOLID) &
         US(0:IBAR,1:JBAR,1:KBAR) = 0._EB
      WHERE(FCVAR(1:IBAR,0:JBAR,1:KBAR,CC_FGSC,JAXIS)==CC_SOLID) &
         VS(1:IBAR,0:JBAR,1:KBAR) = 0._EB
      WHERE(FCVAR(1:IBAR,1:JBAR,0:KBAR,CC_FGSC,KAXIS)==CC_SOLID) &
         WS(1:IBAR,1:JBAR,0:KBAR) = 0._EB

   ELSE PRED_CORR_IF

      ! Update INBOUNDARY faces:
      DO ICF=1,M%N_CUTFACE_MESH
         CF => CUT_FACE(ICF); IF(CF%STATUS /= CC_INBOUNDARY) CYCLE
         CF%VEL(1:CF%NFACE) = 0.5_EB*( CF%VEL(1:CF%NFACE) + &
            CF%VELS(1:CF%NFACE) - DT*CF%FN(1:CF%NFACE) )
      ENDDO

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,IAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(IAXIS,JCF)-CF%XCENLOW(IAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%HS( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%HS( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO)) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(IAXIS,JCF)-CF%XCENLOW(IAXIS,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(HS(I+1,J,K)-HS(I,J,K))) )
                     ENDDO
                  ENDIF
                  U(I,J,K) = DOT_PRODUCT(CF%VEL(1:CF%NFACE), &
                              CF%AREA(1:CF%NFACE)) / (DY(J)*DZ(K))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,JAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(JAXIS,JCF)-CF%XCENLOW(JAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%HS( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%HS( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO)) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(JAXIS,JCF)-CF%XCENLOW(JAXIS,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(HS(I,J+1,K)-HS(I,J,K))) )
                     ENDDO
                  ENDIF
                  V(I,J,K) = DOT_PRODUCT(CF%VEL(1:CF%NFACE), &
                              CF%AREA(1:CF%NFACE)) / (DX(I)*DZ(K))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ICF = FCVAR(I,J,K,CC_IDCF,KAXIS)
               IF (ICF>0) THEN
                  CF => CUT_FACE(ICF); FCTH = 1._EB
                  IF(CF%IWC>0 .AND. &
                     ANY(WALL(CF%IWC)%BOUNDARY_TYPE== &
                     (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
                  IF (ONE_UNKH_PER_CUTCELL) THEN
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(KAXIS,JCF)-CF%XCENLOW(KAXIS,JCF))
                        H_HI = CUT_CELL(CF%CELL_LIST(2,HIGH_IND,JCF))%HS( &
                               CF%CELL_LIST(3,HIGH_IND,JCF))
                        H_LO = CUT_CELL(CF%CELL_LIST(2, LOW_IND,JCF))%HS( &
                               CF%CELL_LIST(3, LOW_IND,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(H_HI-H_LO)) )
                     ENDDO
                  ELSE
                     DO JCF=1,CF%NFACE
                        IDX = 1._EB/(CF%XCENHIGH(KAXIS,JCF)-CF%XCENLOW(KAXIS,JCF))
                        CF%VEL(JCF) = 0.5_EB*( CF%VEL(JCF) + CF%VELS(JCF) - &
                           DT*( CF%FN(JCF) + FCTH*IDX*(HS(I,J,K+1)-HS(I,J,K))) )
                     ENDDO
                  ENDIF
                  W(I,J,K) = DOT_PRODUCT(CF%VEL(1:CF%NFACE), &
                              CF%AREA(1:CF%NFACE)) / (DY(J)*DX(I))
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Regular faces connecting gasphase-gasphase or gasphase-cut-cells:
      DO IFACE=1,M%CC_NRCFACE_H
         RCF => RC_FACE(M%RCF_H(IFACE))
         FCTH = 1._EB; IF(RCF%IWC>0 .AND. &
            ANY(WALL(RCF%IWC)%BOUNDARY_TYPE== &
            (/SOLID_BOUNDARY,NULL_BOUNDARY,MIRROR_BOUNDARY/))) FCTH=0._EB
         I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS)
         K   = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
         IDX = 1._EB / ( RCF%XCEN(X1AXIS,HIGH_IND) - &
                          RCF%XCEN(X1AXIS,LOW_IND) )
         SELECT CASE(X1AXIS)
            CASE(IAXIS)
               U(I,J,K) = 0.5_EB*( M%U_STORE_CC(I,J,K) + US(I,J,K) - &
                  DT*(FVX(I,J,K) + FCTH*IDX*(HS(I+1,J,K)-HS(I,J,K))) )
            CASE(JAXIS)
               V(I,J,K) = 0.5_EB*( M%V_STORE_CC(I,J,K) + VS(I,J,K) - &
                  DT*(FVY(I,J,K) + FCTH*IDX*(HS(I,J+1,K)-HS(I,J,K))) )
            CASE(KAXIS)
               W(I,J,K) = 0.5_EB*( M%W_STORE_CC(I,J,K) + WS(I,J,K) - &
                  DT*(FVZ(I,J,K) + FCTH*IDX*(HS(I,J,K+1)-HS(I,J,K))) )
         END SELECT
      ENDDO

      WALL_CELL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS
         WC => WALL(IW)
         IF(.NOT.(WC%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. &
            (PRES_FLAG==ULMAT_FLAG .AND. &
             WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY)) ) &
            CYCLE WALL_CELL_LOOP_2
         BC  => BOUNDARY_COORD(WC%BC_INDEX)
         I = BC%IIG; J = BC%JJG; K = BC%KKG; IOR = BC%IOR
         SELECT CASE (IOR)
         CASE( 1); I = BC%IIG-1
         CASE( 2); J = BC%JJG-1
         CASE( 3); K = BC%KKG-1
         END SELECT
         IRC = FCVAR(I,J,K,CC_IDRC,ABS(IOR))
         IF(IRC < 1) CYCLE WALL_CELL_LOOP_2
         IDX = 1._EB/( RC_FACE(IRC)%XCEN(ABS(BC%IOR),HIGH_IND) - &
                        RC_FACE(IRC)%XCEN(ABS(BC%IOR),LOW_IND) )
         SELECT CASE (ABS(IOR))
         CASE(1); U(I,J,K) = 0.5_EB*( M%U_STORE_CC(I,J,K) + US(I,J,K) - &
                     DT*(FVX(I,J,K) + IDX*(HS(I+1,J,K)-HS(I,J,K))) )
         CASE(2); V(I,J,K) = 0.5_EB*( M%V_STORE_CC(I,J,K) + VS(I,J,K) - &
                     DT*(FVY(I,J,K) + IDX*(HS(I,J+1,K)-HS(I,J,K))) )
         CASE(3); W(I,J,K) = 0.5_EB*( M%W_STORE_CC(I,J,K) + WS(I,J,K) - &
                     DT*(FVZ(I,J,K) + IDX*(HS(I,J,K+1)-HS(I,J,K))) )
         END SELECT
      ENDDO WALL_CELL_LOOP_2

      DEALLOCATE(M%U_STORE_CC, M%V_STORE_CC, M%W_STORE_CC)

      WHERE(FCVAR(0:IBAR,1:JBAR,1:KBAR,CC_FGSC,IAXIS)==CC_SOLID) &
         U(0:IBAR,1:JBAR,1:KBAR) = 0._EB
      WHERE(FCVAR(1:IBAR,0:JBAR,1:KBAR,CC_FGSC,JAXIS)==CC_SOLID) &
         V(1:IBAR,0:JBAR,1:KBAR) = 0._EB
      WHERE(FCVAR(1:IBAR,1:JBAR,0:KBAR,CC_FGSC,KAXIS)==CC_SOLID) &
         W(1:IBAR,1:JBAR,0:KBAR) = 0._EB

   ENDIF PRED_CORR_IF

ENDIF STORE_IF

END SUBROUTINE CC_PROJECT_VELOCITY_KERNEL


END MODULE CC_VELOCITY_KERNELS
