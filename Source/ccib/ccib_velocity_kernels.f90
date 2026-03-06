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
                            CC_IDCF, CC_UNKZ, NM_START, &
                            CC_VELOCITY_FLUX_TIME_INDEX, &
                            CC_COMPUTE_VISCOSITY_TIME_INDEX, &
                            T_CC_USED
USE CC_SCALARS_DATA, ONLY: TIME_CC_IBM, &
                           UNKZ_IND, RHO_0_CV
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CUTFACE_VELOCITIES, CC_CUTCELL_VELOCITY, CC_COMPUTE_KRES, &
          CC_COMPUTE_VISCOSITY, CC_STORE_FACE_FV, CC_VELOCITY_FLUX

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

END MODULE CC_VELOCITY_KERNELS
