!  +++++++++++++++++++++++ CC_DENSITY_MOD ++++++++++++++++++++++++++

! Density and species transport routines for the
! cut-cell / immersed-boundary method.

MODULE CC_DENSITY_MOD

USE CC_SCALARS_DATA
USE CC_VERIFICATION, ONLY: CC_ROTATED_CUBE_RHS_ZZ
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CC_DENSITY, CC_DENSITY_TS

CONTAINS



! --------------------------------- GET_SHUNN3_QZ --------------------------------

SUBROUTINE GET_SHUNN3_QZ(T,N)

USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_Z_SRC

REAL(EB),INTENT(IN) :: T
INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER I,J,K,NM,IROW,ICC,JCC
REAL(EB) :: FCT,XHAT,ZHAT,Q_Z

FCT=REAL(2*(1-N)+1,EB)

! Mesh Loop:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! First add Q_Z on regular cells to source F_Z:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
            ! divergence from EOS
            XHAT = XC(I) - UF_MMS*T
            ZHAT = ZC(K) - WF_MMS*T
            Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
            F_Z(IROW) = F_Z(IROW) + FCT*Q_Z*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   ! Then add Cut-cell contributions to F_Z:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      DO JCC=1,CUT_CELL(ICC)%NCELL
         IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         ! divergence from EOS
         XHAT = CUT_CELL(ICC)%XYZCEN(IAXIS,JCC) - UF_MMS*T
         ZHAT = CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) - WF_MMS*T
         Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
         F_Z(IROW) = F_Z(IROW) + FCT*Q_Z*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
END SUBROUTINE GET_SHUNN3_QZ


! ------------------------------ GET_SHUNN3_QZ_TS --------------------------
! Thread-safe per-mesh version of GET_SHUNN3_QZ.
! Processes single mesh NM without POINT_TO_MESH.

RECURSIVE SUBROUTINE GET_SHUNN3_QZ_TS(NM,M,T,N)

USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_Z_SRC

INTEGER, INTENT(IN) :: NM, N
REAL(EB),INTENT(IN) :: T
TYPE(MESH_TYPE), INTENT(IN), TARGET :: M

! Local Variables:
INTEGER I,J,K,IROW,ICC,JCC
REAL(EB) :: FCT,XHAT,ZHAT,Q_Z
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: XC, ZC, DX, DY, DZ
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
INTEGER :: IBAR, JBAR, KBAR

FCT=REAL(2*(1-N)+1,EB)

IBAR = M%IBAR; JBAR = M%JBAR; KBAR = M%KBAR
CCVAR => M%CCVAR
XC => M%XC; ZC => M%ZC
DX => M%DX; DY => M%DY; DZ => M%DZ
CUT_CELL => M%CUT_CELL

! Add Q_Z on regular cells to source F_Z:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF(CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         XHAT = XC(I) - UF_MMS*T
         ZHAT = ZC(K) - WF_MMS*T
         Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
         F_Z(IROW) = F_Z(IROW) + FCT*Q_Z*DX(I)*DY(J)*DZ(K)
      ENDDO
   ENDDO
ENDDO

! Then add Cut-cell contributions to F_Z:
DO ICC=1,M%N_CUTCELL_MESH
   DO JCC=1,CUT_CELL(ICC)%NCELL
      IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      XHAT = CUT_CELL(ICC)%XYZCEN(IAXIS,JCC) - UF_MMS*T
      ZHAT = CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) - WF_MMS*T
      Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
      F_Z(IROW) = F_Z(IROW) + FCT*Q_Z*CUT_CELL(ICC)%VOLUME(JCC)
   ENDDO
ENDDO

RETURN
END SUBROUTINE GET_SHUNN3_QZ_TS


! ------------------------------ CC_DENSITY -------------------------------

SUBROUTINE CC_DENSITY(T,DT)

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE MATH_FUNCTIONS, ONLY : EVALUATE_RAMP
USE MPI_F08

REAL(EB), INTENT(IN) :: T,DT

! Local Variables:
INTEGER :: N
INTEGER :: I,J,K,NM,ICC,JCC
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),VCCELL,PBAR_K
! CHARACTER(len=20) :: filename
! LOGICAL, SAVE :: FIRST_CALL = .TRUE.
!

IF (SOLID_PHASE_ONLY) RETURN

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN ! In order to avoid instabilities due to unphysical initial flow fields.
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
      ! CONTINUE
END SELECT

! Advance scalars and density, sanitize results if needed:
CALL CC_DENSITY_EXPLICIT(T,DT)

! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i). Here WBAR=1/SUM(Y_i/W_i).
! Compute temperature in regular and cut-cells, from equation of state:
IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN

   MESHES_LOOP1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)

      ! First Regular Cells:
      ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i).
      ! Extract predicted temperature at next time step from Equation of State
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
               PBAR_K = PBAR_S(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
               TMP(I,J,K) = PBAR_K/(RSUM(I,J,K)*RHOS(I,J,K))
            ENDDO
         ENDDO
      ENDDO

      ! Store RHO*ZZ values at step n:
      IF (.NOT.ALLOCATED(MESHES(NM)%RHO_ZZN)) ALLOCATE(MESHES(NM)%RHO_ZZN(0:IBP1,0:JBP1,0:KBP1,N_TOTAL_SCALARS))
      DO N=1,N_TOTAL_SCALARS
         MESHES(NM)%RHO_ZZN(0:IBP1,0:JBP1,0:KBP1,N) = RHO(0:IBP1,0:JBP1,0:KBP1)*ZZ(0:IBP1,0:JBP1,0:KBP1,N)
      ENDDO

      ! Second cut-cells, these variables being filled are only used for exporting to slices and applying Boundary
      ! conditions on external walls other than NULL or INTERPOLATED in WALL_BC (wall.f90):
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         I  = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
         VCCELL = 0._EB; TMP(I,J,K)=0._EB; RHOS(I,J,K)=0._EB; ZZS(I,J,K,1:N_TRACKED_SPECIES)=0._EB; RSUM(I,J,K)=0._EB
         DO JCC=1,CC%NCELL
            ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i).
            ZZ_GET(1:N_TRACKED_SPECIES) = CC%ZZS(1:N_TRACKED_SPECIES,JCC)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,CC%RSUM(JCC))
            ! Extract predicted temperature at next time step from Equation of State
            ! Use for pressure the height of the underlying cartesian cell centroid:
            PBAR_K = PBAR_S(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CC%UNKZ(JCC)-UNKZ_IND(NM_START))
            CC%TMP(JCC) = PBAR_K/(CC%RSUM(JCC)*CC%RHOS(JCC))
            TMP(I,J,K) = TMP(I,J,K) + CC%TMP(JCC)*CC%VOLUME(JCC)
            RHOS(I,J,K)= RHOS(I,J,K)+ CC%RHOS(JCC)*CC%VOLUME(JCC)
            ZZS(I,J,K,1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*CC%VOLUME(JCC)
            RSUM(I,J,K)= RSUM(I,J,K)+ CC%RSUM(JCC)*CC%VOLUME(JCC)

            VCCELL = VCCELL + CC%VOLUME(JCC)
         ENDDO

         ! Volume average cell variables to underlying cell:
         TMP(I,J,K) = TMP(I,J,K)/VCCELL
         RHOS(I,J,K)= RHOS(I,J,K)/VCCELL
         ZZS(I,J,K,1:N_TRACKED_SPECIES)=ZZS(I,J,K,1:N_TRACKED_SPECIES)/VCCELL
         RSUM(I,J,K)=RSUM(I,J,K)/VCCELL

      ENDDO

      ! Finally set to ambient temperature the temp of SOLID cells:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_CGSC) /= CC_SOLID) CYCLE
               TMP(I,J,K) = TMPA
            ENDDO
         ENDDO
      ENDDO

   ENDDO MESHES_LOOP1

ELSE ! CORRECTOR

   MESHES_LOOP2 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)

      ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)
      ! Extract predicted temperature at next time step from Equation of State
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
               PBAR_K = PBAR(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
               TMP(I,J,K) = PBAR_K/(RSUM(I,J,K)*RHO(I,J,K))
            ENDDO
         ENDDO
      ENDDO

      ! Second cut-cells, these variables being filled are only used for exporting to slices and applying Boundary
      ! conditions on external walls other than NULL or INTERPOLATED in WALL_BC (wall.f90):
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         I  = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
         VCCELL = 0._EB; TMP(I,J,K)=0._EB; RHO(I,J,K)=0._EB; ZZ(I,J,K,1:N_TRACKED_SPECIES)=0._EB; RSUM(I,J,K)=0._EB
         DO JCC=1,CC%NCELL
            ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i).
            ZZ_GET(1:N_TRACKED_SPECIES) = CC%ZZ(1:N_TRACKED_SPECIES,JCC)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,CC%RSUM(JCC))
            ! Extract predicted temperature at next time step from Equation of State
            ! Use for pressure the height of the underlying cartesian cell centroid:
            PBAR_K = PBAR(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CC%UNKZ(JCC)-UNKZ_IND(NM_START))
            CC%TMP(JCC) = PBAR_K/(CC%RSUM(JCC)*CC%RHO(JCC))
            TMP(I,J,K) = TMP(I,J,K) + CC%TMP(JCC)*CC%VOLUME(JCC)
            RHO(I,J,K) = RHO(I,J,K) + CC%RHO(JCC)*CC%VOLUME(JCC)
            ZZ(I,J,K,1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*CC%VOLUME(JCC)
            RSUM(I,J,K)= RSUM(I,J,K)+ CC%RSUM(JCC)*CC%VOLUME(JCC)

            VCCELL = VCCELL + CC%VOLUME(JCC)

         ENDDO

         ! Volume average cell variables to underlying cell:
         TMP(I,J,K) = TMP(I,J,K)/VCCELL
         RHO(I,J,K) = RHO(I,J,K)/VCCELL
         ZZ(I,J,K,1:N_TRACKED_SPECIES)=ZZ(I,J,K,1:N_TRACKED_SPECIES)/VCCELL
         RSUM(I,J,K)=RSUM(I,J,K)/VCCELL

      ENDDO

      ! Finally set to ambient temperature the temp of SOLID cells:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_CGSC) /= CC_SOLID) CYCLE
               TMP(I,J,K) = TMPA
            ENDDO
         ENDDO
      ENDDO

   ENDDO MESHES_LOOP2

ENDIF ! PREDICTOR

RETURN
END SUBROUTINE CC_DENSITY


! ----------------------------- CC_DENSITY_EXPLICIT ------------------------

SUBROUTINE CC_DENSITY_EXPLICIT(T,DT)

REAL(EB), INTENT(IN) :: T,DT

! Local variables:
INTEGER :: N
INTEGER :: IROW_LOC
REAL(EB):: DUMMYT

! Just to avoid compilation warnings: T might be used to define a time dependent source.
DUMMYT = T

! Loop through species:
! This loop performs an either implicit or explicit time advancement of the transport equations for each
! chemical species on the cut-cell implicit region, plus explicit reaction (as done on FDS).
! Scalar bounds are checked on the implicit region regular and cut-cells:
SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

   IF( (MESHES(LOWER_MESH_INDEX)%PREDICTOR.AND.FIRST_PASS) .OR. MESHES(LOWER_MESH_INDEX)%CORRECTOR) THEN
      ! RHS vector (Adv+diff)*zz+F_BC, derived from boundary conditions on immersed and domain Boundaries:
      F_Z(:) = 0._EB
      CALL GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D(N)

      ! Add Advective fluxes for F_Z:
      CALL GET_ADVDIFFVECTOR_SCALAR_3D(N)

      ! Here add the reaction source term M_DOT_PPP, treated explicitly:
      CALL GET_M_DOT_PPP_SCALAR_3D(N)

      IF (PERIODIC_TEST==7) CALL GET_SHUNN3_QZ(T,N)
      IF (PERIODIC_TEST==21 .OR. PERIODIC_TEST==22 .OR. PERIODIC_TEST==23) CALL CC_ROTATED_CUBE_RHS_ZZ(T,N)

      ! Get rho*zz vector at step n:
      CALL GET_RHOZZVECTOR_SCALAR_3D(N)
   ENDIF

   IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN
      IF (FIRST_PASS) THEN
         F_Z0(:,N) = F_Z(:)
         RZ_Z0(:,N) = RZ_Z(:)
      ELSE
         F_Z(:) = F_Z0(:,N)
         RZ_Z(:)= RZ_Z0(:,N)
      ENDIF
   ENDIF

   IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN

      ! Here F_Z: (Adv+Diff)*(rho z)^n + F^n
      ! Advance with Explicit Euler: RZ_Z = RZ_Z - DT*M_MAT_Z^-1*F_Z: where initially
      ! RZ_Z = (rho z)^n, filled in GET_RHOZZVECTOR_SCALAR_3D
      DO IROW_LOC=1,NUNKZ_LOCAL
         RZ_Z(IROW_LOC) = RZ_Z(IROW_LOC) - DT * F_Z(IROW_LOC) / M_MAT_Z(IROW_LOC)
      ENDDO

   ELSE ! CORRECTOR

      ! Here F_Z: (Adv+Diff)*(rho z)^* + F^*
      ! Advance with Corrector SSPRK2: RZ_Z = RZ_Z - DT/2*M_MAT_Z^-1*F_Z: where initially
      ! RZ_Z = 1/2*((rho z)^n + (rho z)^*)
      DO IROW_LOC=1,NUNKZ_LOCAL
         RZ_Z(IROW_LOC) = RZ_Z(IROW_LOC) - 0.5_EB * DT * F_Z(IROW_LOC) / M_MAT_Z(IROW_LOC)
      ENDDO

   ENDIF

   ! Copy back to RHOZZP and CUT_CELL:
   CALL PUT_RHOZZVECTOR_SCALAR_3D(N)

ENDDO SPECIES_LOOP

! Redistribute densities if needed.
CALL CC_CHECK_MASS_DENSITY

! Recompute RHOP, and check for positivity, define mass fraction ZZ and clip if necessary:
CALL GET_RHOZZ_CC_3D

RETURN
END SUBROUTINE CC_DENSITY_EXPLICIT


! ---------------------------- GET_M_DOT_PPP_SCALAR_3D ---------------------------

SUBROUTINE GET_M_DOT_PPP_SCALAR_3D(N)

INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER :: NM,I,J,K,IROW,ICC,JCC,NCELL

! Mesh Loop:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   ! D_SOURCE(:,:,:) and M_DOT_PPP(:,:,:,:) are allocated together,
   ! => if one is not allocated the other also is not allocated.
   IF(.NOT.ALLOCATED(MESHES(NM)%D_SOURCE)) CYCLE MESH_LOOP

   CALL POINT_TO_MESH(NM)

   ! First add M_DOT_PPP on regular cells to source F_Z:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE ! Cycle Reg cells not implicit, cut-cells
                                                  ! underlying Cartesian cells and
                                                  ! solid cells.
            IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
            F_Z(IROW) = F_Z(IROW) - M_DOT_PPP(I,J,K,N)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   ! Then add Cut-cell contributions to F_Z:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      IF (CELL(CELL_INDEX(CUT_CELL(ICC)%IJK(IAXIS),CUT_CELL(ICC)%IJK(JAXIS),CUT_CELL(ICC)%IJK(KAXIS)))%SOLID) CYCLE
      NCELL=CUT_CELL(ICC)%NCELL
      DO JCC=1,NCELL
         IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         F_Z(IROW) = F_Z(IROW) - CUT_CELL(ICC)%M_DOT_PPP(N,JCC)*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO

   ! Finally if Corrector zero out M_DOT_PPP and D_SOURCE:
   IF (MESHES(NM)%CORRECTOR) THEN
      M_DOT_PPP(:,:,:,N) = 0._EB
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         NCELL=CUT_CELL(ICC)%NCELL
         CUT_CELL(ICC)%M_DOT_PPP(N,1:NCELL) = 0._EB
      ENDDO
      IF (N == N_TOTAL_SCALARS) THEN
         D_SOURCE(:,:,:)  = 0._EB
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            NCELL=CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%D_SOURCE(1:NCELL) = 0._EB
         ENDDO
      ENDIF
   ENDIF

ENDDO MESH_LOOP

RETURN
END SUBROUTINE GET_M_DOT_PPP_SCALAR_3D


! ---------------------- GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D --------------------

SUBROUTINE GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D(N)

INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER :: NM,I,J,K
REAL(EB):: PRFCT
INTEGER :: X1AXIS,IFACE,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),ICF
INTEGER :: LOCROW_1,LOCROW_2,ILOC,IROW,ICC,JCC,ISIDE,IW
REAL(EB):: AF,VELC,RHO_Z_PV(-2:1),RHOPV(-2:1),FCT,ZZ_GET_N,FN_ZZ
REAL(EB), POINTER, DIMENSION(:,:,:)  :: RHOP,UP,VP,WP
REAL(EB), POINTER, DIMENSION(:,:,:)  :: UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:):: ZZP
TYPE(CC_REGFACEZ_TYPE),  POINTER, DIMENSION(:) :: REGFACE_Z
LOGICAL :: DO_LO,DO_HI
INTEGER :: IIG,JJG,KKG,IOR
REAL(EB) :: UN

! Initialize workspace pointers (module-level TARGET arrays from CC_SCALARS_DATA)
U_TEMP => U_WORK
F_TEMP => F_WORK
Z_TEMP => Z_WORK

! Mesh Loop:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   UU=>WORK_U
   VV=>WORK_V
   WW=>WORK_W

   IF (MESHES(NM)%PREDICTOR) THEN
      ZZP  => ZZ
      RHOP => RHO
      UU   = U
      VV   = V
      WW   = W
      PRFCT= 1._EB
      WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC=>WALL(IW)
         IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         BC=>BOUNDARY_COORD(WC%BC_INDEX)
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         IOR = BC%IOR
         SELECT CASE(WC%BOUNDARY_TYPE)
            CASE DEFAULT; CYCLE WALL_LOOP
            ! SOLID_BOUNDARY is not currently functional here, but keep for testing
            CASE(SOLID_BOUNDARY);        UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
            CASE(INTERPOLATED_BOUNDARY); UN = UVW_SAVE(IW)
         END SELECT
         SELECT CASE(IOR)
            CASE( 1); UU(IIG-1,JJG,KKG) = UN
            CASE(-1); UU(IIG,JJG,KKG)   = UN
            CASE( 2); VV(IIG,JJG-1,KKG) = UN
            CASE(-2); VV(IIG,JJG,KKG)   = UN
            CASE( 3); WW(IIG,JJG,KKG-1) = UN
            CASE(-3); WW(IIG,JJG,KKG)   = UN
         END SELECT
      ENDDO WALL_LOOP

   ELSE
      ZZP  => ZZS
      RHOP => RHOS
      UU   = US
      VV   = VS
      WW   = WS
      PRFCT= 0._EB
      WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         WC=>WALL(IW)
         IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2
         B1 => BOUNDARY_PROP1(WC%B1_INDEX)
         BC=>BOUNDARY_COORD(WC%BC_INDEX)
         IIG = BC%IIG
         JJG = BC%JJG
         KKG = BC%KKG
         IOR = BC%IOR
         SELECT CASE(WC%BOUNDARY_TYPE)
            CASE DEFAULT; CYCLE WALL_LOOP_2
            ! SOLID_BOUNDARY is not currently functional here, but keep for testing
            CASE(SOLID_BOUNDARY);        UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            CASE(INTERPOLATED_BOUNDARY); UN = UVW_SAVE(IW)
         END SELECT
         SELECT CASE(IOR)
            CASE( 1); UU(IIG-1,JJG,KKG) = UN
            CASE(-1); UU(IIG,JJG,KKG)   = UN
            CASE( 2); VV(IIG,JJG-1,KKG) = UN
            CASE(-2); VV(IIG,JJG,KKG)   = UN
            CASE( 3); WW(IIG,JJG,KKG-1) = UN
            CASE(-3); WW(IIG,JJG,KKG)   = UN
         END SELECT
      ENDDO WALL_LOOP_2
   ENDIF

   ! The use of UU, VV, WW is to maintain the divergence consistent in cells next to INTERPOLATED_BOUNDARY faces, when
   ! The solver being used is the default POISSON solver (i.e. use normal velocities with velocity error).
   UP => UU
   VP => VV
   WP => WW

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.

   ! First add advective fluxes to internal and INTERPOLATED_BOUNDARY regular and cut-cells in the CC region:
   ! IAXIS faces:
   X1AXIS = IAXIS
   REGFACE_Z => CC_REGFACE_IAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC
      IF (IW>0) THEN
         IF (.NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      ENDIF
      I     = REGFACE_Z(IFACE)%IJK(IAXIS)
      J     = REGFACE_Z(IFACE)%IJK(JAXIS)
      K     = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I  ,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I+1,J,K,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DY(J)*DZ(K)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! JAXIS faces:
   X1AXIS = JAXIS
   REGFACE_Z => CC_REGFACE_JAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      I     = REGFACE_Z(IFACE)%IJK(IAXIS)
      J     = REGFACE_Z(IFACE)%IJK(JAXIS)
      K     = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I,J  ,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I,J+1,K,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DX(I)*DZ(K)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! KAXIS faces:
   X1AXIS = KAXIS
   REGFACE_Z => CC_REGFACE_KAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      I     = REGFACE_Z(IFACE)%IJK(IAXIS)
      J     = REGFACE_Z(IFACE)%IJK(JAXIS)
      K     = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I,J,K  ,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I,J,K+1,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DX(I)*DY(J)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z

      IW = MESHES(NM)%RC_FACE(IFACE)%IWC; IF(IW > 0) CYCLE

      I      = MESHES(NM)%RC_FACE(IFACE)%IJK(IAXIS)
      J      = MESHES(NM)%RC_FACE(IFACE)%IJK(JAXIS)
      K      = MESHES(NM)%RC_FACE(IFACE)%IJK(KAXIS)
      X1AXIS = MESHES(NM)%RC_FACE(IFACE)%IJK(KAXIS+1)

      ! Unknowns on related cells:
      IND(LOW_IND)  = MESHES(NM)%RC_FACE(IFACE)%UNKZ(LOW_IND)
      IND(HIGH_IND) = MESHES(NM)%RC_FACE(IFACE)%UNKZ(HIGH_IND)

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            AF = DY(J)*DZ(K)
            RHOPV(-2:1)      = RHOP(I-1:I+2,J,K)
            ! First two cells surrounding face:
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = UU(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)
         CASE(JAXIS)
            AF = DX(I)*DZ(K)
            RHOPV(-2:1)      = RHOP(I,J-1:J+2,K)
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = VV(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)
         CASE(KAXIS)
            AF = DX(I)*DY(J)
            RHOPV(-2:1)      = RHOP(I,J,K-1:K+2)
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = WW(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)
      END SELECT

      DO ILOC=LOCROW_1,LOCROW_2
         IROW=IND_LOC(ILOC)   ! Process Local Unknown number.
         FCT = REAL(3-2*ILOC,EB)
         F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
      ENDDO

   ENDDO

   ! Now Gasphase CUT_FACES:
   DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH

      IF ( MESHES(NM)%CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE

      IW = MESHES(NM)%CUT_FACE(ICF)%IWC; IF(IW > 0) CYCLE

      I = MESHES(NM)%CUT_FACE(ICF)%IJK(IAXIS)
      J = MESHES(NM)%CUT_FACE(ICF)%IJK(JAXIS)
      K = MESHES(NM)%CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = MESHES(NM)%CUT_FACE(ICF)%IJK(KAXIS+1)

      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND
      DO IFACE=1,MESHES(NM)%CUT_FACE(ICF)%NFACE

         ! Unknowns on related cells:
         IND(LOW_IND)  = MESHES(NM)%CUT_FACE(ICF)%UNKZ(LOW_IND,IFACE)
         IND(HIGH_IND) = MESHES(NM)%CUT_FACE(ICF)%UNKZ(HIGH_IND,IFACE)

         IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

         AF = MESHES(NM)%CUT_FACE(ICF)%AREA(IFACE)

         ! Matrix coefficients for advection:
         VELC =        PRFCT *MESHES(NM)%CUT_FACE(ICF)%VEL(IFACE) + &
                (1._EB-PRFCT)*MESHES(NM)%CUT_FACE(ICF)%VELS(IFACE)

         RHOPV(-1:0)    = -1._EB
         RHO_Z_PV(-1:0) =  0._EB
         DO ISIDE=-1,0
            SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
               JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDDO
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
            ENDIF
         CASE(JAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
            ENDIF
         CASE(KAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
            ENDIF
         END SELECT
         VELC  = PRFCT *CUT_FACE(ICF)%VEL(IFACE) + (1._EB-PRFCT)*CUT_FACE(ICF)%VELS(IFACE)
         ! bar{rho*zz}:
         Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_ZZ = F_TEMP(1,1,1)

         DO ILOC=LOCROW_1,LOCROW_2
            IROW=IND_LOC(ILOC)     ! Process Local Unknown number.
            FCT = REAL(3-2*ILOC,EB)
            F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
         ENDDO

      ENDDO

   ENDDO

   ! Then add (Del rho D Del Z)*dv computed on CCDIVERGENCE_PART_1:
   ! Loop over regular cells on CC region:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE ! Cycle Reg cells not implicit, cut-cells
                                                  ! underlying Cartesian cells and
                                                  ! solid cells.
            IROW  = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
            F_Z(IROW) = F_Z(IROW) - DEL_RHO_D_DEL_Z(I,J,K,N)*(DX(I)*DY(J)*DZ(K))
         ENDDO
      ENDDO
   ENDDO

   ! Now cut-cells:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS)
      J = CUT_CELL(ICC)%IJK(JAXIS)
      K = CUT_CELL(ICC)%IJK(KAXIS)
      ! Don't count cut-cells inside an OBST:
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO JCC=1,CUT_CELL(ICC)%NCELL
         IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         F_Z(IROW) = F_Z(IROW) - CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(N,JCC)
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
END SUBROUTINE GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D


! --------------------------- CC_CHECK_MASS_DENSITY ------------------------

SUBROUTINE CC_CHECK_MASS_DENSITY

INTEGER, PARAMETER :: MAX_SURR_CELLS=20
INTEGER :: NM, NCELL, ICC, JCC, I, J, K, IFC, IFACE, IFC1, JFC1, ICC1, JCC1, IC, LOHI, ILH, X1AXIS, &
           II, JJ, KK, IIF, JJF, KKF, IRC, N
REAL(EB), POINTER, DIMENSION(:,:,:)   :: DELTA_RHO,DELTA_RHO_ZZ,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ
REAL(EB) :: MASS_C, MASS_N(1:MAX_SURR_CELLS), RHO_CELL(1:MAX_SURR_CELLS), VOL(1:MAX_SURR_CELLS), &
            RHO_CUT, SIGN_FACTOR, SUM_MASS_N, SUM_RHO_ZZ, CONST, PRFCT, CC_RHOP, CC_RHO_ZZP, RHO_ZZ_MIN, &
            RHO_ZZ_MAX, RHO_ZZ_CUT, RHO_ZZ_TEST
LOGICAL :: CLIP_RHOMIN_CC, CLIP_RHOMAX_CC, CLIP_RHO_ZZ, CLIP_RHO_ZZ_SAVE


! Loop meshes:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   DELTA_RHO => WORK4
   DELTA_RHO =  0._EB
   CLIP_RHOMIN_CC = .FALSE.
   CLIP_RHOMAX_CC = .FALSE.
   IF (MESHES(NM)%PREDICTOR) THEN
      RHOP   => RHOS
      RHO_ZZ => ZZS ! At this stage of the time step, ZZS is actually RHOS*ZZS in cells with UNKZ>0
      PRFCT=  1._EB
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH ! First compute RHOS in cut-cell region:
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID)  CYCLE
         CC%DELTA_RHO(1:CC%NCELL) = 0._EB
         DO JCC=1,CC%NCELL; CC%RHOS(JCC) = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC)); ENDDO
      ENDDO

   ELSE
      RHOP   => RHO
      RHO_ZZ => ZZ ! At this stage of the time step, ZZ is actually RHO*ZZ in cells with UNKZ>0
      PRFCT=  0._EB
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH ! First compute RHO in cut-cell region:
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID)  CYCLE
         CC%DELTA_RHO(1:CC%NCELL) = 0._EB
         DO JCC=1,CC%NCELL; CC%RHO(JCC) = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC)); ENDDO
      ENDDO

   ENDIF

   ! First compute RHOP in cut-cell region regular cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
            IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) THEN
               ! Bring back ZZ*RHO in regular gas cells from mass fractions computed done in DENSITY.
               RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)*RHOP(I,J,K)
            ELSE
               RHOP(I,J,K) = SUM(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES))
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   ! Correct density:
   ! Distribute delta_rho to neighbors, including linked cell:
   ! 1. Compute DELTA_RHO in cut-cells and regular cells.
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IC=CELL_INDEX(I,J,K); IF (CELL(IC)%SOLID) CYCLE
      JCC1_LOOP : DO JCC=1,CC%NCELL
         CC_RHOP = PRFCT*CC%RHOS(JCC) + (1._EB-PRFCT)*CC%RHO(JCC)
         IF (CC_RHOP>=RHOMIN .AND. CC_RHOP<=RHOMAX) CYCLE JCC1_LOOP
         IF (CC_RHOP<RHOMIN) THEN
            RHO_CUT = RHOMIN
            SIGN_FACTOR = 1._EB
            CLIP_RHOMIN_CC = .TRUE.
         ELSE
            RHO_CUT = RHOMAX
            SIGN_FACTOR = -1._EB
            CLIP_RHOMAX_CC = .TRUE.
         ENDIF
         MASS_C = ABS(RHO_CUT-CC_RHOP) * CC%VOLUME(JCC)

         ! Now find connected regular and cut-cells and add their contributions to Delta mass:
         MASS_N = 0._EB; NCELL=0
         DO IFC=2,CC%CCELEM(1,JCC)+1
            IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
            LOHI  = CC%FACE_LIST(2,IFACE)
            ILH   = 2*CC%FACE_LIST(2,IFACE) - 3 ! -1 for LOHI=LOW_IND, 1 for LOHI=HIGH_IND
            X1AXIS= CC%FACE_LIST(3,IFACE)
            IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
            NCELL=NCELL+1
            SELECT CASE(CC%FACE_LIST(1,IFACE))
            CASE(CC_FTYPE_CFGAS)
               IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
               ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
               RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
               VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
            CASE(CC_FTYPE_RCGAS)
               II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
               CASE(CC_FTYPE_RGGAS) ! Regular cell
                  RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell
                  ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                  RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
                  VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
               END SELECT
            END SELECT
            MASS_N(NCELL) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHO_CELL(NCELL)))-RHO_CUT) * VOL(NCELL)
         ENDDO

         ! We have MASS_C, MASS_N and NCELL, now distribute:
         SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)

         CC%DELTA_RHO(JCC) = CC%DELTA_RHO(JCC) + CONST*SUM_MASS_N/CC%VOLUME(JCC)
         NCELL=0
         DO IFC=2,CC%CCELEM(1,JCC)+1
            IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
            LOHI  = CC%FACE_LIST(2,IFACE)
            ILH   = 2*CC%FACE_LIST(2,IFACE) - 3 ! -1 for LOHI=LOW_IND, 1 for LOHI=HIGH_IND
            X1AXIS= CC%FACE_LIST(3,IFACE)
            IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
            NCELL=NCELL+1
            SELECT CASE(CC%FACE_LIST(1,IFACE))
            CASE(CC_FTYPE_CFGAS)
               IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
               ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
               CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
            CASE(CC_FTYPE_RCGAS)
               II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
               CASE(CC_FTYPE_RGGAS) ! Regular cell
                  DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell
                  ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                  CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
               END SELECT
            END SELECT
         ENDDO
      ENDDO JCC1_LOOP
   ENDDO

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IF (RHOP(I,J,K)>=RHOMIN .AND. RHOP(I,J,K)<=RHOMAX) CYCLE
            IF (RHOP(I,J,K)<RHOMIN) THEN
               RHO_CUT = RHOMIN
               SIGN_FACTOR = 1._EB
               CLIP_RHOMIN_CC = .TRUE.
            ELSE
               RHO_CUT = RHOMAX
               SIGN_FACTOR = -1._EB
               CLIP_RHOMAX_CC = .TRUE.
            ENDIF
            MASS_C = ABS(RHO_CUT-RHOP(I,J,K)) * (DX(I)*DY(J)*DZ(K))
            ! Neighbor MASS_N contributions:
            MASS_N = 0._EB; NCELL=0
            DO X1AXIS=IAXIS,KAXIS
               DO LOHI=LOW_IND,HIGH_IND
                  ILH   = 2*LOHI - 3
                  IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
                  NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  IF (IRC>0) THEN ! Regular face and cell on the other side.
                     SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                     CASE(CC_FTYPE_RGGAS) ! Regular cell
                        RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                     CASE(CC_FTYPE_CFGAS) ! Cut-cell
                        ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                        RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
                        VOL(NCELL) = CUT_CELL(ICC1)%VOLUME(JCC1)
                     END SELECT
                  ELSE
                     RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                  ENDIF
                  MASS_N(NCELL) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHO_CELL(NCELL)))-RHO_CUT) * VOL(NCELL)
               ENDDO
            ENDDO
            ! We have MASS_C, MASS_N and NCELL, now distribute:
            SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
            CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
            DELTA_RHO(I,J,K) = DELTA_RHO(I,J,K) + CONST*SUM_MASS_N/(DX(I)*DY(J)*DZ(K))
            ! Neighbor cells:
            NCELL=0
            DO X1AXIS=IAXIS,KAXIS
               DO LOHI=LOW_IND,HIGH_IND
                  ILH   = 2*LOHI - 3
                  IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
                  NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  IF (IRC>0) THEN ! Regular face and cell on the other side.
                     SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                     CASE(CC_FTYPE_RGGAS) ! Regular cell
                        DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                     CASE(CC_FTYPE_CFGAS) ! Cut-cell
                        ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                        CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
                     END SELECT
                  ELSE
                     DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! 2. Assign DELTA_RHO to neighboring cells if clipping has been done.
   IF (CLIP_RHOMIN_CC .OR. CLIP_RHOMAX_CC) THEN
      IF (MESHES(NM)%PREDICTOR) THEN
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%RHOS(JCC) = MIN(RHOMAX,MAX(RHOMIN,CC%RHOS(JCC)+CC%DELTA_RHO(JCC))); ENDDO
         ENDDO
      ELSE
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%RHO(JCC)  = MIN(RHOMAX,MAX(RHOMIN,CC%RHO(JCC) +CC%DELTA_RHO(JCC))); ENDDO
         ENDDO
      ENDIF
      RHOP(1:IBAR,1:JBAR,1:KBAR) = MIN(RHOMAX,MAX(RHOMIN,RHOP(1:IBAR,1:JBAR,1:KBAR)+DELTA_RHO(1:IBAR,1:JBAR,1:KBAR)))
   ENDIF

   ! If there is only one gas species, set rho*Z=rho and return.
   CLIP_RHOMIN = CLIP_RHOMIN .OR. CLIP_RHOMIN_CC
   CLIP_RHOMAX = CLIP_RHOMAX .OR. CLIP_RHOMAX_CC

   IF (N_TRACKED_SPECIES==1) THEN
      IF (CLIP_RHOMIN_CC .OR. CLIP_RHOMAX_CC) THEN
         IF (MESHES(NM)%PREDICTOR) THEN
            DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
               CC => CUT_CELL(ICC)
               IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
               DO JCC=1,CC%NCELL; CC%ZZS(1,JCC) = CC%RHOS(JCC); ENDDO
            ENDDO
         ELSE
            DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
               CC => CUT_CELL(ICC)
               IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
               DO JCC=1,CC%NCELL; CC%ZZ(1,JCC)  = CC%RHO(JCC); ENDDO
            ENDDO
         ENDIF
         RHO_ZZ(1:IBAR,1:JBAR,1:KBAR,1) = RHOP(1:IBAR,1:JBAR,1:KBAR)
      ENDIF
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
               IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) &
                  RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      CYCLE MESH_LOOP
   ENDIF

   ! Correct species mass density
   ! Run through N species
   RHO_ZZ_MIN       =  0._EB
   DELTA_RHO_ZZ     => WORK5
   CLIP_RHO_ZZ_SAVE = .FALSE.
   SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

      DELTA_RHO_ZZ =  0._EB
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CUT_CELL(ICC)%DELTA_RHO_ZZ=0._EB
      ENDDO
      CLIP_RHO_ZZ  = .FALSE.

      ! 1. compute DELTA_RHO_ZZ in cut-cells and regular cells.
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IC=CELL_INDEX(I,J,K); IF (CELL(IC)%SOLID) CYCLE
         JCC2_LOOP : DO JCC=1,CC%NCELL
            CC_RHOP    = PRFCT*CC%RHOS(JCC)  + (1._EB-PRFCT)*CC%RHO(JCC); RHO_ZZ_MAX = CC_RHOP
            CC_RHO_ZZP = PRFCT*CC%ZZS(N,JCC) + (1._EB-PRFCT)*CC%ZZ(N,JCC)
            IF (CC_RHO_ZZP>=RHO_ZZ_MIN .AND. CC_RHO_ZZP<=RHO_ZZ_MAX) CYCLE JCC2_LOOP
            CLIP_RHO_ZZ = .TRUE.
            IF (CC_RHO_ZZP<RHO_ZZ_MIN) THEN
               RHO_ZZ_CUT = RHO_ZZ_MIN
               SIGN_FACTOR = 1._EB
            ELSE
               RHO_ZZ_CUT = RHO_ZZ_MAX
               SIGN_FACTOR = -1._EB
            ENDIF
            MASS_C = ABS(RHO_ZZ_CUT-CC_RHO_ZZP) * CC%VOLUME(JCC)

            ! Now find connected regular and cut-cells and add their contributions to Delta mass ZZ:
            MASS_N = 0._EB; NCELL=0
            DO IFC=2,CC%CCELEM(1,JCC)+1
               IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
               LOHI  = CC%FACE_LIST(2,IFACE)
               ILH   = 2*CC%FACE_LIST(2,IFACE) - 3 ! -1 for LOHI=LOW_IND, 1 for LOHI=HIGH_IND
               X1AXIS= CC%FACE_LIST(3,IFACE)
               IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
               NCELL=NCELL+1
               SELECT CASE(CC%FACE_LIST(1,IFACE))
               CASE(CC_FTYPE_CFGAS)
                  IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
                  ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
                  RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
                  VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
               CASE(CC_FTYPE_RCGAS)
                  II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                  CASE(CC_FTYPE_RGGAS) ! Regular cell
                     RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                  CASE(CC_FTYPE_CFGAS) ! Cut-cell
                     ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                     RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
                     VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
                  END SELECT
               END SELECT
               MASS_N(NCELL) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_CELL(NCELL)))-RHO_ZZ_CUT) * VOL(NCELL)
            ENDDO

            ! We have MASS_C, MASS_N and NCELL, now distribute:
            SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
            CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)

            CC%DELTA_RHO_ZZ(JCC) = CC%DELTA_RHO_ZZ(JCC) + CONST*SUM_MASS_N/CC%VOLUME(JCC)
            NCELL=0
            DO IFC=2,CC%CCELEM(1,JCC)+1
               IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
               LOHI  = CC%FACE_LIST(2,IFACE)
               ILH   = 2*CC%FACE_LIST(2,IFACE) - 3 ! -1 for LOHI=LOW_IND, 1 for LOHI=HIGH_IND
               X1AXIS= CC%FACE_LIST(3,IFACE)
               IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
               NCELL=NCELL+1
               SELECT CASE(CC%FACE_LIST(1,IFACE))
               CASE(CC_FTYPE_CFGAS)
                  IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
                  ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
                  CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
               CASE(CC_FTYPE_RCGAS)
                  II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                  CASE(CC_FTYPE_RGGAS) ! Regular cell
                     DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  CASE(CC_FTYPE_CFGAS) ! Cut-cell
                     ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                     CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  END SELECT
               END SELECT
            ENDDO
         ENDDO JCC2_LOOP
      ENDDO
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               RHO_ZZ_MAX = RHOP(I,J,K)
               IF (RHO_ZZ(I,J,K,N)>=RHO_ZZ_MIN .AND. RHO_ZZ(I,J,K,N)<=RHO_ZZ_MAX) CYCLE
               CLIP_RHO_ZZ = .TRUE.
               IF (RHO_ZZ(I,J,K,N)<RHO_ZZ_MIN) THEN
                  RHO_ZZ_CUT = RHO_ZZ_MIN
                  SIGN_FACTOR = 1._EB
               ELSE
                  RHO_ZZ_CUT = RHO_ZZ_MAX
                  SIGN_FACTOR = -1._EB
               ENDIF
               MASS_C = ABS(RHO_ZZ_CUT-RHO_ZZ(I,J,K,N)) * (DX(I)*DY(J)*DZ(K))
               ! Neighbor MASS_N contributions:
               MASS_N = 0._EB; NCELL=0
               DO X1AXIS=IAXIS,KAXIS
                  DO LOHI=LOW_IND,HIGH_IND
                     ILH   = 2*LOHI - 3
                     IF (CELL(CELL_INDEX(I,J,K))%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
                     NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                     SELECT CASE(X1AXIS)
                     CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                     CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                     CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                     END SELECT
                     IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                     IF (IRC>0) THEN ! Regular face and cell on the other side.
                        SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                        CASE(CC_FTYPE_RGGAS) ! Regular cell
                           RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                        CASE(CC_FTYPE_CFGAS) ! Cut-cell
                           ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                           RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
                           VOL(NCELL) = CUT_CELL(ICC1)%VOLUME(JCC1)
                        END SELECT
                     ELSE
                        RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                     ENDIF
                     MASS_N(NCELL) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_CELL(NCELL)))-RHO_ZZ_CUT) * VOL(NCELL)
                  ENDDO
               ENDDO
               ! We have MASS_C, MASS_N and NCELL, now distribute:
               SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
               CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
               DELTA_RHO_ZZ(I,J,K) = DELTA_RHO_ZZ(I,J,K) + CONST*SUM_MASS_N/(DX(I)*DY(J)*DZ(K))
               ! Neighbor cells:
               NCELL=0
               DO X1AXIS=IAXIS,KAXIS
                  DO LOHI=LOW_IND,HIGH_IND
                     ILH   = 2*LOHI - 3
                     IF (CELL(CELL_INDEX(I,J,K))%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE ! There is a wall cell here.
                     NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                     SELECT CASE(X1AXIS)
                     CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                     CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                     CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                     END SELECT
                     IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                     IF (IRC>0) THEN ! Regular face and cell on the other side.
                        SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                        CASE(CC_FTYPE_RGGAS) ! Regular cell
                           DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                        CASE(CC_FTYPE_CFGAS) ! Cut-cell
                           ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                           CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1)-CONST*MASS_N(NCELL)/VOL(NCELL)
                        END SELECT
                     ELSE
                        DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF (.NOT.CLIP_RHO_ZZ) THEN
         CYCLE SPECIES_LOOP
      ELSE
         CLIP_RHO_ZZ_SAVE = .TRUE.
      ENDIF

      ! 2. Assign excess/deficit RHO_ZZ neighboring cells
      IF (MESHES(NM)%PREDICTOR) THEN
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%ZZS(N,JCC) = MIN(CC%RHOS(JCC),MAX(RHO_ZZ_MIN,CC%ZZS(N,JCC)+CC%DELTA_RHO_ZZ(JCC))); ENDDO
         ENDDO
      ELSE
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%ZZ(N,JCC)  = MIN(CC%RHO(JCC), MAX(RHO_ZZ_MIN,CC%ZZ(N,JCC) +CC%DELTA_RHO_ZZ(JCC))); ENDDO
         ENDDO
      ENDIF
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               RHO_ZZ(I,J,K,N) = MIN(RHOP(I,J,K),MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K,N)+DELTA_RHO_ZZ(I,J,K)))
            ENDDO
         ENDDO
      ENDDO

   ENDDO SPECIES_LOOP

   ! If nothing has been clipped, return

   IF (.NOT.CLIP_RHOMIN_CC .AND. .NOT.CLIP_RHOMAX_CC .AND. .NOT.CLIP_RHO_ZZ_SAVE) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
               IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) &
                  RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      CYCLE MESH_LOOP
   ENDIF

   ! Final check of RHO_ZZ => SUM(RHO_ZZ) = RHO
   IF (MESHES(NM)%PREDICTOR) THEN
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            SUM_RHO_ZZ = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC))
            N = MAXLOC(CC%ZZS(1:N_TRACKED_SPECIES,JCC),1)
            RHO_ZZ_TEST = CC%ZZS(N,JCC) + CC%RHOS(JCC) - SUM_RHO_ZZ
            IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>CC%RHOS(JCC)) THEN  ! Renormalize the original set of RHO_ZZ
               CC%ZZS(1:N_TRACKED_SPECIES,JCC) = CC%RHOS(JCC) * CC%ZZS(1:N_TRACKED_SPECIES,JCC)/SUM_RHO_ZZ
            ELSE  ! Absorb mass deficit/excess into largest RHO_ZZ
               CC%ZZS(N,JCC) = RHO_ZZ_TEST
            ENDIF
         ENDDO
      ENDDO
   ELSE
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            SUM_RHO_ZZ = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC))
            N = MAXLOC(CC%ZZ(1:N_TRACKED_SPECIES,JCC),1)
            RHO_ZZ_TEST = CC%ZZ(N,JCC) + CC%RHO(JCC) - SUM_RHO_ZZ
            IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>CC%RHO(JCC)) THEN  ! Renormalize the original set of RHO_ZZ
               CC%ZZ(1:N_TRACKED_SPECIES,JCC) = CC%RHO(JCC) * CC%ZZ(1:N_TRACKED_SPECIES,JCC)/SUM_RHO_ZZ
            ELSE  ! Absorb mass deficit/excess into largest RHO_ZZ
               CC%ZZ(N,JCC) = RHO_ZZ_TEST
            ENDIF
         ENDDO
      ENDDO
   ENDIF
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE
            SUM_RHO_ZZ = SUM(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES))
            N = MAXLOC(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES),1)
            RHO_ZZ_TEST = RHO_ZZ(I,J,K,N) + RHOP(I,J,K) - SUM_RHO_ZZ
            IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>RHOP(I,J,K)) THEN  ! Renormalize the original set of RHO_ZZ
               RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHOP(I,J,K) * RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/SUM_RHO_ZZ
            ELSE  ! Absorb mass deficit/excess into largest RHO_ZZ
               RHO_ZZ(I,J,K,N) = RHO_ZZ_TEST
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   CALL CC_CV_RHOZZ_AVERAGE

   ! Bring back ZZ in regular gas cells from partial densities.
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE
            IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) &
               RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
CONTAINS

SUBROUTINE CC_CV_RHOZZ_AVERAGE

INTEGER :: IROW_LOC

! CV volumes:
RZ_ZS(:) = 0._EB
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + DX(I)*DY(J)*DZ(K) ! Here store volume in allocated RZ_ZS.
      ENDDO
   ENDDO
ENDDO
DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
   DO JCC=1,CC%NCELL
      IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
      RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + CC%VOLUME(JCC)
   ENDDO
ENDDO

! Loop species:
SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   RZ_Z(:) = 0._EB
   ! Add to CV rhoZZ*Vol:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            RZ_Z( IROW_LOC) = RZ_Z( IROW_LOC) + RHO_ZZ(I,J,K,N)*DX(I)*DY(J)*DZ(K) ! Known rho*zz
         ENDDO
      ENDDO
   ENDDO
   IF (MESHES(NM)%PREDICTOR) THEN
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%ZZS(N,JCC) * CC%VOLUME(JCC) ! Contains rho*zz
         ENDDO
      ENDDO
   ELSE
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%ZZ(N,JCC) * CC%VOLUME(JCC) ! Contains rho*zz
         ENDDO
      ENDDO
   ENDIF

   ! Volume average:
   DO IROW_LOC=UNKZ_IND(NM)-UNKZ_IND(NM_START)+1,UNKZ_IND(NM)-UNKZ_IND(NM_START)+NUNKZ_LOC(NM)
      RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) / RZ_ZS(IROW_LOC)
   ENDDO

   ! Back to cut/reg cell containers of rhoZZ:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            RHO_ZZ(I,J,K,N) = RZ_Z(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
         ENDDO
      ENDDO
   ENDDO
   IF (MESHES(NM)%PREDICTOR) THEN
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZS(N,JCC) = RZ_Z(CC%UNKZ(JCC)-UNKZ_IND(NM_START)); ENDDO
      ENDDO
   ELSE
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZ(N,JCC) = RZ_Z(CC%UNKZ(JCC)-UNKZ_IND(NM_START)); ENDDO
      ENDDO
   ENDIF

ENDDO SPECIES_LOOP

RZ_ZS = 0._EB

END SUBROUTINE CC_CV_RHOZZ_AVERAGE

END SUBROUTINE CC_CHECK_MASS_DENSITY



! ------------------------------ GET_RHOZZ_CC_3D ---------------------------

SUBROUTINE GET_RHOZZ_CC_3D

! Local Variables:
INTEGER :: NM,N,I,J,K,ICC,JCC
REAL(EB), POINTER, DIMENSION(:,:,:)   :: RHOP,UP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: VOLTOT
INTEGER :: NMX

! Loop meshes:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   IF (MESHES(NM)%PREDICTOR) THEN
      RHOP => RHOS
      ZZP  => ZZS
      UP   => U

      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            ! Get rho = sum(rho*z_alpha)
            CC%RHOS(JCC) = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC))

            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               ! Check mass density for positivity
               IF ( (CC%RHOS(JCC)<RHOMIN) .OR. (CC%RHOS(JCC)>RHOMAX) ) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Pred:',ICC,JCC,', Vol Fraction=',&
                  CC%VOLUME(JCC)/(DX(CC%IJK(IAXIS))*DY(CC%IJK(JAXIS))*DZ(CC%IJK(KAXIS))),', UNK_Z=',CC%UNKZ(JCC)
                  WRITE(LU_ERR,*) 'CELL Location=',XC(CC%IJK(IAXIS)),YC(CC%IJK(JAXIS)),ZC(CC%IJK(KAXIS))
                  WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',CC%RHOS(JCC),RHOMIN,RHOMAX
               ENDIF
            ENDIF

            ! Extract z from rho*z
            CC%ZZS(1:N_TOTAL_SCALARS,JCC) = CC%ZZS(1:N_TOTAL_SCALARS,JCC)/CC%RHOS(JCC)

            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               ! Check bounds on z:
               DO N=1,N_TOTAL_SCALARS
                  IF ( (CC%ZZS(N,JCC)<(0._EB-GEOMEPS)) .OR. (CC%ZZS(N,JCC)>(1._EB+GEOMEPS)) ) THEN
                     WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Pred:',ICC,JCC,N
                     WRITE(LU_ERR,*) 'ZZP=',CC%ZZS(N,JCC)
                  ENDIF
               ENDDO
            ELSE
               ! Some z_alpha might be slightly below zero (bounds overrun), assign -ve mass to most abundant species:
               ! Note rho = sum(rho*z_alpha), sum(z_alpha)=1 remain unchanged.
               NMX=MAXLOC(CC%ZZS(1:N_TRACKED_SPECIES,JCC),DIM=1)
               DO N=1,N_TRACKED_SPECIES
                  IF(N==NMX) CYCLE
                  IF ( CC%ZZS(N,JCC) < (0._EB-TWENTY_EPSILON_EB)) THEN
                     CC%ZZS(NMX,JCC) = CC%ZZS(NMX,JCC) + CC%ZZS(N,JCC)
                     CC%ZZS(N,JCC)   = 0._EB
                  ENDIF
               ENDDO
            ENDIF

            ! Clip passive scalars:
            IF (N_PASSIVE_SCALARS==0) CYCLE
            CC%ZZS(ZETA_INDEX,JCC) = MAX(0._EB,MIN(1._EB,CC%ZZS(ZETA_INDEX,JCC)))
         ENDDO

         ! Dump volume average scalar mass fraction and density to Cartesian container:
         I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
         VOLTOT = SUM( CC%VOLUME(1:CC%NCELL) )
         RHOP(I,J,K) = SUM( CC%RHOS(1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
         DO N=1,N_TOTAL_SCALARS
            ZZP(I,J,K,N) = SUM( CC%ZZS(N,1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
         ENDDO

      ENDDO

   ELSE
      RHOP => RHO
      ZZP  => ZZ
      UP   => US

      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            ! Get rho = sum(rho*z_alpha)
            CC%RHO(JCC) = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC))

            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               ! Check mass density for positivity
               IF ( (CC%RHO(JCC)<RHOMIN) .OR. (CC%RHO(JCC)>RHOMAX) ) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Corr:',ICC,JCC,', Vol Fraction=',&
                  CC%VOLUME(JCC)/(DX(CC%IJK(IAXIS))*DY(CC%IJK(JAXIS))*DZ(CC%IJK(KAXIS))),', UNK_Z=',CC%UNKZ(JCC)
                  WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',CC%RHO(JCC),RHOMIN,RHOMAX
               ENDIF
            ENDIF

            ! Extract z from rho*z
            CC%ZZ(1:N_TOTAL_SCALARS,JCC) = CC%ZZ(1:N_TOTAL_SCALARS,JCC)/CC%RHO(JCC)

            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               ! Check bounds on z:
               DO N=1,N_TOTAL_SCALARS
                  IF ( (CC%ZZ(N,JCC)<(0._EB-GEOMEPS)) .OR. (CC%ZZ(N,JCC)>(1._EB+GEOMEPS)) ) THEN
                     WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Corr:',ICC,JCC,N
                     WRITE(LU_ERR,*) 'ZZP=',CC%ZZ(N,JCC)
                  ENDIF
               ENDDO
            ELSE
               ! Some z_alpha might be slightly below zero (bounds overrun), assign -ve mass to most abundant species:
               ! Note rho = sum(rho*z_alpha), sum(z_alpha)=1 remain unchanged.
               NMX=MAXLOC(CC%ZZ(1:N_TRACKED_SPECIES,JCC),DIM=1)
               DO N=1,N_TRACKED_SPECIES
                  IF(N==NMX) CYCLE
                  IF ( CC%ZZ(N,JCC) < (0._EB-TWENTY_EPSILON_EB)) THEN
                     CC%ZZ(NMX,JCC) = CC%ZZ(NMX,JCC) + CC%ZZ(N,JCC)
                     CC%ZZ(N,JCC)   = 0._EB
                  ENDIF
               ENDDO
            ENDIF
            ! Clip passive scalars:
            IF (N_PASSIVE_SCALARS==0) CYCLE
            CC%ZZ(ZETA_INDEX,JCC) = MAX(0._EB,MIN(1._EB,CC%ZZ(ZETA_INDEX,JCC)))
         ENDDO

         ! Dump volume average scalar mass fraction and density to Cartesian container:
         I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
         VOLTOT = SUM( CC%VOLUME(1:CC%NCELL) )
         RHOP(I,J,K) = SUM( CC%RHO(1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
         DO N=1,N_TOTAL_SCALARS
            ZZP(I,J,K,N) = SUM( CC%ZZ(N,1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
         ENDDO

      ENDDO

   ENDIF

   ! Regular Cartesian cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MESHES(NM)%CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE ! Cycle Reg cells not in cc-region, cut-cells
                                                             ! underlying Cartesian cells and
                                                             ! solid cells.

            ! Get rho = sum(rho*z_alpha)
            RHOP(I,J,K) = SUM(ZZP(I,J,K,1:N_TRACKED_SPECIES))

            ! Check mass density for positivity
            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               IF ((RHOP(I,J,K)<RHOMIN) .OR. (RHOP(I,J,K)>RHOMAX) ) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D Cart:',I,J,K
                  WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',RHOP(I,J,K),RHOMIN,RHOMAX
               ENDIF
            ENDIF

            ! Extract z from rho*z
            ZZP(I,J,K,1:N_TOTAL_SCALARS) = ZZP(I,J,K,1:N_TOTAL_SCALARS)/RHOP(I,J,K)

            IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
               ! Check bounds on z:
               DO N=1,N_TOTAL_SCALARS
                  IF (ZZP(I,J,K,N)<(0._EB-GEOMEPS) .OR. ZZP(I,J,K,N)>(1._EB+GEOMEPS)) THEN
                     WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D Cart:',I,J,K,N
                     WRITE(LU_ERR,*) 'ZZP=',ZZP(I,J,K,N)
                  ENDIF
               ENDDO
            ELSE
               ! Some z_alpha might be slightly below zero (bounds overrun), assign -ve mass to most abundant species:
               ! Note rho = sum(rho*z_alpha), sum(z_alpha)=1 remain unchanged.
               NMX=MAXLOC(ZZP(I,J,K,1:N_TRACKED_SPECIES),DIM=1)
               DO N=1,N_TRACKED_SPECIES
                  IF(N==NMX) CYCLE
                  IF ( ZZP(I,J,K,N) < (0._EB-TWENTY_EPSILON_EB)) THEN
                     ZZP(I,J,K,NMX) = ZZP(I,J,K,NMX) + ZZP(I,J,K,N)
                     ZZP(I,J,K,N)   = 0._EB
                  ENDIF
               ENDDO
            ENDIF
            ! Clip passive scalars:
            IF (N_PASSIVE_SCALARS==0) CYCLE
            ZZP(I,J,K,ZETA_INDEX) = MAX(0._EB,MIN(1._EB,ZZP(I,J,K,ZETA_INDEX)))

         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
END SUBROUTINE GET_RHOZZ_CC_3D



! --------------------------- PUT_RHOZZVECTOR_SCALAR_3D --------------------------

SUBROUTINE PUT_RHOZZVECTOR_SCALAR_3D(N)

INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER :: NM,I,J,K,IROW_LOC,ICC,JCC
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

! Loop meshes:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   IF (MESHES(NM)%PREDICTOR) THEN
      ZZP  => ZZS ! Copy rho*z obtained for species N in the end of substep container for z.
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC      = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            CC%ZZS(N,JCC) = RZ_Z(IROW_LOC)
         ENDDO
      ENDDO
   ELSE
      ZZP  => ZZ  ! Copy rho*z obtained for species N in the end of substep container for z.
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC     = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            CC%ZZ(N,JCC) = RZ_Z(IROW_LOC)
         ENDDO
      ENDDO
   ENDIF

   ! Loop on Cartesian Cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC     = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            ZZP(I,J,K,N) = RZ_Z(IROW_LOC)
         ENDDO
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
END SUBROUTINE PUT_RHOZZVECTOR_SCALAR_3D



! --------------------------- GET_RHOZZVECTOR_SCALAR_3D --------------------------

SUBROUTINE GET_RHOZZVECTOR_SCALAR_3D(N)

INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER :: NM,I,J,K,IROW_LOC,ICC,JCC

! Initialize rho*z:
RZ_Z(:) = 0._EB
RZ_ZS(:) = 0._EB

! Loop meshes:
PREDCORR_IF : IF (MESHES(LOWER_MESH_INDEX)%PREDICTOR) THEN
   MESH_LOOP_P : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! Loop on Cartesian Cells:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
               RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + DX(I)*DY(J)*DZ(K) ! Here store volume in allocated RZ_ZS.
               RZ_Z( IROW_LOC) = RZ_Z( IROW_LOC) + RHO(I,J,K)*ZZ(I,J,K,N)*DX(I)*DY(J)*DZ(K) ! Known rho*zz^n
            ENDDO
         ENDDO
      ENDDO
      ! Now loop Cut-cells:
      CUTCELL_LOOP : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + CC%VOLUME(JCC)
            RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%RHO(JCC) * CC%ZZ(N,JCC) * CC%VOLUME(JCC)
         ENDDO
      ENDDO CUTCELL_LOOP
   ENDDO MESH_LOOP_P
   ! Volume average:
   DO IROW_LOC=1,NUNKZ_LOCAL
      RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) / RZ_ZS(IROW_LOC)
   ENDDO
   RZ_ZS = 0._EB

ELSE PREDCORR_IF

   MESH_LOOP_C : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! Loop on Cartesian Cells:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
               RZ_Z(IROW_LOC) = 0.5_EB*(MESHES(NM)%RHO_ZZN(I,J,K,N)+RHOS(I,J,K)*ZZS(I,J,K,N))
               RZ_ZS(IROW_LOC)= RHOS(I,J,K)*ZZS(I,J,K,N)
            ENDDO
         ENDDO
      ENDDO
      ! Now loop Cut-cells:
      CUTCELL_LOOP2 : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC) = 0.5_EB*(CC%RHO(JCC) * CC%ZZ(N,JCC) + CC%RHOS(JCC) * CC%ZZS(N,JCC))
            RZ_ZS(IROW_LOC)= CC%RHOS(JCC) * CC%ZZS(N,JCC)
         ENDDO
      ENDDO CUTCELL_LOOP2
   ENDDO MESH_LOOP_C

ENDIF PREDCORR_IF

RETURN
END SUBROUTINE GET_RHOZZVECTOR_SCALAR_3D



! ------------------------- GET_ADVDIFFVECTOR_SCALAR_3D -------------------------

SUBROUTINE GET_ADVDIFFVECTOR_SCALAR_3D(N)

INTEGER, INTENT(IN) :: N

! Local Variables:
INTEGER :: NM,I,J,K
REAL(EB):: PRFCT
INTEGER :: X1AXIS,IFACE,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),ICF,IND1,IND2,IOR
INTEGER :: LOCROW_1,LOCROW_2,ILOC,IROW,ICC,JCC,IW
REAL(EB):: AF,VELC,RHO_Z,FN_ZZ,FCT
REAL(EB), POINTER, DIMENSION(:,:,:)  :: RHOP,UP,VP,WP
REAL(EB), POINTER, DIMENSION(:,:,:,:)::  ZZP
TYPE(CC_REGFACEZ_TYPE),  POINTER, DIMENSION(:) :: REGFACE_Z
LOGICAL :: DO_LO,DO_HI

! This routine computes RHS due to boundary conditions prescribed in immersed solids
! and domain boundaries.

! First Domain Boundaries:
! Mesh Loop, Advective Fluxes:
MESH_LOOP_DBND : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   IF (MESHES(NM)%PREDICTOR) THEN
      ZZP  => ZZ
      RHOP => RHO
      UP   => U
      VP   => V
      WP   => W
      PRFCT= 1._EB
   ELSE
      ZZP  => ZZS
      RHOP => RHOS
      UP   => US
      VP   => VS
      WP   => WS
      PRFCT= 0._EB
   ENDIF

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.

   ! First add advective fluxes to domain boundary regular and cut-cells:
   ! IAXIS faces:
   X1AXIS = IAXIS
   REGFACE_Z => CC_REGFACE_IAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      I    = REGFACE_Z(IFACE)%IJK(IAXIS)
      J    = REGFACE_Z(IFACE)%IJK(JAXIS)
      K    = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO= REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI= REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I  ,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I+1,J,K,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DY(J)*DZ(K)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! JAXIS faces:
   X1AXIS = JAXIS
   REGFACE_Z => CC_REGFACE_JAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      I     = REGFACE_Z(IFACE)%IJK(IAXIS)
      J     = REGFACE_Z(IFACE)%IJK(JAXIS)
      K     = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I,J  ,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I,J+1,K,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DX(I)*DZ(K)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! KAXIS faces:
   X1AXIS = KAXIS
   REGFACE_Z => CC_REGFACE_KAXIS_Z
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      I     = REGFACE_Z(IFACE)%IJK(IAXIS)
      J     = REGFACE_Z(IFACE)%IJK(JAXIS)
      K     = REGFACE_Z(IFACE)%IJK(KAXIS)
      DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
      DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
      ! Unknowns on related cells:
      IND_LOC(LOW_IND) = CCVAR(I,J,K  ,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= CCVAR(I,J,K+1,CC_UNKZ) - UNKZ_IND(NM_START)

      AF = DX(I)*DY(J)
      IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
      IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   ENDDO

   ! Boundary Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   IFACE_LOOP_RCF1: DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z

      IW=MESHES(NM)%RC_FACE(IFACE)%IWC
      WC=>WALL(IW); IF ( WC%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE IFACE_LOOP_RCF1
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC => BOUNDARY_COORD(WC%BC_INDEX)

      I      = MESHES(NM)%RC_FACE(IFACE)%IJK(IAXIS)
      J      = MESHES(NM)%RC_FACE(IFACE)%IJK(JAXIS)
      K      = MESHES(NM)%RC_FACE(IFACE)%IJK(KAXIS)
      X1AXIS = MESHES(NM)%RC_FACE(IFACE)%IJK(KAXIS+1)

      ! Unknowns on related cells:
      IND(LOW_IND)  = MESHES(NM)%RC_FACE(IFACE)%UNKZ(LOW_IND)
      IND(HIGH_IND) = MESHES(NM)%RC_FACE(IFACE)%UNKZ(HIGH_IND)

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

      ! Row ind(1),ind(2):
      LOCROW_1 = LOW_IND
      LOCROW_2 = HIGH_IND

      IOR = BC%IOR
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ILOC=1,
      !                              when sign of IOR is  1 -> use High Side cell -> ILOC=2.
      ILOC = 1 + (SIGN(1,IOR)+1) / 2
      ! First (rho hs)_i,j,k:
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF = DY(J)*DZ(K)
         VELC = UP(I,J,K)
      CASE(JAXIS)
         AF = DX(I)*DZ(K)
         VELC = VP(I,J,K)
      CASE(KAXIS)
         AF = DX(I)*DY(J)
         VELC = WP(I,J,K)
      END SELECT
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N) ! These have been Flux limited in wall.f90
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT
            ! Already filled in previous X1AXIS select case.
         CASE(SOLID_BOUNDARY)
            IF (MESHES(NM)%PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            IF (MESHES(NM)%CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
         CASE(INTERPOLATED_BOUNDARY)
            VELC = UVW_SAVE(IW)
      END SELECT

      IROW=IND_LOC(ILOC)   ! Process Local Unknown number.
      FCT = REAL(3-2*ILOC,EB)
      F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF

   ENDDO IFACE_LOOP_RCF1

   ! Now Boundary Gasphase CUT_FACES:
   ICF_LOOP1: DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      IF ( MESHES(NM)%CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE ICF_LOOP1
      IW=MESHES(NM)%CUT_FACE(ICF)%IWC
      WC=>WALL(IW); IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE ICF_LOOP1
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC    => BOUNDARY_COORD(WC%BC_INDEX)
      IOR = BC%IOR
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N) ! These have been Flux limited in wall.f90
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ILOC=1,
      !                              when sign of IOR is  1 -> use High Side cell -> ILOC=2 .
      ILOC = 1 + (SIGN(1,IOR)+1) / 2
      FCT  = REAL(3-2*ILOC,EB)
      DO IFACE=1,CUT_FACE(ICF)%NFACE
         AF   = CUT_FACE(ICF)%AREA(IFACE)
         IF(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
            VELC = CUT_FACE(ICF)%VEL_SAVE(IFACE)
         ELSE
            VELC = PRFCT*CUT_FACE(ICF)%VEL(IFACE) + (1._EB-PRFCT)*CUT_FACE(ICF)%VELS(IFACE)
         ENDIF
         ! Unknowns on related cells:
         IND_LOC(ILOC) = MESHES(NM)%CUT_FACE(ICF)%UNKZ(ILOC,IFACE) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         ! First (rho hs)_i,j,k:
         IF (CUT_FACE(ICF)%CELL_LIST(1,ILOC,IFACE) == CC_FTYPE_CFGAS) THEN
            IROW=IND_LOC(ILOC)   ! Process Local Unknown number.
            F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
         ENDIF
      ENDDO ! IFACE

   ENDDO ICF_LOOP1

   ! INBOUNDARY cut-faces, loop on CFACE to add BC defined at SOLID phase:
   DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
      CFA  => CFACE(ICF)
      B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
      IND1 = CFA%CUT_FACE_IND1;                         IND2 = CFA%CUT_FACE_IND2
      ICC  = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2);  JCC  = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
      IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
      IF (MESHES(NM)%PREDICTOR) THEN
         VELC = B1%U_NORMAL
      ELSE
         VELC = B1%U_NORMAL_S
      ENDIF
      IF (VELC>0._EB) THEN
         RHO_Z = PRFCT *CUT_CELL(ICC)% RHO(JCC)*CUT_CELL(ICC)% ZZ(N,JCC) + &
          (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)*CUT_CELL(ICC)%ZZS(N,JCC)
      ELSE
         RHO_Z = B1%RHO_F*B1%ZZ_F(N)
      ENDIF
      F_Z(IROW) = F_Z(IROW) + RHO_Z*VELC*CFA%AREA
   ENDDO

   ! Then add diffusive fluxes through domain boundaries:
   ! Defined in CC_DIVERGENCE_PART_1.

ENDDO MESH_LOOP_DBND

RETURN
END SUBROUTINE GET_ADVDIFFVECTOR_SCALAR_3D


! ==================================================================================
! Thread-safe _TS variants for per-mesh parallelization (Hedgehog integration)
! ==================================================================================


! ------------------------------ GET_M_DOT_PPP_SCALAR_3D_TS ---------------------------

RECURSIVE SUBROUTINE GET_M_DOT_PPP_SCALAR_3D_TS(NM,M,N)

INTEGER, INTENT(IN) :: NM, N
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: I,J,K,IROW,ICC,JCC,NCELL

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: D_SOURCE, M_DOT_PPP_3D
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: M_DOT_PPP
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL

IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
CCVAR => M%CCVAR; CELL_INDEX => M%CELL_INDEX; CELL => M%CELL
DX => M%DX; DY => M%DY; DZ => M%DZ
CUT_CELL => M%CUT_CELL

! D_SOURCE(:,:,:) and M_DOT_PPP(:,:,:,:) are allocated together,
! => if one is not allocated the other also is not allocated.
IF(.NOT.ALLOCATED(M%D_SOURCE)) RETURN

D_SOURCE => M%D_SOURCE
M_DOT_PPP => M%M_DOT_PPP

! First add M_DOT_PPP on regular cells to source F_Z:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         F_Z(IROW) = F_Z(IROW) - M_DOT_PPP(I,J,K,N)*DX(I)*DY(J)*DZ(K)
      ENDDO
   ENDDO
ENDDO

! Then add Cut-cell contributions to F_Z:
DO ICC=1,M%N_CUTCELL_MESH
   IF (CELL(CELL_INDEX(CUT_CELL(ICC)%IJK(IAXIS),CUT_CELL(ICC)%IJK(JAXIS),CUT_CELL(ICC)%IJK(KAXIS)))%SOLID) CYCLE
   NCELL=CUT_CELL(ICC)%NCELL
   DO JCC=1,NCELL
      IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
      F_Z(IROW) = F_Z(IROW) - CUT_CELL(ICC)%M_DOT_PPP(N,JCC)*CUT_CELL(ICC)%VOLUME(JCC)
   ENDDO
ENDDO

! Finally if Corrector zero out M_DOT_PPP and D_SOURCE:
IF (M%CORRECTOR) THEN
   M_DOT_PPP(:,:,:,N) = 0._EB
   DO ICC=1,M%N_CUTCELL_MESH
      NCELL=CUT_CELL(ICC)%NCELL
      CUT_CELL(ICC)%M_DOT_PPP(N,1:NCELL) = 0._EB
   ENDDO
   IF (N == N_TOTAL_SCALARS) THEN
      D_SOURCE(:,:,:)  = 0._EB
      DO ICC=1,M%N_CUTCELL_MESH
         NCELL=CUT_CELL(ICC)%NCELL
         CUT_CELL(ICC)%D_SOURCE(1:NCELL) = 0._EB
      ENDDO
   ENDIF
ENDIF

RETURN
END SUBROUTINE GET_M_DOT_PPP_SCALAR_3D_TS


! --------------------------- PUT_RHOZZVECTOR_SCALAR_3D_TS --------------------------

RECURSIVE SUBROUTINE PUT_RHOZZVECTOR_SCALAR_3D_TS(NM,M,N)

INTEGER, INTENT(IN) :: NM, N
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: I,J,K,IROW_LOC,ICC,JCC

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTCELL_TYPE), POINTER :: CC

IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
CCVAR => M%CCVAR; CELL_INDEX => M%CELL_INDEX; CELL => M%CELL
CUT_CELL => M%CUT_CELL

IF (M%PREDICTOR) THEN
   ZZP  => M%ZZS
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW_LOC      = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         CC%ZZS(N,JCC) = RZ_Z(IROW_LOC)
      ENDDO
   ENDDO
ELSE
   ZZP  => M%ZZ
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW_LOC     = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         CC%ZZ(N,JCC) = RZ_Z(IROW_LOC)
      ENDDO
   ENDDO
ENDIF

! Loop on Cartesian Cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW_LOC     = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         ZZP(I,J,K,N) = RZ_Z(IROW_LOC)
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE PUT_RHOZZVECTOR_SCALAR_3D_TS


! --------------------------- GET_RHOZZVECTOR_SCALAR_3D_TS --------------------------

RECURSIVE SUBROUTINE GET_RHOZZVECTOR_SCALAR_3D_TS(NM,M,N)

INTEGER, INTENT(IN) :: NM, N
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: I,J,K,IROW_LOC,ICC,JCC
INTEGER :: ILC_LO, ILC_HI

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO, RHOS
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ, ZZS
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTCELL_TYPE), POINTER :: CC

IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
CCVAR => M%CCVAR; CELL_INDEX => M%CELL_INDEX; CELL => M%CELL
DX => M%DX; DY => M%DY; DZ => M%DZ
RHO => M%RHO; RHOS => M%RHOS; ZZ => M%ZZ; ZZS => M%ZZS
CUT_CELL => M%CUT_CELL

! Per-mesh range
ILC_LO = UNKZ_ILC(NM) + 1
ILC_HI = UNKZ_ILC(NM) + NUNKZ_LOC(NM)

! Initialize per-mesh range of rho*z:
RZ_Z(ILC_LO:ILC_HI) = 0._EB
RZ_ZS(ILC_LO:ILC_HI) = 0._EB

PREDCORR_IF : IF (M%PREDICTOR) THEN

   ! Loop on Cartesian Cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + DX(I)*DY(J)*DZ(K)
            RZ_Z( IROW_LOC) = RZ_Z( IROW_LOC) + RHO(I,J,K)*ZZ(I,J,K,N)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   ! Now loop Cut-cells:
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + CC%VOLUME(JCC)
         RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%RHO(JCC) * CC%ZZ(N,JCC) * CC%VOLUME(JCC)
      ENDDO
   ENDDO

   ! Volume average (per-mesh range only):
   DO IROW_LOC=ILC_LO,ILC_HI
      RZ_Z(IROW_LOC) = RZ_Z(IROW_LOC) / RZ_ZS(IROW_LOC)
   ENDDO
   RZ_ZS(ILC_LO:ILC_HI) = 0._EB

ELSE PREDCORR_IF

   ! Loop on Cartesian Cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC) = 0.5_EB*(M%RHO_ZZN(I,J,K,N)+RHOS(I,J,K)*ZZS(I,J,K,N))
            RZ_ZS(IROW_LOC)= RHOS(I,J,K)*ZZS(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO
   ! Now loop Cut-cells:
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         RZ_Z(IROW_LOC) = 0.5_EB*(CC%RHO(JCC) * CC%ZZ(N,JCC) + CC%RHOS(JCC) * CC%ZZS(N,JCC))
         RZ_ZS(IROW_LOC)= CC%RHOS(JCC) * CC%ZZS(N,JCC)
      ENDDO
   ENDDO

ENDIF PREDCORR_IF

RETURN
END SUBROUTINE GET_RHOZZVECTOR_SCALAR_3D_TS


! ------------------------------ GET_RHOZZ_CC_3D_TS ---------------------------

RECURSIVE SUBROUTINE GET_RHOZZ_CC_3D_TS(NM,M)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: N,I,J,K,ICC,JCC

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ, XC, YC, ZC
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTCELL_TYPE), POINTER :: CC
REAL(EB) :: VOLTOT
INTEGER :: NMX

IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
CCVAR => M%CCVAR; CELL_INDEX => M%CELL_INDEX; CELL => M%CELL
DX => M%DX; DY => M%DY; DZ => M%DZ
XC => M%XC; YC => M%YC; ZC => M%ZC
CUT_CELL => M%CUT_CELL

IF (M%PREDICTOR) THEN
   RHOP => M%RHOS
   ZZP  => M%ZZS

   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         CC%RHOS(JCC) = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC))

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            IF ( (CC%RHOS(JCC)<RHOMIN) .OR. (CC%RHOS(JCC)>RHOMAX) ) THEN
               WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Pred:',ICC,JCC,', Vol Fraction=',&
               CC%VOLUME(JCC)/(DX(CC%IJK(IAXIS))*DY(CC%IJK(JAXIS))*DZ(CC%IJK(KAXIS))),', UNK_Z=',CC%UNKZ(JCC)
               WRITE(LU_ERR,*) 'CELL Location=',XC(CC%IJK(IAXIS)),YC(CC%IJK(JAXIS)),ZC(CC%IJK(KAXIS))
               WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',CC%RHOS(JCC),RHOMIN,RHOMAX
            ENDIF
         ENDIF

         CC%ZZS(1:N_TOTAL_SCALARS,JCC) = CC%ZZS(1:N_TOTAL_SCALARS,JCC)/CC%RHOS(JCC)

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            DO N=1,N_TOTAL_SCALARS
               IF ( (CC%ZZS(N,JCC)<(0._EB-GEOMEPS)) .OR. (CC%ZZS(N,JCC)>(1._EB+GEOMEPS)) ) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Pred:',ICC,JCC,N
                  WRITE(LU_ERR,*) 'ZZP=',CC%ZZS(N,JCC)
               ENDIF
            ENDDO
         ELSE
            NMX=MAXLOC(CC%ZZS(1:N_TRACKED_SPECIES,JCC),DIM=1)
            DO N=1,N_TRACKED_SPECIES
               IF(N==NMX) CYCLE
               IF ( CC%ZZS(N,JCC) < (0._EB-TWENTY_EPSILON_EB)) THEN
                  CC%ZZS(NMX,JCC) = CC%ZZS(NMX,JCC) + CC%ZZS(N,JCC)
                  CC%ZZS(N,JCC)   = 0._EB
               ENDIF
            ENDDO
         ENDIF

         IF (N_PASSIVE_SCALARS==0) CYCLE
         CC%ZZS(ZETA_INDEX,JCC) = MAX(0._EB,MIN(1._EB,CC%ZZS(ZETA_INDEX,JCC)))
      ENDDO

      I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      VOLTOT = SUM( CC%VOLUME(1:CC%NCELL) )
      RHOP(I,J,K) = SUM( CC%RHOS(1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
      DO N=1,N_TOTAL_SCALARS
         ZZP(I,J,K,N) = SUM( CC%ZZS(N,1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
      ENDDO
   ENDDO

ELSE
   RHOP => M%RHO
   ZZP  => M%ZZ

   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         CC%RHO(JCC) = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC))

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            IF ( (CC%RHO(JCC)<RHOMIN) .OR. (CC%RHO(JCC)>RHOMAX) ) THEN
               WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Corr:',ICC,JCC,', Vol Fraction=',&
               CC%VOLUME(JCC)/(DX(CC%IJK(IAXIS))*DY(CC%IJK(JAXIS))*DZ(CC%IJK(KAXIS))),', UNK_Z=',CC%UNKZ(JCC)
               WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',CC%RHO(JCC),RHOMIN,RHOMAX
            ENDIF
         ENDIF

         CC%ZZ(1:N_TOTAL_SCALARS,JCC) = CC%ZZ(1:N_TOTAL_SCALARS,JCC)/CC%RHO(JCC)

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            DO N=1,N_TOTAL_SCALARS
               IF ( (CC%ZZ(N,JCC)<(0._EB-GEOMEPS)) .OR. (CC%ZZ(N,JCC)>(1._EB+GEOMEPS)) ) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D CC Corr:',ICC,JCC,N
                  WRITE(LU_ERR,*) 'ZZP=',CC%ZZ(N,JCC)
               ENDIF
            ENDDO
         ELSE
            NMX=MAXLOC(CC%ZZ(1:N_TRACKED_SPECIES,JCC),DIM=1)
            DO N=1,N_TRACKED_SPECIES
               IF(N==NMX) CYCLE
               IF ( CC%ZZ(N,JCC) < (0._EB-TWENTY_EPSILON_EB)) THEN
                  CC%ZZ(NMX,JCC) = CC%ZZ(NMX,JCC) + CC%ZZ(N,JCC)
                  CC%ZZ(N,JCC)   = 0._EB
               ENDIF
            ENDDO
         ENDIF
         IF (N_PASSIVE_SCALARS==0) CYCLE
         CC%ZZ(ZETA_INDEX,JCC) = MAX(0._EB,MIN(1._EB,CC%ZZ(ZETA_INDEX,JCC)))
      ENDDO

      I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      VOLTOT = SUM( CC%VOLUME(1:CC%NCELL) )
      RHOP(I,J,K) = SUM( CC%RHO(1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
      DO N=1,N_TOTAL_SCALARS
         ZZP(I,J,K,N) = SUM( CC%ZZ(N,1:CC%NCELL)*CC%VOLUME(1:CC%NCELL) )/VOLTOT
      ENDDO
   ENDDO

ENDIF

! Regular Cartesian cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE

         RHOP(I,J,K) = SUM(ZZP(I,J,K,1:N_TRACKED_SPECIES))

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            IF ((RHOP(I,J,K)<RHOMIN) .OR. (RHOP(I,J,K)>RHOMAX) ) THEN
               WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D Cart:',I,J,K
               WRITE(LU_ERR,*) 'RHOP,MIN,MAX=',RHOP(I,J,K),RHOMIN,RHOMAX
            ENDIF
         ENDIF

         ZZP(I,J,K,1:N_TOTAL_SCALARS) = ZZP(I,J,K,1:N_TOTAL_SCALARS)/RHOP(I,J,K)

         IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
            DO N=1,N_TOTAL_SCALARS
               IF (ZZP(I,J,K,N)<(0._EB-GEOMEPS) .OR. ZZP(I,J,K,N)>(1._EB+GEOMEPS)) THEN
                  WRITE(LU_ERR,*) 'GET_RHOZZ_CC_3D Cart:',I,J,K,N
                  WRITE(LU_ERR,*) 'ZZP=',ZZP(I,J,K,N)
               ENDIF
            ENDDO
         ELSE
            NMX=MAXLOC(ZZP(I,J,K,1:N_TRACKED_SPECIES),DIM=1)
            DO N=1,N_TRACKED_SPECIES
               IF(N==NMX) CYCLE
               IF ( ZZP(I,J,K,N) < (0._EB-TWENTY_EPSILON_EB)) THEN
                  ZZP(I,J,K,NMX) = ZZP(I,J,K,NMX) + ZZP(I,J,K,N)
                  ZZP(I,J,K,N)   = 0._EB
               ENDIF
            ENDDO
         ENDIF
         IF (N_PASSIVE_SCALARS==0) CYCLE
         ZZP(I,J,K,ZETA_INDEX) = MAX(0._EB,MIN(1._EB,ZZP(I,J,K,ZETA_INDEX)))
      ENDDO
   ENDDO
ENDDO

RETURN
END SUBROUTINE GET_RHOZZ_CC_3D_TS


! ---------------------- GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS --------------------

RECURSIVE SUBROUTINE GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS(NM,M,N)

INTEGER, INTENT(IN) :: NM, N
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: I,J,K
REAL(EB):: PRFCT
INTEGER :: X1AXIS,IFACE,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),ICF
INTEGER :: LOCROW_1,LOCROW_2,ILOC,IROW,ICC,JCC,ISIDE,IW
REAL(EB):: AF,VELC,RHO_Z_PV(-2:1),RHOPV(-2:1),FCT,ZZ_GET_N,FN_ZZ
REAL(EB), POINTER, DIMENSION(:,:,:)  :: RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:):: ZZP
LOGICAL :: DO_LO,DO_HI
INTEGER :: IIG,JJG,KKG,IOR
REAL(EB) :: UN
INTEGER :: ILO_FACE_L, IHI_FACE_L, JLO_FACE_L, JHI_FACE_L, KLO_FACE_L, KHI_FACE_L

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER :: N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ, UVW_SAVE
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO, RHOS
REAL(EB), POINTER, DIMENSION(:,:,:) :: U, US, V, VS, W, WS
REAL(EB), POINTER, DIMENSION(:,:,:) :: WORK_U, WORK_V, WORK_W
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: DEL_RHO_D_DEL_Z, ZZ, ZZS
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
TYPE(WALL_TYPE), POINTER, DIMENSION(:) :: WALL
TYPE(BOUNDARY_COORD_TYPE), POINTER, DIMENSION(:) :: BOUNDARY_COORD
TYPE(BOUNDARY_PROP1_TYPE), POINTER, DIMENSION(:) :: BOUNDARY_PROP1
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTFACE_TYPE), POINTER, DIMENSION(:) :: CUT_FACE
TYPE(CC_RCFACE_TYPE), POINTER, DIMENSION(:) :: RC_FACE
TYPE(CC_REGFACEZ_TYPE), POINTER, DIMENSION(:) :: CC_REGFACE_IAXIS_Z, CC_REGFACE_JAXIS_Z, CC_REGFACE_KAXIS_Z
! CC_SCALARS_DATA convenience pointer shadows
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
! Workspace shadows (thread-local)
REAL(EB), TARGET, DIMENSION(0:3,0:3,0:3) :: F_WORK
REAL(EB), TARGET, DIMENSION(-1:3,-1:3,-1:3) :: U_WORK, Z_WORK
REAL(EB), POINTER, DIMENSION(:,:,:) :: U_TEMP, Z_TEMP, F_TEMP
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU, VV, WW, UP, VP, WP
TYPE(CC_REGFACEZ_TYPE), POINTER, DIMENSION(:) :: REGFACE_Z

! Thread-safe alias setup
IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
N_EXTERNAL_WALL_CELLS => M%N_EXTERNAL_WALL_CELLS
N_INTERNAL_WALL_CELLS => M%N_INTERNAL_WALL_CELLS
CELL_INDEX => M%CELL_INDEX; CCVAR => M%CCVAR
DX => M%DX; DY => M%DY; DZ => M%DZ
UVW_SAVE => M%UVW_SAVE
RHO => M%RHO; RHOS => M%RHOS
U => M%U; US => M%US; V => M%V; VS => M%VS; W => M%W; WS => M%WS
WORK_U => M%WORK_U; WORK_V => M%WORK_V; WORK_W => M%WORK_W
DEL_RHO_D_DEL_Z => M%DEL_RHO_D_DEL_Z
ZZ => M%ZZ; ZZS => M%ZZS
CELL => M%CELL; WALL => M%WALL
BOUNDARY_COORD => M%BOUNDARY_COORD; BOUNDARY_PROP1 => M%BOUNDARY_PROP1
CUT_CELL => M%CUT_CELL; CUT_FACE => M%CUT_FACE; RC_FACE => M%RC_FACE
CC_REGFACE_IAXIS_Z => M%CC_REGFACE_IAXIS_Z
CC_REGFACE_JAXIS_Z => M%CC_REGFACE_JAXIS_Z
CC_REGFACE_KAXIS_Z => M%CC_REGFACE_KAXIS_Z

! Initialize thread-local workspace pointers
U_TEMP => U_WORK
F_TEMP => F_WORK
Z_TEMP => Z_WORK

UU=>WORK_U
VV=>WORK_V
WW=>WORK_W

IF (M%PREDICTOR) THEN
   ZZP  => ZZ
   RHOP => RHO
   UU   = U
   VV   = V
   WW   = W
   PRFCT= 1._EB
   WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      IOR = BC%IOR
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT; CYCLE WALL_LOOP
         CASE(SOLID_BOUNDARY);        UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
         CASE(INTERPOLATED_BOUNDARY); UN = UVW_SAVE(IW)
      END SELECT
      SELECT CASE(IOR)
         CASE( 1); UU(IIG-1,JJG,KKG) = UN
         CASE(-1); UU(IIG,JJG,KKG)   = UN
         CASE( 2); VV(IIG,JJG-1,KKG) = UN
         CASE(-2); VV(IIG,JJG,KKG)   = UN
         CASE( 3); WW(IIG,JJG,KKG-1) = UN
         CASE(-3); WW(IIG,JJG,KKG)   = UN
      END SELECT
   ENDDO WALL_LOOP

ELSE
   ZZP  => ZZS
   RHOP => RHOS
   UU   = US
   VV   = VS
   WW   = WS
   PRFCT= 0._EB
   WALL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      IOR = BC%IOR
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT; CYCLE WALL_LOOP_2
         CASE(SOLID_BOUNDARY);        UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
         CASE(INTERPOLATED_BOUNDARY); UN = UVW_SAVE(IW)
      END SELECT
      SELECT CASE(IOR)
         CASE( 1); UU(IIG-1,JJG,KKG) = UN
         CASE(-1); UU(IIG,JJG,KKG)   = UN
         CASE( 2); VV(IIG,JJG-1,KKG) = UN
         CASE(-2); VV(IIG,JJG,KKG)   = UN
         CASE( 3); WW(IIG,JJG,KKG-1) = UN
         CASE(-3); WW(IIG,JJG,KKG)   = UN
      END SELECT
   ENDDO WALL_LOOP_2
ENDIF

UP => UU
VP => VV
WP => WW

! Face bounds (local, not module SAVE)
ILO_FACE_L = 0;    IHI_FACE_L = IBAR
JLO_FACE_L = 0;    JHI_FACE_L = JBAR
KLO_FACE_L = 0;    KHI_FACE_L = KBAR

! First add advective fluxes to internal and INTERPOLATED_BOUNDARY regular and cut-cells in the CC region:
! IAXIS faces:
X1AXIS = IAXIS
REGFACE_Z => CC_REGFACE_IAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF (IW>0) THEN
      IF (.NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   ENDIF
   I     = REGFACE_Z(IFACE)%IJK(IAXIS)
   J     = REGFACE_Z(IFACE)%IJK(JAXIS)
   K     = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I  ,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I+1,J,K,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DY(J)*DZ(K)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
REGFACE_Z => CC_REGFACE_JAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I     = REGFACE_Z(IFACE)%IJK(IAXIS)
   J     = REGFACE_Z(IFACE)%IJK(JAXIS)
   K     = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I,J  ,K,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I,J+1,K,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DX(I)*DZ(K)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
REGFACE_Z => CC_REGFACE_KAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I     = REGFACE_Z(IFACE)%IJK(IAXIS)
   J     = REGFACE_Z(IFACE)%IJK(JAXIS)
   K     = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I,J,K  ,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I,J,K+1,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DX(I)*DY(J)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
DO IFACE=1,M%CC_NRCFACE_Z

   IW = M%RC_FACE(IFACE)%IWC; IF(IW > 0) CYCLE

   I      = M%RC_FACE(IFACE)%IJK(IAXIS)
   J      = M%RC_FACE(IFACE)%IJK(JAXIS)
   K      = M%RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = M%RC_FACE(IFACE)%IJK(KAXIS+1)

   IND(LOW_IND)  = M%RC_FACE(IFACE)%UNKZ(LOW_IND)
   IND(HIGH_IND) = M%RC_FACE(IFACE)%UNKZ(HIGH_IND)

   IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

   LOCROW_1 = LOW_IND
   LOCROW_2 = HIGH_IND
   SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF = DY(J)*DZ(K)
         RHOPV(-2:1)      = RHOP(I-1:I+2,J,K)
         DO ISIDE=-1,0
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS)
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
            CASE(CC_FTYPE_CFGAS)
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDDO
         ISIDE=-2
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         VELC = UU(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_ZZ = F_TEMP(1,1,1)
      CASE(JAXIS)
         AF = DX(I)*DZ(K)
         RHOPV(-2:1)      = RHOP(I,J-1:J+2,K)
         DO ISIDE=-1,0
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS)
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
            CASE(CC_FTYPE_CFGAS)
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDDO
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         VELC = VV(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_ZZ = F_TEMP(1,1,1)
      CASE(KAXIS)
         AF = DX(I)*DY(J)
         RHOPV(-2:1)      = RHOP(I,J,K-1:K+2)
         DO ISIDE=-1,0
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS)
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
            CASE(CC_FTYPE_CFGAS)
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDDO
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDIF
         VELC = WW(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_ZZ = F_TEMP(1,1,1)
   END SELECT

   DO ILOC=LOCROW_1,LOCROW_2
      IROW=IND_LOC(ILOC)
      FCT = REAL(3-2*ILOC,EB)
      F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
   ENDDO

ENDDO

! Now Gasphase CUT_FACES:
DO ICF = 1,M%N_CUTFACE_MESH

   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE

   IW = CUT_FACE(ICF)%IWC; IF(IW > 0) CYCLE

   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)

   LOCROW_1 = LOW_IND
   LOCROW_2 = HIGH_IND
   DO IFACE=1,CUT_FACE(ICF)%NFACE

      IND(LOW_IND)  = CUT_FACE(ICF)%UNKZ(LOW_IND,IFACE)
      IND(HIGH_IND) = CUT_FACE(ICF)%UNKZ(HIGH_IND,IFACE)

      IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START)
      IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

      AF = CUT_FACE(ICF)%AREA(IFACE)

      VELC =        PRFCT *CUT_FACE(ICF)%VEL(IFACE) + &
             (1._EB-PRFCT)*CUT_FACE(ICF)%VELS(IFACE)

      RHOPV(-1:0)    = -1._EB
      RHO_Z_PV(-1:0) =  0._EB
      DO ISIDE=-1,0
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS)
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
            ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
         END SELECT
         RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
      ENDDO
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         ISIDE=-2
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
            RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
            RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
         ENDIF
      CASE(JAXIS)
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
            RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
            RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
         ENDIF
      CASE(KAXIS)
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1)
         ELSE
            ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
            RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
         ENDIF
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1)
         ELSE
            ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
            RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
         ENDIF
      END SELECT
      VELC  = PRFCT *CUT_FACE(ICF)%VEL(IFACE) + (1._EB-PRFCT)*CUT_FACE(ICF)%VELS(IFACE)
      Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
      U_TEMP(1,1,1) = VELC
      CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
      FN_ZZ = F_TEMP(1,1,1)

      DO ILOC=LOCROW_1,LOCROW_2
         IROW=IND_LOC(ILOC)
         FCT = REAL(3-2*ILOC,EB)
         F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
      ENDDO

   ENDDO

ENDDO

! Then add (Del rho D Del Z)*dv computed on CCDIVERGENCE_PART_1:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW  = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         F_Z(IROW) = F_Z(IROW) - DEL_RHO_D_DEL_Z(I,J,K,N)*(DX(I)*DY(J)*DZ(K))
      ENDDO
   ENDDO
ENDDO

! Now cut-cells:
DO ICC=1,M%N_CUTCELL_MESH
   I = CUT_CELL(ICC)%IJK(IAXIS)
   J = CUT_CELL(ICC)%IJK(JAXIS)
   K = CUT_CELL(ICC)%IJK(KAXIS)
   IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
   DO JCC=1,CUT_CELL(ICC)%NCELL
      IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
      F_Z(IROW) = F_Z(IROW) - CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(N,JCC)
   ENDDO
ENDDO

RETURN
END SUBROUTINE GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS


! ------------------------- GET_ADVDIFFVECTOR_SCALAR_3D_TS -------------------------

RECURSIVE SUBROUTINE GET_ADVDIFFVECTOR_SCALAR_3D_TS(NM,M,N)

INTEGER, INTENT(IN) :: NM, N
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

! Local Variables:
INTEGER :: I,J,K
REAL(EB):: PRFCT
INTEGER :: X1AXIS,IFACE,IND(LOW_IND:HIGH_IND),IND_LOC(LOW_IND:HIGH_IND),ICF,IND1,IND2,IOR
INTEGER :: LOCROW_1,LOCROW_2,ILOC,IROW,ICC,JCC,IW
REAL(EB):: AF,VELC,RHO_Z,FN_ZZ,FCT
LOGICAL :: DO_LO,DO_HI

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER :: INTERNAL_CFACE_CELLS_LB, N_INTERNAL_CFACE_CELLS
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ, UVW_SAVE
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP, RHO, RHOS
REAL(EB), POINTER, DIMENSION(:,:,:) :: UP, VP, WP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP, ZZ, ZZS
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
TYPE(WALL_TYPE), POINTER, DIMENSION(:) :: WALL
TYPE(BOUNDARY_COORD_TYPE), POINTER, DIMENSION(:) :: BOUNDARY_COORD
TYPE(BOUNDARY_PROP1_TYPE), POINTER, DIMENSION(:) :: BOUNDARY_PROP1
TYPE(CFACE_TYPE), POINTER, DIMENSION(:) :: CFACE
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTFACE_TYPE), POINTER, DIMENSION(:) :: CUT_FACE
TYPE(CC_RCFACE_TYPE), POINTER, DIMENSION(:) :: RC_FACE
TYPE(CC_REGFACEZ_TYPE), POINTER, DIMENSION(:) :: CC_REGFACE_IAXIS_Z, CC_REGFACE_JAXIS_Z, CC_REGFACE_KAXIS_Z
! CC_SCALARS_DATA convenience pointer shadows
TYPE(WALL_TYPE), POINTER :: WC
TYPE(CFACE_TYPE), POINTER :: CFA
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(CC_REGFACEZ_TYPE), POINTER, DIMENSION(:) :: REGFACE_Z

! Thread-safe alias setup
IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
INTERNAL_CFACE_CELLS_LB => M%INTERNAL_CFACE_CELLS_LB
N_INTERNAL_CFACE_CELLS => M%N_INTERNAL_CFACE_CELLS
CELL_INDEX => M%CELL_INDEX; CCVAR => M%CCVAR
DX => M%DX; DY => M%DY; DZ => M%DZ
UVW_SAVE => M%UVW_SAVE
RHO => M%RHO; RHOS => M%RHOS
ZZ => M%ZZ; ZZS => M%ZZS
CELL => M%CELL; WALL => M%WALL
BOUNDARY_COORD => M%BOUNDARY_COORD; BOUNDARY_PROP1 => M%BOUNDARY_PROP1
CFACE => M%CFACE
CUT_CELL => M%CUT_CELL; CUT_FACE => M%CUT_FACE; RC_FACE => M%RC_FACE
CC_REGFACE_IAXIS_Z => M%CC_REGFACE_IAXIS_Z
CC_REGFACE_JAXIS_Z => M%CC_REGFACE_JAXIS_Z
CC_REGFACE_KAXIS_Z => M%CC_REGFACE_KAXIS_Z

IF (M%PREDICTOR) THEN
   ZZP  => ZZ
   RHOP => RHO
   UP   => M%U
   VP   => M%V
   WP   => M%W
   PRFCT= 1._EB
ELSE
   ZZP  => ZZS
   RHOP => RHOS
   UP   => M%US
   VP   => M%VS
   WP   => M%WS
   PRFCT= 0._EB
ENDIF

! First add advective fluxes to domain boundary regular and cut-cells:
! IAXIS faces:
X1AXIS = IAXIS
REGFACE_Z => CC_REGFACE_IAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I    = REGFACE_Z(IFACE)%IJK(IAXIS)
   J    = REGFACE_Z(IFACE)%IJK(JAXIS)
   K    = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO= REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI= REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I  ,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I+1,J,K,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DY(J)*DZ(K)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
REGFACE_Z => CC_REGFACE_JAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I     = REGFACE_Z(IFACE)%IJK(IAXIS)
   J     = REGFACE_Z(IFACE)%IJK(JAXIS)
   K     = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I,J  ,K,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I,J+1,K,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DX(I)*DZ(K)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
REGFACE_Z => CC_REGFACE_KAXIS_Z
DO IFACE=1,M%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WC%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I     = REGFACE_Z(IFACE)%IJK(IAXIS)
   J     = REGFACE_Z(IFACE)%IJK(JAXIS)
   K     = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO = REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI = REGFACE_Z(IFACE)%DO_HI_IND
   IND_LOC(LOW_IND) = CCVAR(I,J,K  ,CC_UNKZ) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= CCVAR(I,J,K+1,CC_UNKZ) - UNKZ_IND(NM_START)

   AF = DX(I)*DY(J)
   IF (DO_LO) F_Z(IND_LOC( LOW_IND)) = F_Z(IND_LOC( LOW_IND)) + REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
   IF (DO_HI) F_Z(IND_LOC(HIGH_IND)) = F_Z(IND_LOC(HIGH_IND)) - REGFACE_Z(IFACE)%RHOZZ_U(N)*AF
ENDDO

! Boundary Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
IFACE_LOOP_RCF1: DO IFACE=1,M%CC_NRCFACE_Z

   IW=M%RC_FACE(IFACE)%IWC
   WC=>WALL(IW); IF ( WC%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE IFACE_LOOP_RCF1
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)

   I      = M%RC_FACE(IFACE)%IJK(IAXIS)
   J      = M%RC_FACE(IFACE)%IJK(JAXIS)
   K      = M%RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = M%RC_FACE(IFACE)%IJK(KAXIS+1)

   IND(LOW_IND)  = M%RC_FACE(IFACE)%UNKZ(LOW_IND)
   IND(HIGH_IND) = M%RC_FACE(IFACE)%UNKZ(HIGH_IND)

   IND_LOC(LOW_IND) = IND(LOW_IND) - UNKZ_IND(NM_START)
   IND_LOC(HIGH_IND)= IND(HIGH_IND)- UNKZ_IND(NM_START)

   LOCROW_1 = LOW_IND
   LOCROW_2 = HIGH_IND

   IOR = BC%IOR
   ILOC = 1 + (SIGN(1,IOR)+1) / 2
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      AF = DY(J)*DZ(K)
      VELC = UP(I,J,K)
   CASE(JAXIS)
      AF = DX(I)*DZ(K)
      VELC = VP(I,J,K)
   CASE(KAXIS)
      AF = DX(I)*DY(J)
      VELC = WP(I,J,K)
   END SELECT
   FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
   SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
      CASE(SOLID_BOUNDARY)
         IF (M%PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
         IF (M%CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      CASE(INTERPOLATED_BOUNDARY)
         VELC = UVW_SAVE(IW)
   END SELECT

   IROW=IND_LOC(ILOC)
   FCT = REAL(3-2*ILOC,EB)
   F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF

ENDDO IFACE_LOOP_RCF1

! Now Boundary Gasphase CUT_FACES:
ICF_LOOP1: DO ICF = 1,M%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE ICF_LOOP1
   IW=CUT_FACE(ICF)%IWC
   WC=>WALL(IW); IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE ICF_LOOP1
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC    => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
   ILOC = 1 + (SIGN(1,IOR)+1) / 2
   FCT  = REAL(3-2*ILOC,EB)
   DO IFACE=1,CUT_FACE(ICF)%NFACE
      AF   = CUT_FACE(ICF)%AREA(IFACE)
      IF(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
         VELC = CUT_FACE(ICF)%VEL_SAVE(IFACE)
      ELSE
         VELC = PRFCT*CUT_FACE(ICF)%VEL(IFACE) + (1._EB-PRFCT)*CUT_FACE(ICF)%VELS(IFACE)
      ENDIF
      IND_LOC(ILOC) = CUT_FACE(ICF)%UNKZ(ILOC,IFACE) - UNKZ_IND(NM_START)
      IF (CUT_FACE(ICF)%CELL_LIST(1,ILOC,IFACE) == CC_FTYPE_CFGAS) THEN
         IROW=IND_LOC(ILOC)
         F_Z(IROW) = F_Z(IROW) + FCT*FN_ZZ*VELC*AF
      ENDIF
   ENDDO

ENDDO ICF_LOOP1

! INBOUNDARY cut-faces, loop on CFACE to add BC defined at SOLID phase:
DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA  => CFACE(ICF)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   IND1 = CFA%CUT_FACE_IND1;                         IND2 = CFA%CUT_FACE_IND2
   ICC  = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2);  JCC  = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
   IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
   IF (M%PREDICTOR) THEN
      VELC = B1%U_NORMAL
   ELSE
      VELC = B1%U_NORMAL_S
   ENDIF
   IF (VELC>0._EB) THEN
      RHO_Z = PRFCT *CUT_CELL(ICC)% RHO(JCC)*CUT_CELL(ICC)% ZZ(N,JCC) + &
       (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)*CUT_CELL(ICC)%ZZS(N,JCC)
   ELSE
      RHO_Z = B1%RHO_F*B1%ZZ_F(N)
   ENDIF
   F_Z(IROW) = F_Z(IROW) + RHO_Z*VELC*CFA%AREA
ENDDO

RETURN
END SUBROUTINE GET_ADVDIFFVECTOR_SCALAR_3D_TS


! ------------------------ CC_CHECK_MASS_DENSITY_TS ------------------------

RECURSIVE SUBROUTINE CC_CHECK_MASS_DENSITY_TS(NM,M)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

INTEGER, PARAMETER :: MAX_SURR_CELLS=20
INTEGER :: NCELL, ICC, JCC, I, J, K, IFC, IFACE, IFC1, JFC1, ICC1, JCC1, IC, LOHI, ILH, X1AXIS, &
           II, JJ, KK, IIF, JJF, KKF, IRC, N
REAL(EB) :: MASS_C, MASS_N(1:MAX_SURR_CELLS), RHO_CELL(1:MAX_SURR_CELLS), VOL(1:MAX_SURR_CELLS), &
            RHO_CUT, SIGN_FACTOR, SUM_MASS_N, SUM_RHO_ZZ, CONST, PRFCT, CC_RHOP, CC_RHO_ZZP, RHO_ZZ_MIN, &
            RHO_ZZ_MAX, RHO_ZZ_CUT, RHO_ZZ_TEST
LOGICAL :: CLIP_RHOMIN_CC, CLIP_RHOMAX_CC, CLIP_RHO_ZZ, CLIP_RHO_ZZ_SAVE

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
INTEGER, POINTER, DIMENSION(:,:,:,:,:) :: FCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP, RHO, RHOS
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ, ZZ, ZZS
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_RHO, DELTA_RHO_ZZ
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTFACE_TYPE), POINTER, DIMENSION(:) :: CUT_FACE
TYPE(CC_RCFACE_TYPE), POINTER, DIMENSION(:) :: RC_FACE
TYPE(CC_CUTCELL_TYPE), POINTER :: CC

! Thread-safe alias setup
IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
CELL_INDEX => M%CELL_INDEX; CCVAR => M%CCVAR; FCVAR => M%FCVAR
DX => M%DX; DY => M%DY; DZ => M%DZ
RHO => M%RHO; RHOS => M%RHOS; ZZ => M%ZZ; ZZS => M%ZZS
CELL => M%CELL
CUT_CELL => M%CUT_CELL; CUT_FACE => M%CUT_FACE; RC_FACE => M%RC_FACE

DELTA_RHO => M%WORK4
DELTA_RHO =  0._EB
CLIP_RHOMIN_CC = .FALSE.
CLIP_RHOMAX_CC = .FALSE.
IF (M%PREDICTOR) THEN
   RHOP   => RHOS
   RHO_ZZ => ZZS
   PRFCT=  1._EB
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID)  CYCLE
      CC%DELTA_RHO(1:CC%NCELL) = 0._EB
      DO JCC=1,CC%NCELL; CC%RHOS(JCC) = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC)); ENDDO
   ENDDO

ELSE
   RHOP   => RHO
   RHO_ZZ => ZZ
   PRFCT=  0._EB
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID)  CYCLE
      CC%DELTA_RHO(1:CC%NCELL) = 0._EB
      DO JCC=1,CC%NCELL; CC%RHO(JCC) = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC)); ENDDO
   ENDDO

ENDIF

! First compute RHOP in cut-cell region regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) THEN
            RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)*RHOP(I,J,K)
         ELSE
            RHOP(I,J,K) = SUM(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES))
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Correct density:
! 1. Compute DELTA_RHO in cut-cells and regular cells.
DO ICC=1,M%N_CUTCELL_MESH
   CC => CUT_CELL(ICC)
   I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IC=CELL_INDEX(I,J,K); IF (CELL(IC)%SOLID) CYCLE
   JCC1_LOOP : DO JCC=1,CC%NCELL
      CC_RHOP = PRFCT*CC%RHOS(JCC) + (1._EB-PRFCT)*CC%RHO(JCC)
      IF (CC_RHOP>=RHOMIN .AND. CC_RHOP<=RHOMAX) CYCLE JCC1_LOOP
      IF (CC_RHOP<RHOMIN) THEN
         RHO_CUT = RHOMIN
         SIGN_FACTOR = 1._EB
         CLIP_RHOMIN_CC = .TRUE.
      ELSE
         RHO_CUT = RHOMAX
         SIGN_FACTOR = -1._EB
         CLIP_RHOMAX_CC = .TRUE.
      ENDIF
      MASS_C = ABS(RHO_CUT-CC_RHOP) * CC%VOLUME(JCC)

      MASS_N = 0._EB; NCELL=0
      DO IFC=2,CC%CCELEM(1,JCC)+1
         IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
         LOHI  = CC%FACE_LIST(2,IFACE)
         ILH   = 2*CC%FACE_LIST(2,IFACE) - 3
         X1AXIS= CC%FACE_LIST(3,IFACE)
         IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
         NCELL=NCELL+1
         SELECT CASE(CC%FACE_LIST(1,IFACE))
         CASE(CC_FTYPE_CFGAS)
            IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
            ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
            RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
            VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
         CASE(CC_FTYPE_RCGAS)
            II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
            SELECT CASE(X1AXIS)
            CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
            CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
            CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
            END SELECT
            IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
            SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
            CASE(CC_FTYPE_RGGAS)
               RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
            CASE(CC_FTYPE_CFGAS)
               ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
               RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
               VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
            END SELECT
         END SELECT
         MASS_N(NCELL) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHO_CELL(NCELL)))-RHO_CUT) * VOL(NCELL)
      ENDDO

      SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
      CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)

      CC%DELTA_RHO(JCC) = CC%DELTA_RHO(JCC) + CONST*SUM_MASS_N/CC%VOLUME(JCC)
      NCELL=0
      DO IFC=2,CC%CCELEM(1,JCC)+1
         IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
         LOHI  = CC%FACE_LIST(2,IFACE)
         ILH   = 2*CC%FACE_LIST(2,IFACE) - 3
         X1AXIS= CC%FACE_LIST(3,IFACE)
         IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
         NCELL=NCELL+1
         SELECT CASE(CC%FACE_LIST(1,IFACE))
         CASE(CC_FTYPE_CFGAS)
            IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
            ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
            CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
         CASE(CC_FTYPE_RCGAS)
            II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
            SELECT CASE(X1AXIS)
            CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
            CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
            CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
            END SELECT
            IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
            SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
            CASE(CC_FTYPE_RGGAS)
               DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
            CASE(CC_FTYPE_CFGAS)
               ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
               CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
            END SELECT
         END SELECT
      ENDDO
   ENDDO JCC1_LOOP
ENDDO

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IF (RHOP(I,J,K)>=RHOMIN .AND. RHOP(I,J,K)<=RHOMAX) CYCLE
         IF (RHOP(I,J,K)<RHOMIN) THEN
            RHO_CUT = RHOMIN
            SIGN_FACTOR = 1._EB
            CLIP_RHOMIN_CC = .TRUE.
         ELSE
            RHO_CUT = RHOMAX
            SIGN_FACTOR = -1._EB
            CLIP_RHOMAX_CC = .TRUE.
         ENDIF
         MASS_C = ABS(RHO_CUT-RHOP(I,J,K)) * (DX(I)*DY(J)*DZ(K))
         MASS_N = 0._EB; NCELL=0
         IC = CELL_INDEX(I,J,K)
         DO X1AXIS=IAXIS,KAXIS
            DO LOHI=LOW_IND,HIGH_IND
               ILH   = 2*LOHI - 3
               IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
               NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               IF (IRC>0) THEN
                  SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                  CASE(CC_FTYPE_RGGAS)
                     RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                  CASE(CC_FTYPE_CFGAS)
                     ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                     RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%RHOS(JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%RHO(JCC1)
                     VOL(NCELL) = CUT_CELL(ICC1)%VOLUME(JCC1)
                  END SELECT
               ELSE
                  RHO_CELL(NCELL) = RHOP(II,JJ,KK); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
               ENDIF
               MASS_N(NCELL) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHO_CELL(NCELL)))-RHO_CUT) * VOL(NCELL)
            ENDDO
         ENDDO
         SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
         DELTA_RHO(I,J,K) = DELTA_RHO(I,J,K) + CONST*SUM_MASS_N/(DX(I)*DY(J)*DZ(K))
         NCELL=0
         DO X1AXIS=IAXIS,KAXIS
            DO LOHI=LOW_IND,HIGH_IND
               ILH   = 2*LOHI - 3
               IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
               NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               IF (IRC>0) THEN
                  SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                  CASE(CC_FTYPE_RGGAS)
                     DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  CASE(CC_FTYPE_CFGAS)
                     ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                     CUT_CELL(ICC1)%DELTA_RHO(JCC1) = CUT_CELL(ICC1)%DELTA_RHO(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  END SELECT
               ELSE
                  DELTA_RHO(II,JJ,KK) = DELTA_RHO(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

! 2. Assign DELTA_RHO to neighboring cells if clipping has been done.
IF (CLIP_RHOMIN_CC .OR. CLIP_RHOMAX_CC) THEN
   IF (M%PREDICTOR) THEN
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%RHOS(JCC) = MIN(RHOMAX,MAX(RHOMIN,CC%RHOS(JCC)+CC%DELTA_RHO(JCC))); ENDDO
      ENDDO
   ELSE
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%RHO(JCC)  = MIN(RHOMAX,MAX(RHOMIN,CC%RHO(JCC) +CC%DELTA_RHO(JCC))); ENDDO
      ENDDO
   ENDIF
   RHOP(1:IBAR,1:JBAR,1:KBAR) = MIN(RHOMAX,MAX(RHOMIN,RHOP(1:IBAR,1:JBAR,1:KBAR)+DELTA_RHO(1:IBAR,1:JBAR,1:KBAR)))
ENDIF

! Thread-safe per-mesh clip flag update:
M%CLIP_RHOMIN = M%CLIP_RHOMIN .OR. CLIP_RHOMIN_CC
M%CLIP_RHOMAX = M%CLIP_RHOMAX .OR. CLIP_RHOMAX_CC

IF (N_TRACKED_SPECIES==1) THEN
   IF (CLIP_RHOMIN_CC .OR. CLIP_RHOMAX_CC) THEN
      IF (M%PREDICTOR) THEN
         DO ICC=1,M%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%ZZS(1,JCC) = CC%RHOS(JCC); ENDDO
         ENDDO
      ELSE
         DO ICC=1,M%N_CUTCELL_MESH
            CC => CUT_CELL(ICC)
            IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL; CC%ZZ(1,JCC)  = CC%RHO(JCC); ENDDO
         ENDDO
      ENDIF
      RHO_ZZ(1:IBAR,1:JBAR,1:KBAR,1) = RHOP(1:IBAR,1:JBAR,1:KBAR)
   ENDIF
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) &
               RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   RETURN
ENDIF

! Correct species mass density
RHO_ZZ_MIN       =  0._EB
DELTA_RHO_ZZ     => M%WORK5
CLIP_RHO_ZZ_SAVE = .FALSE.
SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   DELTA_RHO_ZZ =  0._EB
   DO ICC=1,M%N_CUTCELL_MESH
      CUT_CELL(ICC)%DELTA_RHO_ZZ=0._EB
   ENDDO
   CLIP_RHO_ZZ  = .FALSE.

   ! 1. compute DELTA_RHO_ZZ in cut-cells and regular cells.
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      I=CC%IJK(IAXIS); J=CC%IJK(JAXIS); K=CC%IJK(KAXIS); IC=CELL_INDEX(I,J,K); IF (CELL(IC)%SOLID) CYCLE
      JCC2_LOOP : DO JCC=1,CC%NCELL
         CC_RHOP    = PRFCT*CC%RHOS(JCC)  + (1._EB-PRFCT)*CC%RHO(JCC); RHO_ZZ_MAX = CC_RHOP
         CC_RHO_ZZP = PRFCT*CC%ZZS(N,JCC) + (1._EB-PRFCT)*CC%ZZ(N,JCC)
         IF (CC_RHO_ZZP>=RHO_ZZ_MIN .AND. CC_RHO_ZZP<=RHO_ZZ_MAX) CYCLE JCC2_LOOP
         CLIP_RHO_ZZ = .TRUE.
         IF (CC_RHO_ZZP<RHO_ZZ_MIN) THEN
            RHO_ZZ_CUT = RHO_ZZ_MIN
            SIGN_FACTOR = 1._EB
         ELSE
            RHO_ZZ_CUT = RHO_ZZ_MAX
            SIGN_FACTOR = -1._EB
         ENDIF
         MASS_C = ABS(RHO_ZZ_CUT-CC_RHO_ZZP) * CC%VOLUME(JCC)

         MASS_N = 0._EB; NCELL=0
         DO IFC=2,CC%CCELEM(1,JCC)+1
            IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
            LOHI  = CC%FACE_LIST(2,IFACE)
            ILH   = 2*CC%FACE_LIST(2,IFACE) - 3
            X1AXIS= CC%FACE_LIST(3,IFACE)
            IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
            NCELL=NCELL+1
            SELECT CASE(CC%FACE_LIST(1,IFACE))
            CASE(CC_FTYPE_CFGAS)
               IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
               ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
               RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
               VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
            CASE(CC_FTYPE_RCGAS)
               II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
               CASE(CC_FTYPE_RGGAS)
                  RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
               CASE(CC_FTYPE_CFGAS)
                  ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                  RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
                  VOL(NCELL)      = CUT_CELL(ICC1)%VOLUME(JCC1)
               END SELECT
            END SELECT
            MASS_N(NCELL) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_CELL(NCELL)))-RHO_ZZ_CUT) * VOL(NCELL)
         ENDDO

         SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)

         CC%DELTA_RHO_ZZ(JCC) = CC%DELTA_RHO_ZZ(JCC) + CONST*SUM_MASS_N/CC%VOLUME(JCC)
         NCELL=0
         DO IFC=2,CC%CCELEM(1,JCC)+1
            IFACE = CC%CCELEM(IFC,JCC); IF(CC%FACE_LIST(1,IFACE)==CC_FTYPE_CFINB) CYCLE
            LOHI  = CC%FACE_LIST(2,IFACE)
            ILH   = 2*CC%FACE_LIST(2,IFACE) - 3
            X1AXIS= CC%FACE_LIST(3,IFACE)
            IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
            NCELL=NCELL+1
            SELECT CASE(CC%FACE_LIST(1,IFACE))
            CASE(CC_FTYPE_CFGAS)
               IFC1 = CC%FACE_LIST(4,IFACE); JFC1 = CC%FACE_LIST(5,IFACE)
               ICC1 = CUT_FACE(IFC1)%CELL_LIST(2,LOHI,JFC1); JCC1 = CUT_FACE(IFC1)%CELL_LIST(3,LOHI,JFC1)
               CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
            CASE(CC_FTYPE_RCGAS)
               II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
               SELECT CASE(X1AXIS)
               CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
               CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
               CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
               END SELECT
               IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
               SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
               CASE(CC_FTYPE_RGGAS)
                  DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
               CASE(CC_FTYPE_CFGAS)
                  ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                  CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) - CONST*MASS_N(NCELL)/VOL(NCELL)
               END SELECT
            END SELECT
         ENDDO
      ENDDO JCC2_LOOP
   ENDDO
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            RHO_ZZ_MAX = RHOP(I,J,K)
            IF (RHO_ZZ(I,J,K,N)>=RHO_ZZ_MIN .AND. RHO_ZZ(I,J,K,N)<=RHO_ZZ_MAX) CYCLE
            CLIP_RHO_ZZ = .TRUE.
            IF (RHO_ZZ(I,J,K,N)<RHO_ZZ_MIN) THEN
               RHO_ZZ_CUT = RHO_ZZ_MIN
               SIGN_FACTOR = 1._EB
            ELSE
               RHO_ZZ_CUT = RHO_ZZ_MAX
               SIGN_FACTOR = -1._EB
            ENDIF
            MASS_C = ABS(RHO_ZZ_CUT-RHO_ZZ(I,J,K,N)) * (DX(I)*DY(J)*DZ(K))
            MASS_N = 0._EB; NCELL=0
            IC = CELL_INDEX(I,J,K)
            DO X1AXIS=IAXIS,KAXIS
               DO LOHI=LOW_IND,HIGH_IND
                  ILH   = 2*LOHI - 3
                  IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
                  NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  IF (IRC>0) THEN
                     SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                     CASE(CC_FTYPE_RGGAS)
                        RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                     CASE(CC_FTYPE_CFGAS)
                        ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                        RHO_CELL(NCELL) = PRFCT*CUT_CELL(ICC1)%ZZS(N,JCC1)+(1._EB-PRFCT)*CUT_CELL(ICC1)%ZZ(N,JCC1)
                        VOL(NCELL) = CUT_CELL(ICC1)%VOLUME(JCC1)
                     END SELECT
                  ELSE
                     RHO_CELL(NCELL) = RHO_ZZ(II,JJ,KK,N); VOL(NCELL) = DX(II)*DY(JJ)*DZ(KK)
                  ENDIF
                  MASS_N(NCELL) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_CELL(NCELL)))-RHO_ZZ_CUT) * VOL(NCELL)
               ENDDO
            ENDDO
            SUM_MASS_N = SUM(MASS_N(1:NCELL)); IF (SUM_MASS_N<=TWENTY_EPSILON_EB) CYCLE
            CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
            DELTA_RHO_ZZ(I,J,K) = DELTA_RHO_ZZ(I,J,K) + CONST*SUM_MASS_N/(DX(I)*DY(J)*DZ(K))
            NCELL=0
            DO X1AXIS=IAXIS,KAXIS
               DO LOHI=LOW_IND,HIGH_IND
                  ILH   = 2*LOHI - 3
                  IF (CELL(IC)%WALL_INDEX(ILH*X1AXIS)/=0) CYCLE
                  NCELL=NCELL+1; II=I; JJ=J; KK=K; IIF=I; JJF=J; KKF=K
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS); IIF = IIF+LOHI-2; II = II+ILH
                  CASE(JAXIS); JJF = JJF+LOHI-2; JJ = JJ+ILH
                  CASE(KAXIS); KKF = KKF+LOHI-2; KK = KK+ILH
                  END SELECT
                  IRC = FCVAR(IIF,JJF,KKF,CC_IDRC,X1AXIS)
                  IF (IRC>0) THEN
                     SELECT CASE(RC_FACE(IRC)%CELL_LIST(1,LOHI))
                     CASE(CC_FTYPE_RGGAS)
                        DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                     CASE(CC_FTYPE_CFGAS)
                        ICC1 = RC_FACE(IRC)%CELL_LIST(2,LOHI);  JCC1 = RC_FACE(IRC)%CELL_LIST(3,LOHI)
                        CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1) = CUT_CELL(ICC1)%DELTA_RHO_ZZ(JCC1)-CONST*MASS_N(NCELL)/VOL(NCELL)
                     END SELECT
                  ELSE
                     DELTA_RHO_ZZ(II,JJ,KK) = DELTA_RHO_ZZ(II,JJ,KK) - CONST*MASS_N(NCELL)/VOL(NCELL)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   IF (.NOT.CLIP_RHO_ZZ) THEN
      CYCLE SPECIES_LOOP
   ELSE
      CLIP_RHO_ZZ_SAVE = .TRUE.
   ENDIF

   ! 2. Assign excess/deficit RHO_ZZ neighboring cells
   IF (M%PREDICTOR) THEN
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZS(N,JCC) = MIN(CC%RHOS(JCC),MAX(RHO_ZZ_MIN,CC%ZZS(N,JCC)+CC%DELTA_RHO_ZZ(JCC))); ENDDO
      ENDDO
   ELSE
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC)
         IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZ(N,JCC)  = MIN(CC%RHO(JCC), MAX(RHO_ZZ_MIN,CC%ZZ(N,JCC) +CC%DELTA_RHO_ZZ(JCC))); ENDDO
      ENDDO
   ENDIF
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            RHO_ZZ(I,J,K,N) = MIN(RHOP(I,J,K),MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K,N)+DELTA_RHO_ZZ(I,J,K)))
         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

! If nothing has been clipped, return
IF (.NOT.CLIP_RHOMIN_CC .AND. .NOT.CLIP_RHOMAX_CC .AND. .NOT.CLIP_RHO_ZZ_SAVE) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE)  CYCLE
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) &
               RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   RETURN
ENDIF

! Final check of RHO_ZZ => SUM(RHO_ZZ) = RHO
IF (M%PREDICTOR) THEN
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         SUM_RHO_ZZ = SUM(CC%ZZS(1:N_TRACKED_SPECIES,JCC))
         N = MAXLOC(CC%ZZS(1:N_TRACKED_SPECIES,JCC),1)
         RHO_ZZ_TEST = CC%ZZS(N,JCC) + CC%RHOS(JCC) - SUM_RHO_ZZ
         IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>CC%RHOS(JCC)) THEN
            CC%ZZS(1:N_TRACKED_SPECIES,JCC) = CC%RHOS(JCC) * CC%ZZS(1:N_TRACKED_SPECIES,JCC)/SUM_RHO_ZZ
         ELSE
            CC%ZZS(N,JCC) = RHO_ZZ_TEST
         ENDIF
      ENDDO
   ENDDO
ELSE
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC)
      IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         SUM_RHO_ZZ = SUM(CC%ZZ(1:N_TRACKED_SPECIES,JCC))
         N = MAXLOC(CC%ZZ(1:N_TRACKED_SPECIES,JCC),1)
         RHO_ZZ_TEST = CC%ZZ(N,JCC) + CC%RHO(JCC) - SUM_RHO_ZZ
         IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>CC%RHO(JCC)) THEN
            CC%ZZ(1:N_TRACKED_SPECIES,JCC) = CC%RHO(JCC) * CC%ZZ(1:N_TRACKED_SPECIES,JCC)/SUM_RHO_ZZ
         ELSE
            CC%ZZ(N,JCC) = RHO_ZZ_TEST
         ENDIF
      ENDDO
   ENDDO
ENDIF
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE
         SUM_RHO_ZZ = SUM(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES))
         N = MAXLOC(RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES),1)
         RHO_ZZ_TEST = RHO_ZZ(I,J,K,N) + RHOP(I,J,K) - SUM_RHO_ZZ
         IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>RHOP(I,J,K)) THEN
            RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHOP(I,J,K) * RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/SUM_RHO_ZZ
         ELSE
            RHO_ZZ(I,J,K,N) = RHO_ZZ_TEST
         ENDIF
      ENDDO
   ENDDO
ENDDO

CALL CC_CV_RHOZZ_AVERAGE_TS

! Bring back ZZ in regular gas cells from partial densities.
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) &
            RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES) = RHO_ZZ(I,J,K,1:N_TRACKED_SPECIES)/RHOP(I,J,K)
      ENDDO
   ENDDO
ENDDO

RETURN
CONTAINS

SUBROUTINE CC_CV_RHOZZ_AVERAGE_TS

INTEGER :: IROW_LOC
INTEGER :: ILC_LO, ILC_HI

ILC_LO = UNKZ_ILC(NM) + 1
ILC_HI = UNKZ_ILC(NM) + NUNKZ_LOC(NM)

! CV volumes:
RZ_ZS(ILC_LO:ILC_HI) = 0._EB
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + DX(I)*DY(J)*DZ(K)
      ENDDO
   ENDDO
ENDDO
DO ICC=1,M%N_CUTCELL_MESH
   CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
   DO JCC=1,CC%NCELL
      IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
      RZ_ZS(IROW_LOC) = RZ_ZS(IROW_LOC) + CC%VOLUME(JCC)
   ENDDO
ENDDO

! Loop species:
SPECIES_LOOP_AVG: DO N=1,N_TRACKED_SPECIES

   RZ_Z(ILC_LO:ILC_HI) = 0._EB
   ! Add to CV rhoZZ*Vol:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            RZ_Z( IROW_LOC) = RZ_Z( IROW_LOC) + RHO_ZZ(I,J,K,N)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   IF (M%PREDICTOR) THEN
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%ZZS(N,JCC) * CC%VOLUME(JCC)
         ENDDO
      ENDDO
   ELSE
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW_LOC = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) + CC%ZZ(N,JCC) * CC%VOLUME(JCC)
         ENDDO
      ENDDO
   ENDIF

   ! Volume average (per-mesh range):
   DO IROW_LOC=UNKZ_IND(NM)-UNKZ_IND(NM_START)+1,UNKZ_IND(NM)-UNKZ_IND(NM_START)+NUNKZ_LOC(NM)
      RZ_Z(IROW_LOC)  = RZ_Z( IROW_LOC) / RZ_ZS(IROW_LOC)
   ENDDO

   ! Back to cut/reg cell containers of rhoZZ:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            RHO_ZZ(I,J,K,N) = RZ_Z(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
         ENDDO
      ENDDO
   ENDDO
   IF (M%PREDICTOR) THEN
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZS(N,JCC) = RZ_Z(CC%UNKZ(JCC)-UNKZ_IND(NM_START)); ENDDO
      ENDDO
   ELSE
      DO ICC=1,M%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
         DO JCC=1,CC%NCELL; CC%ZZ(N,JCC) = RZ_Z(CC%UNKZ(JCC)-UNKZ_IND(NM_START)); ENDDO
      ENDDO
   ENDIF

ENDDO SPECIES_LOOP_AVG

RZ_ZS(ILC_LO:ILC_HI) = 0._EB

END SUBROUTINE CC_CV_RHOZZ_AVERAGE_TS

END SUBROUTINE CC_CHECK_MASS_DENSITY_TS


! ----------------------------- CC_DENSITY_EXPLICIT_TS ------------------------

RECURSIVE SUBROUTINE CC_DENSITY_EXPLICIT_TS(NM,M,T,DT)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: T,DT

! Local variables:
INTEGER :: N
INTEGER :: IROW_LOC
REAL(EB):: DUMMYT
INTEGER :: ILC_LO, ILC_HI

DUMMYT = T

! Per-mesh range in global arrays
ILC_LO = UNKZ_ILC(NM) + 1
ILC_HI = UNKZ_ILC(NM) + NUNKZ_LOC(NM)

SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

   IF( (M%PREDICTOR.AND.FIRST_PASS) .OR. M%CORRECTOR) THEN
      ! Zero per-mesh range of F_Z:
      F_Z(ILC_LO:ILC_HI) = 0._EB
      CALL GET_EXPLICIT_ADVDIFFVECTOR_SCALAR_3D_TS(NM,M,N)

      CALL GET_ADVDIFFVECTOR_SCALAR_3D_TS(NM,M,N)

      CALL GET_M_DOT_PPP_SCALAR_3D_TS(NM,M,N)

      IF (PERIODIC_TEST==7) CALL GET_SHUNN3_QZ_TS(NM,M,T,N)

      CALL GET_RHOZZVECTOR_SCALAR_3D_TS(NM,M,N)
   ENDIF

   IF (M%PREDICTOR) THEN
      IF (FIRST_PASS) THEN
         F_Z0(ILC_LO:ILC_HI,N)  = F_Z(ILC_LO:ILC_HI)
         RZ_Z0(ILC_LO:ILC_HI,N) = RZ_Z(ILC_LO:ILC_HI)
      ELSE
         F_Z(ILC_LO:ILC_HI)  = F_Z0(ILC_LO:ILC_HI,N)
         RZ_Z(ILC_LO:ILC_HI) = RZ_Z0(ILC_LO:ILC_HI,N)
      ENDIF
   ENDIF

   IF (M%PREDICTOR) THEN
      DO IROW_LOC=ILC_LO,ILC_HI
         RZ_Z(IROW_LOC) = RZ_Z(IROW_LOC) - DT * F_Z(IROW_LOC) / M_MAT_Z(IROW_LOC)
      ENDDO
   ELSE
      DO IROW_LOC=ILC_LO,ILC_HI
         RZ_Z(IROW_LOC) = RZ_Z(IROW_LOC) - 0.5_EB * DT * F_Z(IROW_LOC) / M_MAT_Z(IROW_LOC)
      ENDDO
   ENDIF

   CALL PUT_RHOZZVECTOR_SCALAR_3D_TS(NM,M,N)

ENDDO SPECIES_LOOP

CALL CC_CHECK_MASS_DENSITY_TS(NM,M)

CALL GET_RHOZZ_CC_3D_TS(NM,M)

RETURN
END SUBROUTINE CC_DENSITY_EXPLICIT_TS


! ------------------------------ CC_DENSITY_TS -------------------------------

RECURSIVE SUBROUTINE CC_DENSITY_TS(NM,M,T,DT)

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: T,DT

! Local Variables:
INTEGER :: N
INTEGER :: I,J,K,ICC,JCC
REAL(EB) :: ZZ_GET(1:N_TRACKED_SPECIES),VCCELL,PBAR_K

! Thread-safe local shadows
INTEGER, POINTER :: IBAR, JBAR, KBAR, IBP1, JBP1, KBP1
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX, PRESSURE_ZONE
INTEGER, POINTER, DIMENSION(:,:,:,:) :: CCVAR
REAL(EB), POINTER, DIMENSION(:) :: DX, DY, DZ
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR, PBAR_S
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO, RHOS, TMP, RSUM
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ, ZZS
TYPE(CELL_TYPE), POINTER, DIMENSION(:) :: CELL
TYPE(CC_CUTCELL_TYPE), POINTER, DIMENSION(:) :: CUT_CELL
TYPE(CC_CUTCELL_TYPE), POINTER :: CC

IF (SOLID_PHASE_ONLY) RETURN

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
      ! CONTINUE
END SELECT

! Thread-safe alias setup
IBAR => M%IBAR; JBAR => M%JBAR; KBAR => M%KBAR
IBP1 => M%IBP1; JBP1 => M%JBP1; KBP1 => M%KBP1
CELL_INDEX => M%CELL_INDEX; CCVAR => M%CCVAR; PRESSURE_ZONE => M%PRESSURE_ZONE
DX => M%DX; DY => M%DY; DZ => M%DZ
PBAR => M%PBAR; PBAR_S => M%PBAR_S
RHO => M%RHO; RHOS => M%RHOS; TMP => M%TMP; RSUM => M%RSUM
ZZ => M%ZZ; ZZS => M%ZZS
CELL => M%CELL
CUT_CELL => M%CUT_CELL

! Advance scalars and density:
CALL CC_DENSITY_EXPLICIT_TS(NM,M,T,DT)

! Compute molecular weight term and temperature from EOS:
IF (M%PREDICTOR) THEN

   ! First Regular Cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            PBAR_K = PBAR_S(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
            TMP(I,J,K) = PBAR_K/(RSUM(I,J,K)*RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Store RHO*ZZ values at step n:
   IF (.NOT.ALLOCATED(M%RHO_ZZN)) ALLOCATE(M%RHO_ZZN(0:IBP1,0:JBP1,0:KBP1,N_TOTAL_SCALARS))
   DO N=1,N_TOTAL_SCALARS
      M%RHO_ZZN(0:IBP1,0:JBP1,0:KBP1,N) = RHO(0:IBP1,0:JBP1,0:KBP1)*ZZ(0:IBP1,0:JBP1,0:KBP1,N)
   ENDDO

   ! Cut-cells:
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      I  = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      VCCELL = 0._EB; TMP(I,J,K)=0._EB; RHOS(I,J,K)=0._EB; ZZS(I,J,K,1:N_TRACKED_SPECIES)=0._EB; RSUM(I,J,K)=0._EB
      DO JCC=1,CC%NCELL
         ZZ_GET(1:N_TRACKED_SPECIES) = CC%ZZS(1:N_TRACKED_SPECIES,JCC)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,CC%RSUM(JCC))
         PBAR_K = PBAR_S(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CC%UNKZ(JCC)-UNKZ_IND(NM_START))
         CC%TMP(JCC) = PBAR_K/(CC%RSUM(JCC)*CC%RHOS(JCC))
         TMP(I,J,K) = TMP(I,J,K) + CC%TMP(JCC)*CC%VOLUME(JCC)
         RHOS(I,J,K)= RHOS(I,J,K)+ CC%RHOS(JCC)*CC%VOLUME(JCC)
         ZZS(I,J,K,1:N_TRACKED_SPECIES) = ZZS(I,J,K,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*CC%VOLUME(JCC)
         RSUM(I,J,K)= RSUM(I,J,K)+ CC%RSUM(JCC)*CC%VOLUME(JCC)
         VCCELL = VCCELL + CC%VOLUME(JCC)
      ENDDO
      TMP(I,J,K) = TMP(I,J,K)/VCCELL
      RHOS(I,J,K)= RHOS(I,J,K)/VCCELL
      ZZS(I,J,K,1:N_TRACKED_SPECIES)=ZZS(I,J,K,1:N_TRACKED_SPECIES)/VCCELL
      RSUM(I,J,K)=RSUM(I,J,K)/VCCELL
   ENDDO

   ! Set to ambient temperature the temp of SOLID cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_CGSC) /= CC_SOLID) CYCLE
            TMP(I,J,K) = TMPA
         ENDDO
      ENDDO
   ENDDO

ELSE ! CORRECTOR

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUM(I,J,K))
            PBAR_K = PBAR(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START))
            TMP(I,J,K) = PBAR_K/(RSUM(I,J,K)*RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Cut-cells:
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      I  = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      VCCELL = 0._EB; TMP(I,J,K)=0._EB; RHO(I,J,K)=0._EB; ZZ(I,J,K,1:N_TRACKED_SPECIES)=0._EB; RSUM(I,J,K)=0._EB
      DO JCC=1,CC%NCELL
         ZZ_GET(1:N_TRACKED_SPECIES) = CC%ZZ(1:N_TRACKED_SPECIES,JCC)
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,CC%RSUM(JCC))
         PBAR_K = PBAR(K,PRESSURE_ZONE(I,J,K)) - P_0(K) + P_0_CV(CC%UNKZ(JCC)-UNKZ_IND(NM_START))
         CC%TMP(JCC) = PBAR_K/(CC%RSUM(JCC)*CC%RHO(JCC))
         TMP(I,J,K) = TMP(I,J,K) + CC%TMP(JCC)*CC%VOLUME(JCC)
         RHO(I,J,K) = RHO(I,J,K) + CC%RHO(JCC)*CC%VOLUME(JCC)
         ZZ(I,J,K,1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES) + ZZ_GET(1:N_TRACKED_SPECIES)*CC%VOLUME(JCC)
         RSUM(I,J,K)= RSUM(I,J,K)+ CC%RSUM(JCC)*CC%VOLUME(JCC)
         VCCELL = VCCELL + CC%VOLUME(JCC)
      ENDDO
      TMP(I,J,K) = TMP(I,J,K)/VCCELL
      RHO(I,J,K) = RHO(I,J,K)/VCCELL
      ZZ(I,J,K,1:N_TRACKED_SPECIES)=ZZ(I,J,K,1:N_TRACKED_SPECIES)/VCCELL
      RSUM(I,J,K)=RSUM(I,J,K)/VCCELL
   ENDDO

   ! Set to ambient temperature the temp of SOLID cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_CGSC) /= CC_SOLID) CYCLE
            TMP(I,J,K) = TMPA
         ENDDO
      ENDDO
   ENDDO

ENDIF ! PREDICTOR

RETURN
END SUBROUTINE CC_DENSITY_TS


END MODULE CC_DENSITY_MOD
