!> \brief Pure computation kernels extracted from MASS module.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE MASS_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE,EXTERNAL_WALL_TYPE,SPECIES_MIXTURE_TYPE,SPECIES_MIXTURE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC MASS_FINITE_DIFFERENCES_NEW_KERNEL,DENSITY_KERNEL, &
       DENSITY_BLOCK_PREPROCESSING,DENSITY_BLOCK_KERNEL_COMPUTE,DENSITY_BLOCK_POSTPROCESSING

CONTAINS


!> \brief Compute spatial differences for mass transport equations
!> \param M Mesh data structure

SUBROUTINE MASS_FINITE_DIFFERENCES_NEW_KERNEL(M)

USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE
USE PHYSICAL_FUNCTIONS, ONLY: GET_MOLECULAR_WEIGHT

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB) :: MW_F,MW_G,ZZ_GET(1:N_TRACKED_SPECIES)
REAL(EB), TARGET, DIMENSION(0:3,0:3,0:3) :: F_WORK
REAL(EB), TARGET, DIMENSION(-1:3,-1:3,-1:3) :: U_WORK,Z_WORK
REAL(EB), POINTER, DIMENSION(:,:,:) :: F_TEMP,U_TEMP,Z_TEMP,FX_P,FY_P,FZ_P
REAL(EB), PARAMETER :: DUMMY=0._EB
INTEGER  :: I,J,K,N,IOR,IW,IIG,JJG,KKG,II,JJ,KK,IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHO_Z_P,RHO_RMW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

IF (SOLID_PHASE_ONLY) RETURN

IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
   RHOP => M%RHO
   ZZP => M%ZZ
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
   RHOP => M%RHOS
   ZZP => M%ZZS
ENDIF

! Reset counter for CLIP_RHOMIN, CLIP_RHOMAX
! Done here so DT_RESTRICT_COUNT will persist until WRITE_DIAGNOSTICS is called

IF (PREDICTOR) M%DT_RESTRICT_COUNT = 0

! Species face values

SPECIES_LOOP: DO N=1,N_TOTAL_SCALARS

   RHO_Z_P=>M%WORK_PAD

   DO K=-1,M%KBP1+1
      DO J=-1,M%JBP1+1
         DO I=-1,M%IBP1+1
            RHO_Z_P(I,J,K) = RHOP(I,J,K)*ZZP(I,J,K,N)
         ENDDO
      ENDDO
   ENDDO

   ! Compute scalar face values

   FX_P(LBOUND(M%FX,1):,LBOUND(M%FX,2):,LBOUND(M%FX,3):) => M%FX(:,:,:,N)
   FY_P(LBOUND(M%FY,1):,LBOUND(M%FY,2):,LBOUND(M%FY,3):) => M%FY(:,:,:,N)
   FZ_P(LBOUND(M%FZ,1):,LBOUND(M%FZ,2):,LBOUND(M%FZ,3):) => M%FZ(:,:,:,N)
   CALL GET_SCALAR_FACE_VALUE(UU,RHO_Z_P,FX_P,0,M%IBAR,1,M%JBAR,1,M%KBAR,1,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(VV,RHO_Z_P,FY_P,1,M%IBAR,0,M%JBAR,1,M%KBAR,2,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(WW,RHO_Z_P,FZ_P,1,M%IBAR,1,M%JBAR,0,M%KBAR,3,I_FLUX_LIMITER)

   U_TEMP => U_WORK
   F_TEMP => F_WORK
   Z_TEMP => Z_WORK
   WALL_LOOP_2: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_2
      BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
      B1=>M%BOUNDARY_PROP1(WC%B1_INDEX)

      II  = BC%II
      JJ  = BC%JJ
      KK  = BC%KK
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      IOR = BC%IOR
      IC  = M%CELL_INDEX(II,JJ,KK)

      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. .NOT.M%CELL(IC)%SOLID .AND. .NOT.M%CELL(IC)%EXTERIOR) THEN
         ! thin obstruction
         SELECT CASE(IOR)
            CASE( 1); M%FX(IIG-1,JJG,KKG,N) = 0._EB
            CASE(-1); M%FX(IIG,JJG,KKG,N)   = 0._EB
            CASE( 2); M%FY(IIG,JJG-1,KKG,N) = 0._EB
            CASE(-2); M%FY(IIG,JJG,KKG,N)   = 0._EB
            CASE( 3); M%FZ(IIG,JJG,KKG-1,N) = 0._EB
            CASE(-3); M%FZ(IIG,JJG,KKG,N)   = 0._EB
         END SELECT
      ELSE
         IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN
            SELECT CASE(IOR)
               CASE( 1); M%FX(IIG-1,JJG,KKG,N) = B1%RHO_F*B1%ZZ_F(N)
               CASE(-1); M%FX(IIG,JJG,KKG,N)   = B1%RHO_F*B1%ZZ_F(N)
               CASE( 2); M%FY(IIG,JJG-1,KKG,N) = B1%RHO_F*B1%ZZ_F(N)
               CASE(-2); M%FY(IIG,JJG,KKG,N)   = B1%RHO_F*B1%ZZ_F(N)
               CASE( 3); M%FZ(IIG,JJG,KKG-1,N) = B1%RHO_F*B1%ZZ_F(N)
               CASE(-3); M%FZ(IIG,JJG,KKG,N)   = B1%RHO_F*B1%ZZ_F(N)
            END SELECT
         ENDIF
      ENDIF

      ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

      OFF_WALL_IF_2: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

         OFF_WALL_SELECT_2: SELECT CASE(IOR)
            CASE( 1) OFF_WALL_SELECT_2
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II+1,JJ,KK))%WALL_INDEX(+1)>0)) THEN
                  Z_TEMP(0:3,1,1) = (/RHO_Z_P(II+1,JJ,KK),RHO_Z_P(II+1:II+2,JJ,KK),DUMMY/)
                  U_TEMP(1,1,1) = UU(II+1,JJ,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  M%FX(II+1,JJ,KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_2
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II-1,JJ,KK))%WALL_INDEX(-1)>0)) THEN
                  Z_TEMP(0:3,1,1) = (/DUMMY,RHO_Z_P(II-2:II-1,JJ,KK),RHO_Z_P(II-1,JJ,KK)/)
                  U_TEMP(1,1,1) = UU(II-2,JJ,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  M%FX(II-2,JJ,KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ+1,KK))%WALL_INDEX(+2)>0)) THEN
                  Z_TEMP(1,0:3,1) = (/RHO_Z_P(II,JJ+1,KK),RHO_Z_P(II,JJ+1:JJ+2,KK),DUMMY/)
                  U_TEMP(1,1,1) = VV(II,JJ+1,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  M%FY(II,JJ+1,KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_2
               IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ-1,KK))%WALL_INDEX(-2)>0)) THEN
                  Z_TEMP(1,0:3,1) = (/DUMMY,RHO_Z_P(II,JJ-2:JJ-1,KK),RHO_Z_P(II,JJ-1,KK)/)
                  U_TEMP(1,1,1) = VV(II,JJ-2,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  M%FY(II,JJ-2,KK,N) = F_TEMP(1,1,1)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ,KK+1))%WALL_INDEX(+3)>0)) THEN
                  Z_TEMP(1,1,0:3) = (/RHO_Z_P(II,JJ,KK+1),RHO_Z_P(II,JJ,KK+1:KK+2),DUMMY/)
                  U_TEMP(1,1,1) = WW(II,JJ,KK+1)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  M%FZ(II,JJ,KK+1,N) = F_TEMP(1,1,1)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_2
               IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ,KK-1))%WALL_INDEX(-3)>0)) THEN
                  Z_TEMP(1,1,0:3) = (/DUMMY,RHO_Z_P(II,JJ,KK-2:KK-1),RHO_Z_P(II,JJ,KK-1)/)
                  U_TEMP(1,1,1) = WW(II,JJ,KK-2)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  M%FZ(II,JJ,KK-2,N) = F_TEMP(1,1,1)
               ENDIF
         END SELECT OFF_WALL_SELECT_2

      ENDIF OFF_WALL_IF_2

   ENDDO WALL_LOOP_2

ENDDO SPECIES_LOOP

FACE_CORRECTION_IF: IF (FLUX_LIMITER_MW_CORRECTION) THEN

   ! Repeat the above for DENSITY

   RHO_RMW=>M%WORK_PAD

   DO K=-1,M%KBP1+1
      DO J=-1,M%JBP1+1
         DO I=-1,M%IBP1+1
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW_G)
            RHO_RMW(I,J,K) = RHOP(I,J,K)/MW_G
         ENDDO
      ENDDO
   ENDDO

   FX_P(LBOUND(M%FX,1):,LBOUND(M%FX,2):,LBOUND(M%FX,3):) => M%FX(:,:,:,0)
   FY_P(LBOUND(M%FY,1):,LBOUND(M%FY,2):,LBOUND(M%FY,3):) => M%FY(:,:,:,0)
   FZ_P(LBOUND(M%FZ,1):,LBOUND(M%FZ,2):,LBOUND(M%FZ,3):) => M%FZ(:,:,:,0)
   CALL GET_SCALAR_FACE_VALUE(UU,RHO_RMW,FX_P,0,M%IBAR,1,M%JBAR,1,M%KBAR,1,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(VV,RHO_RMW,FY_P,1,M%IBAR,0,M%JBAR,1,M%KBAR,2,I_FLUX_LIMITER)
   CALL GET_SCALAR_FACE_VALUE(WW,RHO_RMW,FZ_P,1,M%IBAR,1,M%JBAR,0,M%KBAR,3,I_FLUX_LIMITER)

   U_TEMP => U_WORK
   F_TEMP => F_WORK
   Z_TEMP => Z_WORK
   WALL_LOOP_3: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP_3
      BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
      B1=>M%BOUNDARY_PROP1(WC%B1_INDEX)

      II  = BC%II
      JJ  = BC%JJ
      KK  = BC%KK
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      IOR = BC%IOR
      IC  = M%CELL_INDEX(II,JJ,KK)

      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. .NOT.M%CELL(IC)%SOLID .AND. .NOT.M%CELL(IC)%EXTERIOR) THEN
         SELECT CASE(IOR)
            CASE( 1); M%FX(IIG-1,JJG,KKG,0) = 0._EB
            CASE(-1); M%FX(IIG,JJG,KKG,0)   = 0._EB
            CASE( 2); M%FY(IIG,JJG-1,KKG,0) = 0._EB
            CASE(-2); M%FY(IIG,JJG,KKG,0)   = 0._EB
            CASE( 3); M%FZ(IIG,JJG,KKG-1,0) = 0._EB
            CASE(-3); M%FZ(IIG,JJG,KKG,0)   = 0._EB
         END SELECT
      ELSE
         IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN
            ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
            CALL GET_MOLECULAR_WEIGHT(ZZ_GET,MW_F)
            SELECT CASE(IOR)
               CASE( 1); M%FX(IIG-1,JJG,KKG,0) = B1%RHO_F/MW_F
               CASE(-1); M%FX(IIG,JJG,KKG,0)   = B1%RHO_F/MW_F
               CASE( 2); M%FY(IIG,JJG-1,KKG,0) = B1%RHO_F/MW_F
               CASE(-2); M%FY(IIG,JJG,KKG,0)   = B1%RHO_F/MW_F
               CASE( 3); M%FZ(IIG,JJG,KKG-1,0) = B1%RHO_F/MW_F
               CASE(-3); M%FZ(IIG,JJG,KKG,0)   = B1%RHO_F/MW_F
            END SELECT
         ENDIF
      ENDIF

      ! Overwrite first off-wall advective flux if flow is away from the wall and if the face is not also a wall cell

      OFF_WALL_IF_3: IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .AND. WC%BOUNDARY_TYPE/=OPEN_BOUNDARY) THEN

         OFF_WALL_SELECT_3: SELECT CASE(IOR)
            CASE( 1) OFF_WALL_SELECT_3
               !      ghost          FX/UU(II+1)
               ! ///   II   ///  II+1  |  II+2  | ...
               !                       ^ WALL_INDEX(II+1,+1)
               IF ((UU(II+1,JJ,KK)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II+1,JJ,KK))%WALL_INDEX(+1)>0)) THEN
                  Z_TEMP(0:3,1,1) = (/RHO_RMW(II+1,JJ,KK),RHO_RMW(II+1:II+2,JJ,KK),DUMMY/)
                  U_TEMP(1,1,1) = UU(II+1,JJ,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  M%FX(II+1,JJ,KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-1) OFF_WALL_SELECT_3
               !            FX/UU(II-2)     ghost
               ! ... |  II-2  |  II-1  ///   II   ///
               !              ^ WALL_INDEX(II-1,-1)
               IF ((UU(II-2,JJ,KK)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II-1,JJ,KK))%WALL_INDEX(-1)>0)) THEN
                  Z_TEMP(0:3,1,1) = (/DUMMY,RHO_RMW(II-2:II-1,JJ,KK),RHO_RMW(II-1,JJ,KK)/)
                  U_TEMP(1,1,1) = UU(II-2,JJ,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
                  M%FX(II-2,JJ,KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE( 2) OFF_WALL_SELECT_3
               IF ((VV(II,JJ+1,KK)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ+1,KK))%WALL_INDEX(+2)>0)) THEN
                  Z_TEMP(1,0:3,1) = (/RHO_RMW(II,JJ+1,KK),RHO_RMW(II,JJ+1:JJ+2,KK),DUMMY/)
                  U_TEMP(1,1,1) = VV(II,JJ+1,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  M%FY(II,JJ+1,KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-2) OFF_WALL_SELECT_3
               IF ((VV(II,JJ-2,KK)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ-1,KK))%WALL_INDEX(-2)>0)) THEN
                  Z_TEMP(1,0:3,1) = (/DUMMY,RHO_RMW(II,JJ-2:JJ-1,KK),RHO_RMW(II,JJ-1,KK)/)
                  U_TEMP(1,1,1) = VV(II,JJ-2,KK)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,2,I_FLUX_LIMITER)
                  M%FY(II,JJ-2,KK,0) = F_TEMP(1,1,1)
               ENDIF
            CASE( 3) OFF_WALL_SELECT_3
               IF ((WW(II,JJ,KK+1)>0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ,KK+1))%WALL_INDEX(+3)>0)) THEN
                  Z_TEMP(1,1,0:3) = (/RHO_RMW(II,JJ,KK+1),RHO_RMW(II,JJ,KK+1:KK+2),DUMMY/)
                  U_TEMP(1,1,1) = WW(II,JJ,KK+1)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  M%FZ(II,JJ,KK+1,0) = F_TEMP(1,1,1)
               ENDIF
            CASE(-3) OFF_WALL_SELECT_3
               IF ((WW(II,JJ,KK-2)<0._EB) .AND. .NOT.(M%CELL(M%CELL_INDEX(II,JJ,KK-1))%WALL_INDEX(-3)>0)) THEN
                  Z_TEMP(1,1,0:3) = (/DUMMY,RHO_RMW(II,JJ,KK-2:KK-1),RHO_RMW(II,JJ,KK-1)/)
                  U_TEMP(1,1,1) = WW(II,JJ,KK-2)
                  CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,3,I_FLUX_LIMITER)
                  M%FZ(II,JJ,KK-2,0) = F_TEMP(1,1,1)
               ENDIF
         END SELECT OFF_WALL_SELECT_3

      ENDIF OFF_WALL_IF_3

   ENDDO WALL_LOOP_3

   ! Now correct the max face value of (RHO*ZZ) such that SUM(RHO*ZZ/MW)_FACE = RHO_FACE/MW_FACE
   ! (necessary condition to preserve isothermal flow)

   DO K=0,M%KBAR
      DO J=0,M%JBAR
         DO I=0,M%IBAR
            N=MAXLOC(M%FX(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            M%FX(I,J,K,N) = MW_G*MAX( 0._EB, M%FX(I,J,K,0) &
                                           - SUM(M%FX(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                           - SUM(M%FX(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )

            N=MAXLOC(M%FY(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            M%FY(I,J,K,N) = MW_G*MAX( 0._EB, M%FY(I,J,K,0) &
                                           - SUM(M%FY(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                           - SUM(M%FY(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )

            N=MAXLOC(M%FZ(I,J,K,1:N_TRACKED_SPECIES),1)
            MW_G = SPECIES_MIXTURE(N)%MW
            M%FZ(I,J,K,N) = MW_G*MAX( 0._EB, M%FZ(I,J,K,0) &
                                           - SUM(M%FZ(I,J,K,1:(N-1))/SPECIES_MIXTURE(1:(N-1))%MW) &
                                           - SUM(M%FZ(I,J,K,(N+1):N_TRACKED_SPECIES)/SPECIES_MIXTURE((N+1):N_TRACKED_SPECIES)%MW) )
         ENDDO
      ENDDO
   ENDDO

ENDIF FACE_CORRECTION_IF

END SUBROUTINE MASS_FINITE_DIFFERENCES_NEW_KERNEL


!> \brief Update the species mass fractions and density
!> \param M Mesh data structure
!> \param T Simulation time (s)
!> \param DT Time step (s)
!> \param NM Mesh index

SUBROUTINE DENSITY_KERNEL(M,T,DT,NM,WORK_BRANCH)

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_Z_OF_RHO,VD2D_MMS_Z_SRC,UF_MMS,WF_MMS,VD2D_MMS_RHO_OF_Z,VD2D_MMS_Z_SRC
USE SOOT_ROUTINES, ONLY: SETTLING_VELOCITY
USE CC_VERIFICATION, ONLY : ROTATED_CUBE_RHS_ZZ
USE CC_DIVERGENCE_KERNELS, ONLY : SET_EXIMADVFLX_3D

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN), OPTIONAL :: WORK_BRANCH
REAL(EB) :: RHS,Q_Z,XHAT,ZHAT
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
INTEGER :: I,J,K,N,IW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: DEL_RHO_D_DEL_Z__0
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
REAL(EB), POINTER, DIMENSION(:,:,:) :: WK4,WK5

IF (SOLID_PHASE_ONLY) RETURN

! Select WORK scratch arrays based on branch
IF (PRESENT(WORK_BRANCH) .AND. WORK_BRANCH==2) THEN
   WK4=>M%WORK4_B; WK5=>M%WORK5_B
ELSE
   WK4=>M%WORK4; WK5=>M%WORK5
ENDIF

! If the RHS of the continuity equation does not yet satisfy the divergence constraint, return.
! This is typical of the case where an initial velocity field is specified by the user.

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
END SELECT

UU=>M%WORK_U
VV=>M%WORK_V
WW=>M%WORK_W
DEL_RHO_D_DEL_Z__0=>M%SWORK4

PREDICTOR_STEP: SELECT CASE (PREDICTOR)

CASE(.TRUE.) PREDICTOR_STEP

   IF (FIRST_PASS) THEN
      ! This IF is required because DEL_RHO_D_DEL_Z is updated to the next time level in divg within
      ! the CHANGE_TIME_STEP loop in main while we are determining the appropriate stable DT.
      IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. (GRAVITATIONAL_SETTLING .OR. THERMOPHORETIC_SETTLING)) CALL SETTLING_VELOCITY(NM)
      DEL_RHO_D_DEL_Z__0 = M%DEL_RHO_D_DEL_Z
   ENDIF

   ! Correct boundary velocity at wall cells

   UU=M%U
   VV=M%V
   WW=M%W


   WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP
      BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
      SELECT CASE(BC%IOR)
         CASE( 1); UU(BC%IIG-1,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-1); UU(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 2); VV(BC%IIG  ,BC%JJG-1,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-2); VV(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 3); WW(BC%IIG  ,BC%JJG  ,BC%KKG-1) = M%UVW_SAVE(IW)
         CASE(-3); WW(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
      END SELECT
   ENDDO WALL_LOOP

   ! Predictor step for mass density

   DO N=1,N_TOTAL_SCALARS
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
               RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) &
                   + (M%FX(I,J,K,N)*UU(I,J,K)*M%R(I) - M%FX(I-1,J,K,N)*UU(I-1,J,K)*M%R(I-1))*M%RDX(I)*M%RRN(I) &
                   + (M%FY(I,J,K,N)*VV(I,J,K)        - M%FY(I,J-1,K,N)*VV(I,J-1,K)           )*M%RDY(J)           &
                   + (M%FZ(I,J,K,N)*WW(I,J,K)        - M%FZ(I,J,K-1,N)*WW(I,J,K-1)           )*M%RDZ(K)
               M%ZZS(I,J,K,N) = M%RHO(I,J,K)*M%ZZ(I,J,K,N) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO


   IF (CC_IBM) CALL SET_EXIMADVFLX_3D(M,UU,VV,WW)
   IF (STORE_SPECIES_FLUX) THEN
      DO N=1,N_TOTAL_SCALARS
         DO K=0,M%KBAR
            DO J=0,M%JBAR
               DO I=0,M%IBAR
                  M%ADV_FX(I,J,K,N) = M%FX(I,J,K,N)*UU(I,J,K)
                  M%ADV_FY(I,J,K,N) = M%FY(I,J,K,N)*VV(I,J,K)
                  M%ADV_FZ(I,J,K,N) = M%FZ(I,J,K,N)*WW(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Add gas production source term

   IF (ALLOCATED(M%M_DOT_PPP)) M%ZZS(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES) = &
                                        M%ZZS(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES) + &
                                        DT*M%M_DOT_PPP(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES)

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               ! divergence from EOS
               XHAT = M%XC(I) - UF_MMS*T
               ZHAT = M%ZC(K) - WF_MMS*T
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               M%ZZS(I,J,K,1) = M%ZZS(I,J,K,1) - DT*Q_Z
               M%ZZS(I,J,K,2) = M%ZZS(I,J,K,2) + DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ELSEIF(PERIODIC_TEST==21 .OR. PERIODIC_TEST==22 .OR. PERIODIC_TEST==23) THEN
      CALL ROTATED_CUBE_RHS_ZZ(T,DT,NM)
   ENDIF


   ! Get rho = sum(rho*Y_alpha)

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%RHOS(I,J,K) = SUM(M%ZZS(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

   ! Check mass density for positivity

   CALL CHECK_MASS_DENSITY

   ALLOCATE(ZZ_GET(1:N_TOTAL_SCALARS))

   ! Extract mass fraction from RHO * ZZ

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%ZZS(I,J,K,1:N_TOTAL_SCALARS) = M%ZZS(I,J,K,1:N_TOTAL_SCALARS)/M%RHOS(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   ! Passive scalars

   CALL CLIP_PASSIVE_SCALARS

   ! Predict background pressure at next time step

   DO I=1,N_ZONE
      M%PBAR_S(:,I) = M%PBAR(:,I) + M%D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = M%ZZS(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,M%RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%TMP(I,J,K) = M%PBAR_S(K,M%PRESSURE_ZONE(I,J,K))/(M%RSUM(I,J,K)*M%RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(ZZ_GET)


CASE(.FALSE.) PREDICTOR_STEP  ! CORRECTOR step

   ! Correct boundary velocity at wall cells

   UU=M%US
   VV=M%VS
   WW=M%WS


   WALL_LOOP_2: DO IW=1,M%N_EXTERNAL_WALL_CELLS
      WC => M%WALL(IW)
      EWC => M%EXTERNAL_WALL(IW)
      IF (EWC%BOUNDARY_TYPE_PREVIOUS/=INTERPOLATED_BOUNDARY) CYCLE WALL_LOOP_2
      BC => M%BOUNDARY_COORD(WC%BC_INDEX)
      SELECT CASE(BC%IOR)
         CASE( 1); UU(BC%IIG-1,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-1); UU(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 2); VV(BC%IIG  ,BC%JJG-1,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-2); VV(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 3); WW(BC%IIG  ,BC%JJG  ,BC%KKG-1) = M%UVW_SAVE(IW)
         CASE(-3); WW(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
      END SELECT
   ENDDO WALL_LOOP_2

   IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. (GRAVITATIONAL_SETTLING .OR. THERMOPHORETIC_SETTLING)) CALL SETTLING_VELOCITY(NM)

   ! Compute species mass density at the next time step

   DO N=1,N_TOTAL_SCALARS
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
               RHS = - M%DEL_RHO_D_DEL_Z(I,J,K,N) &
                   + (M%FX(I,J,K,N)*UU(I,J,K)*M%R(I) - M%FX(I-1,J,K,N)*UU(I-1,J,K)*M%R(I-1))*M%RDX(I)*M%RRN(I) &
                   + (M%FY(I,J,K,N)*VV(I,J,K)        - M%FY(I,J-1,K,N)*VV(I,J-1,K)           )*M%RDY(J)           &
                   + (M%FZ(I,J,K,N)*WW(I,J,K)        - M%FZ(I,J,K-1,N)*WW(I,J,K-1)           )*M%RDZ(K)
               M%ZZ(I,J,K,N) = .5_EB*( M%RHO(I,J,K)*M%ZZ(I,J,K,N) + M%RHOS(I,J,K)*M%ZZS(I,J,K,N) - DT*RHS )
            ENDDO
         ENDDO
      ENDDO
   ENDDO


   ! Add gas production source term

   IF (ALLOCATED(M%M_DOT_PPP)) THEN
      M%ZZ(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES) = M%ZZ(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES) + &
                                                     0.5_EB*DT*M%M_DOT_PPP(0:M%IBP1,0:M%JBP1,0:M%KBP1,1:N_TRACKED_SPECIES)
      IF (.NOT. CC_IBM) THEN ! We will use these for Regular cells in cut-cell region in CC_DENSITY.
         M%M_DOT_PPP = 0._EB
         M%D_SOURCE  = 0._EB
      ENDIF
   ENDIF

   IF (CC_IBM) CALL SET_EXIMADVFLX_3D(M,UU,VV,WW)
   IF (STORE_SPECIES_FLUX) THEN
      DO N=1,N_TOTAL_SCALARS
         DO K=0,M%KBAR
            DO J=0,M%JBAR
               DO I=0,M%IBAR
                  M%ADV_FX(I,J,K,N) = 0.5_EB*( M%ADV_FX(I,J,K,N) + M%FX(I,J,K,N)*UU(I,J,K) )
                  M%ADV_FY(I,J,K,N) = 0.5_EB*( M%ADV_FY(I,J,K,N) + M%FY(I,J,K,N)*VV(I,J,K) )
                  M%ADV_FZ(I,J,K,N) = 0.5_EB*( M%ADV_FZ(I,J,K,N) + M%FZ(I,J,K,N)*WW(I,J,K) )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Manufactured solution

   IF (PERIODIC_TEST==7) THEN
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               ! divergence from EOS
               XHAT = M%XC(I) - UF_MMS*T
               ZHAT = M%ZC(K) - WF_MMS*T
               Q_Z = VD2D_MMS_Z_SRC(XHAT,ZHAT,T)
               M%ZZ(I,J,K,1) = M%ZZ(I,J,K,1) - .5_EB*DT*Q_Z
               M%ZZ(I,J,K,2) = M%ZZ(I,J,K,2) + .5_EB*DT*Q_Z
            ENDDO
         ENDDO
      ENDDO
   ELSEIF(PERIODIC_TEST==21 .OR. PERIODIC_TEST==22 .OR. PERIODIC_TEST==23) THEN
      CALL ROTATED_CUBE_RHS_ZZ(T,DT,NM)
   ENDIF

   ! Get rho = sum(rho*Y_alpha)


   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%RHO(I,J,K) = SUM(M%ZZ(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

   ! Check mass density for positivity

   CALL CHECK_MASS_DENSITY

   ALLOCATE(ZZ_GET(1:N_TOTAL_SCALARS))

   ! Extract Y_n from rho*Y_n

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%ZZ(I,J,K,1:N_TOTAL_SCALARS) = M%ZZ(I,J,K,1:N_TOTAL_SCALARS)/M%RHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   ! Passive scalars

   CALL CLIP_PASSIVE_SCALARS

   ! Correct background pressure

   DO I=1,N_ZONE
      M%PBAR(:,I) = 0.5_EB*(M%PBAR(:,I) + M%PBAR_S(:,I) + M%D_PBAR_DT_S(I)*DT)
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = M%ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,M%RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%TMP(I,J,K) = M%PBAR(K,M%PRESSURE_ZONE(I,J,K))/(M%RSUM(I,J,K)*M%RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(ZZ_GET)


END SELECT PREDICTOR_STEP

CONTAINS

!> \brief Redistribute mass from cells below or above the density cut-off limits
!> \details Do not apply OpenMP to this routine

SUBROUTINE CHECK_MASS_DENSITY

REAL(EB) :: MASS_N(-3:3),CONST,MASS_C,RHO_ZZ_CUT,RHO_CUT,VC(-3:3),SIGN_FACTOR,SUM_MASS_N,VC1(-3:3),&
            RHO_ZZ_MIN,RHO_ZZ_MAX,SUM_RHO_ZZ,RHO_ZZ_TEST
INTEGER :: IC
LOGICAL :: CLIP_RHO_ZZ(N_TRACKED_SPECIES)
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_RHO,DELTA_RHO_ZZ,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ

DELTA_RHO => WK4
DELTA_RHO =  0._EB
M%CLIP_RHOMIN = .FALSE.
M%CLIP_RHOMAX = .FALSE.

IF (PREDICTOR) THEN
   RHO_ZZ => M%ZZS  ! At this stage of the time step, ZZS is actually RHOS*ZZS
   RHOP   => M%RHOS
ELSE
   RHO_ZZ => M%ZZ   ! At this stage of the time step, ZZ is actually RHO*ZZ
   RHOP   => M%RHO
ENDIF

! Correct density

DO K=1,M%KBAR
   DO J=1,M%JBAR
      VC1( 0)  = M%DY(J)  *M%DZ(K)
      VC1(-1)  = VC1( 0)
      VC1( 1)  = VC1( 0)
      VC1(-2)  = M%DY(J-1)*M%DZ(K)
      VC1( 2)  = M%DY(J+1)*M%DZ(K)
      VC1(-3)  = M%DY(J)  *M%DZ(K-1)
      VC1( 3)  = M%DY(J)  *M%DZ(K+1)
      DO I=1,M%IBAR
         IF (RHOP(I,J,K)>=RHOMIN .AND. RHOP(I,J,K)<=RHOMAX) CYCLE
         IC = M%CELL_INDEX(I,J,K)
         IF (M%CELL(IC)%SOLID) CYCLE
         IF (RHOP(I,J,K)<RHOMIN) THEN
            RHO_CUT = RHOMIN
            SIGN_FACTOR = 1._EB
            M%CLIP_RHOMIN = .TRUE.
         ELSE
            RHO_CUT = RHOMAX
            SIGN_FACTOR = -1._EB
            M%CLIP_RHOMAX = .TRUE.
         ENDIF
         MASS_N = 0._EB
         VC( 0)  = M%DX(I)  * VC1( 0)
         VC(-1)  = M%DX(I-1)* VC1(-1)
         VC( 1)  = M%DX(I+1)* VC1( 1)
         VC(-2)  = M%DX(I)  * VC1(-2)
         VC( 2)  = M%DX(I)  * VC1( 2)
         VC(-3)  = M%DX(I)  * VC1(-3)
         VC( 3)  = M%DX(I)  * VC1( 3)

         MASS_C = ABS(RHO_CUT-RHOP(I,J,K))*VC(0)
         IF (M%CELL(IC)%WALL_INDEX(-1)==0) MASS_N(-1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I-1,J,K)))-RHO_CUT)*VC(-1)
         IF (M%CELL(IC)%WALL_INDEX( 1)==0) MASS_N( 1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I+1,J,K)))-RHO_CUT)*VC( 1)
         IF (M%CELL(IC)%WALL_INDEX(-2)==0) MASS_N(-2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J-1,K)))-RHO_CUT)*VC(-2)
         IF (M%CELL(IC)%WALL_INDEX( 2)==0) MASS_N( 2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J+1,K)))-RHO_CUT)*VC( 2)
         IF (M%CELL(IC)%WALL_INDEX(-3)==0) MASS_N(-3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K-1)))-RHO_CUT)*VC(-3)
         IF (M%CELL(IC)%WALL_INDEX( 3)==0) MASS_N( 3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K+1)))-RHO_CUT)*VC( 3)
         SUM_MASS_N = SUM(MASS_N)
         IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
         DELTA_RHO(I,J,K)   = DELTA_RHO(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
         DELTA_RHO(I-1,J,K) = DELTA_RHO(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
         DELTA_RHO(I+1,J,K) = DELTA_RHO(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
         DELTA_RHO(I,J-1,K) = DELTA_RHO(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
         DELTA_RHO(I,J+1,K) = DELTA_RHO(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
         DELTA_RHO(I,J,K-1) = DELTA_RHO(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
         DELTA_RHO(I,J,K+1) = DELTA_RHO(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
      ENDDO
   ENDDO
ENDDO

! Assign excess/deficit mass to neighboring cells if clipping has been done.

IF (M%CLIP_RHOMIN .OR. M%CLIP_RHOMAX) &
   RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR) = MIN(RHOMAX,MAX(RHOMIN,RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR)+DELTA_RHO(1:M%IBAR,1:M%JBAR,1:M%KBAR)))

! If there is only one gas species, set rho*Z=rho and return.

IF (N_TRACKED_SPECIES==1) THEN
   IF (M%CLIP_RHOMIN .OR. M%CLIP_RHOMAX) RHO_ZZ(1:M%IBAR,1:M%JBAR,1:M%KBAR,1) = RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR)
   RETURN
ENDIF

! Correct species mass density

RHO_ZZ_MIN = 0._EB
CLIP_RHO_ZZ = .FALSE.

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   DELTA_RHO_ZZ => WK5
   DELTA_RHO_ZZ = 0._EB

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         VC1( 0)  = M%DY(J)  *M%DZ(K)
         VC1(-1)  = VC1( 0)
         VC1( 1)  = VC1( 0)
         VC1(-2)  = M%DY(J-1)*M%DZ(K)
         VC1( 2)  = M%DY(J+1)*M%DZ(K)
         VC1(-3)  = M%DY(J)  *M%DZ(K-1)
         VC1( 3)  = M%DY(J)  *M%DZ(K+1)
         DO I=1,M%IBAR

            IC = M%CELL_INDEX(I,J,K)
            IF (M%CELL(IC)%SOLID) CYCLE

            RHO_ZZ_MAX = RHOP(I,J,K)
            IF (RHO_ZZ(I,J,K,N)>=RHO_ZZ_MIN .AND. RHO_ZZ(I,J,K,N)<=RHO_ZZ_MAX) CYCLE
            CLIP_RHO_ZZ(N) = .TRUE.
            IF (RHO_ZZ(I,J,K,N)<RHO_ZZ_MIN) THEN
               RHO_ZZ_CUT = RHO_ZZ_MIN
               SIGN_FACTOR = 1._EB
            ELSE
               RHO_ZZ_CUT = RHO_ZZ_MAX
               SIGN_FACTOR = -1._EB
            ENDIF
            MASS_N = 0._EB
            VC( 0)  = M%DX(I)  * VC1( 0)
            VC(-1)  = M%DX(I-1)* VC1(-1)
            VC( 1)  = M%DX(I+1)* VC1( 1)
            VC(-2)  = M%DX(I)  * VC1(-2)
            VC( 2)  = M%DX(I)  * VC1( 2)
            VC(-3)  = M%DX(I)  * VC1(-3)
            VC( 3)  = M%DX(I)  * VC1( 3)

            MASS_C = ABS(RHO_ZZ_CUT-RHO_ZZ(I,J,K,N))*VC(0)
            IF (M%CELL(IC)%WALL_INDEX(-1)==0) MASS_N(-1) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I-1,J,K,N)))-RHO_ZZ_CUT)*VC(-1)
            IF (M%CELL(IC)%WALL_INDEX( 1)==0) MASS_N( 1) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I+1,J,K,N)))-RHO_ZZ_CUT)*VC( 1)
            IF (M%CELL(IC)%WALL_INDEX(-2)==0) MASS_N(-2) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J-1,K,N)))-RHO_ZZ_CUT)*VC(-2)
            IF (M%CELL(IC)%WALL_INDEX( 2)==0) MASS_N( 2) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J+1,K,N)))-RHO_ZZ_CUT)*VC( 2)
            IF (M%CELL(IC)%WALL_INDEX(-3)==0) MASS_N(-3) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K-1,N)))-RHO_ZZ_CUT)*VC(-3)
            IF (M%CELL(IC)%WALL_INDEX( 3)==0) MASS_N( 3) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K+1,N)))-RHO_ZZ_CUT)*VC( 3)
            SUM_MASS_N = SUM(MASS_N)
            IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
            CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
            DELTA_RHO_ZZ(I,J,K)   = DELTA_RHO_ZZ(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
            DELTA_RHO_ZZ(I-1,J,K) = DELTA_RHO_ZZ(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
            DELTA_RHO_ZZ(I+1,J,K) = DELTA_RHO_ZZ(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
            DELTA_RHO_ZZ(I,J-1,K) = DELTA_RHO_ZZ(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
            DELTA_RHO_ZZ(I,J+1,K) = DELTA_RHO_ZZ(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
            DELTA_RHO_ZZ(I,J,K-1) = DELTA_RHO_ZZ(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
            DELTA_RHO_ZZ(I,J,K+1) = DELTA_RHO_ZZ(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
         ENDDO
      ENDDO
   ENDDO

   IF (.NOT.CLIP_RHO_ZZ(N)) CYCLE

   ! Assign excess/deficit RHO_ZZ neighboring cells

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RHO_ZZ(I,J,K,N) = MIN(RHOP(I,J,K),MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K,N)+DELTA_RHO_ZZ(I,J,K)))
         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

! If nothing has been clipped, return

IF (.NOT.M%CLIP_RHOMIN .AND. .NOT.M%CLIP_RHOMAX .AND. .NOT. ANY(CLIP_RHO_ZZ)) RETURN

! Final check of RHO_ZZ to ensure that ZZ(:,:,:,1:N_TRACKED_SPECIES) sums to 1

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
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

END SUBROUTINE CHECK_MASS_DENSITY


SUBROUTINE CLIP_PASSIVE_SCALARS

! Currently only set up for unmixed fraction, zeta

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

IF (N_PASSIVE_SCALARS==0) RETURN

IF (PREDICTOR) THEN
   ZZP=>M%ZZS
ELSE
   ZZP=>M%ZZ
ENDIF

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
         ZZP(I,J,K,ZETA_INDEX) = MAX(0._EB,MIN(1._EB,ZZP(I,J,K,ZETA_INDEX)))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CLIP_PASSIVE_SCALARS

END SUBROUTINE DENSITY_KERNEL


!> \brief Sequential preprocessing for density block decomposition.
!> Sets up work arrays (UU/VV/WW), handles wall boundary corrections,
!> and calls SETTLING_VELOCITY. Must run sequentially per mesh before
!> K-block parallel execution.
!> \param M Mesh data structure
!> \param T Current time
!> \param DT Time step
!> \param NM Mesh index

RECURSIVE SUBROUTINE DENSITY_BLOCK_PREPROCESSING(M,T,DT,NM)

USE SOOT_ROUTINES, ONLY: SETTLING_VELOCITY

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
INTEGER :: IW
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: DEL_RHO_D_DEL_Z__0
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

IF (SOLID_PHASE_ONLY) RETURN

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
END SELECT

UU=>M%WORK_U
VV=>M%WORK_V
WW=>M%WORK_W
DEL_RHO_D_DEL_Z__0=>M%SWORK4

IF (PREDICTOR) THEN

   IF (FIRST_PASS) THEN
      IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. (GRAVITATIONAL_SETTLING .OR. THERMOPHORETIC_SETTLING)) CALL SETTLING_VELOCITY(NM)
      DEL_RHO_D_DEL_Z__0 = M%DEL_RHO_D_DEL_Z
   ENDIF

   UU=M%U
   VV=M%V
   WW=M%W

   DO IW=1,M%N_EXTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
      BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
      SELECT CASE(BC%IOR)
         CASE( 1); UU(BC%IIG-1,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-1); UU(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 2); VV(BC%IIG  ,BC%JJG-1,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-2); VV(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 3); WW(BC%IIG  ,BC%JJG  ,BC%KKG-1) = M%UVW_SAVE(IW)
         CASE(-3); WW(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
      END SELECT
   ENDDO

ELSE ! CORRECTOR

   IF (ANY(SPECIES_MIXTURE%DEPOSITING) .AND. (GRAVITATIONAL_SETTLING .OR. THERMOPHORETIC_SETTLING)) CALL SETTLING_VELOCITY(NM)

   UU=M%US
   VV=M%VS
   WW=M%WS

   DO IW=1,M%N_EXTERNAL_WALL_CELLS
      WC => M%WALL(IW)
      EWC => M%EXTERNAL_WALL(IW)
      IF (EWC%BOUNDARY_TYPE_PREVIOUS/=INTERPOLATED_BOUNDARY) CYCLE
      BC => M%BOUNDARY_COORD(WC%BC_INDEX)
      SELECT CASE(BC%IOR)
         CASE( 1); UU(BC%IIG-1,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-1); UU(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 2); VV(BC%IIG  ,BC%JJG-1,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE(-2); VV(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
         CASE( 3); WW(BC%IIG  ,BC%JJG  ,BC%KKG-1) = M%UVW_SAVE(IW)
         CASE(-3); WW(BC%IIG  ,BC%JJG  ,BC%KKG  ) = M%UVW_SAVE(IW)
      END SELECT
   ENDDO

ENDIF

END SUBROUTINE DENSITY_BLOCK_PREPROCESSING


!> \brief K-block parallel kernel for density computation.
!> Computes species mass density and total density for cells in [K1,K2].
!> Must be called after DENSITY_BLOCK_PREPROCESSING sets up work arrays.
!> \param M Mesh data structure
!> \param T Current time
!> \param DT Time step
!> \param NM Mesh index
!> \param K1 Start of K-range (1-based inclusive)
!> \param K2 End of K-range (1-based inclusive)

RECURSIVE SUBROUTINE DENSITY_BLOCK_KERNEL_COMPUTE(M,T,DT,NM,K1,K2)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM,K1,K2
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: RHS
INTEGER :: I,J,K,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: DEL_RHO_D_DEL_Z__0

IF (SOLID_PHASE_ONLY) RETURN

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
END SELECT

UU=>M%WORK_U
VV=>M%WORK_V
WW=>M%WORK_W
DEL_RHO_D_DEL_Z__0=>M%SWORK4

PREDICTOR_STEP: IF (PREDICTOR) THEN

   ! Species mass density (K1:K2)
   DO N=1,N_TOTAL_SCALARS
      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
               RHS = - DEL_RHO_D_DEL_Z__0(I,J,K,N) &
                   + (M%FX(I,J,K,N)*UU(I,J,K)*M%R(I) - M%FX(I-1,J,K,N)*UU(I-1,J,K)*M%R(I-1))*M%RDX(I)*M%RRN(I) &
                   + (M%FY(I,J,K,N)*VV(I,J,K)        - M%FY(I,J-1,K,N)*VV(I,J-1,K)           )*M%RDY(J)           &
                   + (M%FZ(I,J,K,N)*WW(I,J,K)        - M%FZ(I,J,K-1,N)*WW(I,J,K-1)           )*M%RDZ(K)
               M%ZZS(I,J,K,N) = M%RHO(I,J,K)*M%ZZ(I,J,K,N) - DT*RHS
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Add gas production source term (K1:K2, interior only)
   IF (ALLOCATED(M%M_DOT_PPP)) THEN
      DO N=1,N_TRACKED_SPECIES
         DO K=K1,K2
            DO J=1,M%JBAR
               DO I=1,M%IBAR
                  M%ZZS(I,J,K,N) = M%ZZS(I,J,K,N) + DT*M%M_DOT_PPP(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Get rho = sum(rho*Y_alpha) (K1:K2)
   DO K=K1,K2
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%RHOS(I,J,K) = SUM(M%ZZS(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

ELSE PREDICTOR_STEP ! CORRECTOR

   ! Species mass density (K1:K2)
   DO N=1,N_TOTAL_SCALARS
      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
               RHS = - M%DEL_RHO_D_DEL_Z(I,J,K,N) &
                   + (M%FX(I,J,K,N)*UU(I,J,K)*M%R(I) - M%FX(I-1,J,K,N)*UU(I-1,J,K)*M%R(I-1))*M%RDX(I)*M%RRN(I) &
                   + (M%FY(I,J,K,N)*VV(I,J,K)        - M%FY(I,J-1,K,N)*VV(I,J-1,K)           )*M%RDY(J)           &
                   + (M%FZ(I,J,K,N)*WW(I,J,K)        - M%FZ(I,J,K-1,N)*WW(I,J,K-1)           )*M%RDZ(K)
               M%ZZ(I,J,K,N) = .5_EB*( M%RHO(I,J,K)*M%ZZ(I,J,K,N) + M%RHOS(I,J,K)*M%ZZS(I,J,K,N) - DT*RHS )
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Add gas production source term (K1:K2, interior only)
   IF (ALLOCATED(M%M_DOT_PPP)) THEN
      DO N=1,N_TRACKED_SPECIES
         DO K=K1,K2
            DO J=1,M%JBAR
               DO I=1,M%IBAR
                  M%ZZ(I,J,K,N) = M%ZZ(I,J,K,N) + 0.5_EB*DT*M%M_DOT_PPP(I,J,K,N)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Get rho = sum(rho*Y_alpha) (K1:K2)
   DO K=K1,K2
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%RHO(I,J,K) = SUM(M%ZZ(I,J,K,1:N_TRACKED_SPECIES))
         ENDDO
      ENDDO
   ENDDO

ENDIF PREDICTOR_STEP

END SUBROUTINE DENSITY_BLOCK_KERNEL_COMPUTE


!> \brief Sequential postprocessing for density block decomposition.
!> Runs CHECK_MASS_DENSITY, extracts mass fractions, updates pressure/temperature.
!> Must run sequentially per mesh after all K-blocks complete.
!> \param M Mesh data structure
!> \param T Current time
!> \param DT Time step
!> \param NM Mesh index

RECURSIVE SUBROUTINE DENSITY_BLOCK_POSTPROCESSING(M,T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY : GET_SPECIFIC_GAS_CONSTANT
USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_Z_OF_RHO,VD2D_MMS_Z_SRC,UF_MMS,WF_MMS,VD2D_MMS_RHO_OF_Z,VD2D_MMS_Z_SRC

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
INTEGER :: I,J,K,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW

IF (SOLID_PHASE_ONLY) RETURN

SELECT CASE (PERIODIC_TEST)
   CASE DEFAULT
      IF (ICYC<=1) RETURN
   CASE (5,8)
      RETURN
   CASE (4,7,11,21,22)
END SELECT

UU=>M%WORK_U
VV=>M%WORK_V
WW=>M%WORK_W

IF (PREDICTOR) THEN

   IF (STORE_SPECIES_FLUX) THEN
      DO N=1,N_TOTAL_SCALARS
         DO K=0,M%KBAR
            DO J=0,M%JBAR
               DO I=0,M%IBAR
                  M%ADV_FX(I,J,K,N) = M%FX(I,J,K,N)*UU(I,J,K)
                  M%ADV_FY(I,J,K,N) = M%FY(I,J,K,N)*VV(I,J,K)
                  M%ADV_FZ(I,J,K,N) = M%FZ(I,J,K,N)*WW(I,J,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   CALL CHECK_MASS_DENSITY_POST

   ALLOCATE(ZZ_GET(1:N_TOTAL_SCALARS))

   ! Extract mass fraction from RHO * ZZ
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%ZZS(I,J,K,1:N_TOTAL_SCALARS) = M%ZZS(I,J,K,1:N_TOTAL_SCALARS)/M%RHOS(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   CALL CLIP_PASSIVE_SCALARS_POST

   ! Predict background pressure at next time step
   DO I=1,N_ZONE
      M%PBAR_S(:,I) = M%PBAR(:,I) + M%D_PBAR_DT(I)*DT
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = M%ZZS(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,M%RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%TMP(I,J,K) = M%PBAR_S(K,M%PRESSURE_ZONE(I,J,K))/(M%RSUM(I,J,K)*M%RHOS(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(ZZ_GET)

ELSE ! CORRECTOR

   ! Clear M_DOT_PPP and D_SOURCE after corrector addition
   IF (ALLOCATED(M%M_DOT_PPP)) THEN
      M%M_DOT_PPP = 0._EB
      M%D_SOURCE  = 0._EB
   ENDIF

   IF (STORE_SPECIES_FLUX) THEN
      DO N=1,N_TOTAL_SCALARS
         DO K=0,M%KBAR
            DO J=0,M%JBAR
               DO I=0,M%IBAR
                  M%ADV_FX(I,J,K,N) = 0.5_EB*( M%ADV_FX(I,J,K,N) + M%FX(I,J,K,N)*UU(I,J,K) )
                  M%ADV_FY(I,J,K,N) = 0.5_EB*( M%ADV_FY(I,J,K,N) + M%FY(I,J,K,N)*VV(I,J,K) )
                  M%ADV_FZ(I,J,K,N) = 0.5_EB*( M%ADV_FZ(I,J,K,N) + M%FZ(I,J,K,N)*WW(I,J,K) )
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   CALL CHECK_MASS_DENSITY_POST

   ALLOCATE(ZZ_GET(1:N_TOTAL_SCALARS))

   ! Extract Y_n from rho*Y_n
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%ZZ(I,J,K,1:N_TOTAL_SCALARS) = M%ZZ(I,J,K,1:N_TOTAL_SCALARS)/M%RHO(I,J,K)
         ENDDO
      ENDDO
   ENDDO

   CALL CLIP_PASSIVE_SCALARS_POST

   ! Correct background pressure
   DO I=1,N_ZONE
      M%PBAR(:,I) = 0.5_EB*(M%PBAR(:,I) + M%PBAR_S(:,I) + M%D_PBAR_DT_S(I)*DT)
   ENDDO

   ! Compute molecular weight term RSUM=R0*SUM(Y_i/W_i)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = M%ZZ(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,M%RSUM(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   ! Extract predicted temperature at next time step from Equation of State
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            M%TMP(I,J,K) = M%PBAR(K,M%PRESSURE_ZONE(I,J,K))/(M%RSUM(I,J,K)*M%RHO(I,J,K))
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(ZZ_GET)

ENDIF

CONTAINS

!> \brief Redistribute mass from cells below or above the density cut-off limits
!> \details Do not apply OpenMP to this routine. Cross-K scatter prevents K-blocking.

SUBROUTINE CHECK_MASS_DENSITY_POST

REAL(EB) :: MASS_N(-3:3),CONST,MASS_C,RHO_ZZ_CUT,RHO_CUT,VC(-3:3),SIGN_FACTOR,SUM_MASS_N,VC1(-3:3),&
            RHO_ZZ_MIN,RHO_ZZ_MAX,SUM_RHO_ZZ,RHO_ZZ_TEST
INTEGER :: IC
LOGICAL :: CLIP_RHO_ZZ(N_TRACKED_SPECIES)
REAL(EB), POINTER, DIMENSION(:,:,:) :: DELTA_RHO,DELTA_RHO_ZZ,RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_ZZ

DELTA_RHO => M%WORK4
DELTA_RHO =  0._EB
M%CLIP_RHOMIN = .FALSE.
M%CLIP_RHOMAX = .FALSE.

IF (PREDICTOR) THEN
   RHO_ZZ => M%ZZS
   RHOP   => M%RHOS
ELSE
   RHO_ZZ => M%ZZ
   RHOP   => M%RHO
ENDIF

DO K=1,M%KBAR
   DO J=1,M%JBAR
      VC1( 0)  = M%DY(J)  *M%DZ(K)
      VC1(-1)  = VC1( 0)
      VC1( 1)  = VC1( 0)
      VC1(-2)  = M%DY(J-1)*M%DZ(K)
      VC1( 2)  = M%DY(J+1)*M%DZ(K)
      VC1(-3)  = M%DY(J)  *M%DZ(K-1)
      VC1( 3)  = M%DY(J)  *M%DZ(K+1)
      DO I=1,M%IBAR
         IF (RHOP(I,J,K)>=RHOMIN .AND. RHOP(I,J,K)<=RHOMAX) CYCLE
         IC = M%CELL_INDEX(I,J,K)
         IF (M%CELL(IC)%SOLID) CYCLE
         IF (RHOP(I,J,K)<RHOMIN) THEN
            RHO_CUT = RHOMIN
            SIGN_FACTOR = 1._EB
            M%CLIP_RHOMIN = .TRUE.
         ELSE
            RHO_CUT = RHOMAX
            SIGN_FACTOR = -1._EB
            M%CLIP_RHOMAX = .TRUE.
         ENDIF
         MASS_N = 0._EB
         VC( 0)  = M%DX(I)  * VC1( 0)
         VC(-1)  = M%DX(I-1)* VC1(-1)
         VC( 1)  = M%DX(I+1)* VC1( 1)
         VC(-2)  = M%DX(I)  * VC1(-2)
         VC( 2)  = M%DX(I)  * VC1( 2)
         VC(-3)  = M%DX(I)  * VC1(-3)
         VC( 3)  = M%DX(I)  * VC1( 3)

         MASS_C = ABS(RHO_CUT-RHOP(I,J,K))*VC(0)
         IF (M%CELL(IC)%WALL_INDEX(-1)==0) MASS_N(-1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I-1,J,K)))-RHO_CUT)*VC(-1)
         IF (M%CELL(IC)%WALL_INDEX( 1)==0) MASS_N( 1) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I+1,J,K)))-RHO_CUT)*VC( 1)
         IF (M%CELL(IC)%WALL_INDEX(-2)==0) MASS_N(-2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J-1,K)))-RHO_CUT)*VC(-2)
         IF (M%CELL(IC)%WALL_INDEX( 2)==0) MASS_N( 2) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J+1,K)))-RHO_CUT)*VC( 2)
         IF (M%CELL(IC)%WALL_INDEX(-3)==0) MASS_N(-3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K-1)))-RHO_CUT)*VC(-3)
         IF (M%CELL(IC)%WALL_INDEX( 3)==0) MASS_N( 3) = ABS(MIN(RHOMAX,MAX(RHOMIN,RHOP(I,J,K+1)))-RHO_CUT)*VC( 3)
         SUM_MASS_N = SUM(MASS_N)
         IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
         CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
         DELTA_RHO(I,J,K)   = DELTA_RHO(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
         DELTA_RHO(I-1,J,K) = DELTA_RHO(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
         DELTA_RHO(I+1,J,K) = DELTA_RHO(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
         DELTA_RHO(I,J-1,K) = DELTA_RHO(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
         DELTA_RHO(I,J+1,K) = DELTA_RHO(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
         DELTA_RHO(I,J,K-1) = DELTA_RHO(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
         DELTA_RHO(I,J,K+1) = DELTA_RHO(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
      ENDDO
   ENDDO
ENDDO

IF (M%CLIP_RHOMIN .OR. M%CLIP_RHOMAX) &
   RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR) = MIN(RHOMAX,MAX(RHOMIN,RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR) &
                                        +DELTA_RHO(1:M%IBAR,1:M%JBAR,1:M%KBAR)))

IF (N_TRACKED_SPECIES==1) THEN
   IF (M%CLIP_RHOMIN .OR. M%CLIP_RHOMAX) RHO_ZZ(1:M%IBAR,1:M%JBAR,1:M%KBAR,1) = RHOP(1:M%IBAR,1:M%JBAR,1:M%KBAR)
   RETURN
ENDIF

RHO_ZZ_MIN = 0._EB
CLIP_RHO_ZZ = .FALSE.

SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

   DELTA_RHO_ZZ => M%WORK5
   DELTA_RHO_ZZ = 0._EB

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         VC1( 0)  = M%DY(J)  *M%DZ(K)
         VC1(-1)  = VC1( 0)
         VC1( 1)  = VC1( 0)
         VC1(-2)  = M%DY(J-1)*M%DZ(K)
         VC1( 2)  = M%DY(J+1)*M%DZ(K)
         VC1(-3)  = M%DY(J)  *M%DZ(K-1)
         VC1( 3)  = M%DY(J)  *M%DZ(K+1)
         DO I=1,M%IBAR

            IC = M%CELL_INDEX(I,J,K)
            IF (M%CELL(IC)%SOLID) CYCLE

            RHO_ZZ_MAX = RHOP(I,J,K)
            IF (RHO_ZZ(I,J,K,N)>=RHO_ZZ_MIN .AND. RHO_ZZ(I,J,K,N)<=RHO_ZZ_MAX) CYCLE
            CLIP_RHO_ZZ(N) = .TRUE.
            IF (RHO_ZZ(I,J,K,N)<RHO_ZZ_MIN) THEN
               RHO_ZZ_CUT = RHO_ZZ_MIN
               SIGN_FACTOR = 1._EB
            ELSE
               RHO_ZZ_CUT = RHO_ZZ_MAX
               SIGN_FACTOR = -1._EB
            ENDIF
            MASS_N = 0._EB
            VC( 0)  = M%DX(I)  * VC1( 0)
            VC(-1)  = M%DX(I-1)* VC1(-1)
            VC( 1)  = M%DX(I+1)* VC1( 1)
            VC(-2)  = M%DX(I)  * VC1(-2)
            VC( 2)  = M%DX(I)  * VC1( 2)
            VC(-3)  = M%DX(I)  * VC1(-3)
            VC( 3)  = M%DX(I)  * VC1( 3)

            MASS_C = ABS(RHO_ZZ_CUT-RHO_ZZ(I,J,K,N))*VC(0)
            IF (M%CELL(IC)%WALL_INDEX(-1)==0) &
               MASS_N(-1) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I-1,J,K,N)))-RHO_ZZ_CUT)*VC(-1)
            IF (M%CELL(IC)%WALL_INDEX( 1)==0) &
               MASS_N( 1) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I+1,J,K,N)))-RHO_ZZ_CUT)*VC( 1)
            IF (M%CELL(IC)%WALL_INDEX(-2)==0) &
               MASS_N(-2) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J-1,K,N)))-RHO_ZZ_CUT)*VC(-2)
            IF (M%CELL(IC)%WALL_INDEX( 2)==0) &
               MASS_N( 2) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J+1,K,N)))-RHO_ZZ_CUT)*VC( 2)
            IF (M%CELL(IC)%WALL_INDEX(-3)==0) &
               MASS_N(-3) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K-1,N)))-RHO_ZZ_CUT)*VC(-3)
            IF (M%CELL(IC)%WALL_INDEX( 3)==0) &
               MASS_N( 3) = ABS(MIN(RHO_ZZ_MAX,MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K+1,N)))-RHO_ZZ_CUT)*VC( 3)
            SUM_MASS_N = SUM(MASS_N)
            IF (SUM_MASS_N<=TWO_EPSILON_EB) CYCLE
            CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
            DELTA_RHO_ZZ(I,J,K)   = DELTA_RHO_ZZ(I,J,K)   + CONST*SUM_MASS_N/VC( 0)
            DELTA_RHO_ZZ(I-1,J,K) = DELTA_RHO_ZZ(I-1,J,K) - CONST*MASS_N(-1)/VC(-1)
            DELTA_RHO_ZZ(I+1,J,K) = DELTA_RHO_ZZ(I+1,J,K) - CONST*MASS_N( 1)/VC( 1)
            DELTA_RHO_ZZ(I,J-1,K) = DELTA_RHO_ZZ(I,J-1,K) - CONST*MASS_N(-2)/VC(-2)
            DELTA_RHO_ZZ(I,J+1,K) = DELTA_RHO_ZZ(I,J+1,K) - CONST*MASS_N( 2)/VC( 2)
            DELTA_RHO_ZZ(I,J,K-1) = DELTA_RHO_ZZ(I,J,K-1) - CONST*MASS_N(-3)/VC(-3)
            DELTA_RHO_ZZ(I,J,K+1) = DELTA_RHO_ZZ(I,J,K+1) - CONST*MASS_N( 3)/VC( 3)
         ENDDO
      ENDDO
   ENDDO

   IF (.NOT.CLIP_RHO_ZZ(N)) CYCLE

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RHO_ZZ(I,J,K,N) = MIN(RHOP(I,J,K),MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K,N)+DELTA_RHO_ZZ(I,J,K)))
         ENDDO
      ENDDO
   ENDDO

ENDDO SPECIES_LOOP

IF (.NOT.M%CLIP_RHOMIN .AND. .NOT.M%CLIP_RHOMAX .AND. .NOT. ANY(CLIP_RHO_ZZ)) RETURN

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
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

END SUBROUTINE CHECK_MASS_DENSITY_POST


SUBROUTINE CLIP_PASSIVE_SCALARS_POST

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

IF (N_PASSIVE_SCALARS==0) RETURN

IF (PREDICTOR) THEN
   ZZP=>M%ZZS
ELSE
   ZZP=>M%ZZ
ENDIF

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
         ZZP(I,J,K,ZETA_INDEX) = MAX(0._EB,MIN(1._EB,ZZP(I,J,K,ZETA_INDEX)))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CLIP_PASSIVE_SCALARS_POST

END SUBROUTINE DENSITY_BLOCK_POSTPROCESSING


END MODULE MASS_KERNELS
