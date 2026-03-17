!> \brief Pure computation kernels extracted from VELO module.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE VELO_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE,BOUNDARY_PROP2_TYPE,SURFACE_TYPE,SURFACE, &
                 RAMPS_TYPE,RAMPS,OMESH_TYPE,VENTS_TYPE,EDGE_TYPE,EXTERNAL_WALL_TYPE,OBSTRUCTION_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC BAROCLINIC_CORRECTION_KERNEL,VELOCITY_PREDICTOR_KERNEL,VELOCITY_PREDICTOR_BLOCK_KERNEL, &
       VELOCITY_CORRECTOR_KERNEL,VELOCITY_CORRECTOR_BLOCK_KERNEL,VELOCITY_FLUX_KERNEL, &
       VELOCITY_FLUX_BLOCK_KERNEL, &
       COMPUTE_VISCOSITY_KERNEL,COMPUTE_VISCOSITY_BLOCK_KERNEL,COMPUTE_VISCOSITY_POST_BLOCK, &
       CHECK_STABILITY_KERNEL,VELOCITY_BC_PROCESS_EDGES_KERNEL,VISCOSITY_BC_KERNEL, &
       MATCH_VELOCITY_KERNEL,NO_FLUX_KERNEL,MATCH_VELOCITY_FLUX_KERNEL

CONTAINS


!> \brief Compute the baroclinic torque correction terms.
!> \param M Mesh data structure
!> \param T Current simulation time (s)

RECURSIVE SUBROUTINE BAROCLINIC_CORRECTION_KERNEL(M,T)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: T
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,HP,P,RRHO
INTEGER  :: I,J,K

! If the baroclinic torque term has been added to the momentum equation RHS, subtract it off.

IF (M%BAROCLINIC_TERMS_ATTACHED) THEN
   M%FVX = M%FVX - M%FVX_B
   M%FVY = M%FVY - M%FVY_B
   M%FVZ = M%FVZ - M%FVZ_B
ENDIF

P    => M%WORK1 ! p=rho*(H-K)
RRHO => M%WORK2 ! reciprocal of rho

IF (PREDICTOR) THEN
   RHOP => M%RHO
   HP   => M%H
ELSE
   RHOP => M%RHOS
   HP   => M%HS
ENDIF

! Compute pressure and 1/rho in each grid cell

DO K=0,M%KBP1
   DO J=0,M%JBP1
      DO I=0,M%IBP1
         P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-M%KRES(I,J,K))
         RRHO(I,J,K) = 1._EB/RHOP(I,J,K)
      ENDDO
   ENDDO
ENDDO

! Compute baroclinic term in the x momentum equation, p*d/dx(1/rho)

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%FVX_B(I,J,K) = -(P(I,J,K)*RHOP(I+1,J,K)+P(I+1,J,K)*RHOP(I,J,K))*(RRHO(I+1,J,K)-RRHO(I,J,K))*M%RDXN(I)/ &
                         (RHOP(I+1,J,K)+RHOP(I,J,K))
         M%FVX(I,J,K) = M%FVX(I,J,K) + M%FVX_B(I,J,K)
      ENDDO
   ENDDO
ENDDO

! Compute baroclinic term in the y momentum equation, p*d/dy(1/rho)

IF (.NOT.TWO_D) THEN
   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            M%FVY_B(I,J,K) = -(P(I,J,K)*RHOP(I,J+1,K)+P(I,J+1,K)*RHOP(I,J,K))*(RRHO(I,J+1,K)-RRHO(I,J,K))*M%RDYN(J)/ &
                            (RHOP(I,J+1,K)+RHOP(I,J,K))
            M%FVY(I,J,K) = M%FVY(I,J,K) + M%FVY_B(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute baroclinic term in the z momentum equation, p*d/dz(1/rho)

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%FVZ_B(I,J,K) = -(P(I,J,K)*RHOP(I,J,K+1)+P(I,J,K+1)*RHOP(I,J,K))*(RRHO(I,J,K+1)-RRHO(I,J,K))*M%RDZN(K)/ &
                         (RHOP(I,J,K+1)+RHOP(I,J,K))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + M%FVZ_B(I,J,K)
      ENDDO
   ENDDO
ENDDO

M%BAROCLINIC_TERMS_ATTACHED = .TRUE.

END SUBROUTINE BAROCLINIC_CORRECTION_KERNEL


!> \brief Predict the velocity components at the next time step.
!> \details Core OMP loops that compute US/VS/WS from U/V/W, FVX/FVY/FVZ, and H.
!> \param M Mesh data structure
!> \param DT Time step (s)

RECURSIVE SUBROUTINE VELOCITY_PREDICTOR_KERNEL(M,DT)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER :: I,J,K

IF (FREEZE_VELOCITY) THEN
   M%US = M%U
   M%VS = M%V
   M%WS = M%W
   RETURN
ENDIF


DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%US(I,J,K) = M%U(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(M%H(I+1,J,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%VS(I,J,K) = M%V(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(M%H(I,J+1,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%WS(I,J,K) = M%W(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(M%H(I,J,K+1)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO


END SUBROUTINE VELOCITY_PREDICTOR_KERNEL


!> \brief Block-decomposed velocity predictor: updates US, VS, WS for K-range [K1, K2].
!> \param M Mesh data structure
!> \param DT Time step (s)
!> \param K1 Start of K cell range (1-based inclusive)
!> \param K2 End of K cell range (1-based inclusive)

RECURSIVE SUBROUTINE VELOCITY_PREDICTOR_BLOCK_KERNEL(M,DT,K1,K2)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: K1,K2
INTEGER :: I,J,K,K1_W,K2_W

IF (FREEZE_VELOCITY) THEN
   M%US(0:M%IBAR,1:M%JBAR,K1:K2) = M%U(0:M%IBAR,1:M%JBAR,K1:K2)
   M%VS(1:M%IBAR,0:M%JBAR,K1:K2) = M%V(1:M%IBAR,0:M%JBAR,K1:K2)
   K1_W = K1 - 1
   K2_W = K2 - 1
   IF (K2==M%KBAR) K2_W = M%KBAR
   M%WS(1:M%IBAR,1:M%JBAR,K1_W:K2_W) = M%W(1:M%IBAR,1:M%JBAR,K1_W:K2_W)
   RETURN
ENDIF

DO K=K1,K2
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%US(I,J,K) = M%U(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(M%H(I+1,J,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

DO K=K1,K2
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%VS(I,J,K) = M%V(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(M%H(I,J+1,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

K1_W = K1 - 1
K2_W = K2 - 1
IF (K2==M%KBAR) K2_W = M%KBAR

DO K=K1_W,K2_W
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%WS(I,J,K) = M%W(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(M%H(I,J,K+1)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE VELOCITY_PREDICTOR_BLOCK_KERNEL


!> \brief Correct the velocity components at the next time step.
!> \details Core OMP loops that compute U/V/W from U/V/W, US/VS/WS, FVX/FVY/FVZ, and HS.
!> \param M Mesh data structure
!> \param DT Time step (s)

RECURSIVE SUBROUTINE VELOCITY_CORRECTOR_KERNEL(M,DT)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER :: I,J,K

IF (FREEZE_VELOCITY) THEN
   M%U = M%US
   M%V = M%VS
   M%W = M%WS
   RETURN
ENDIF


DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%U(I,J,K) = 0.5_EB*( M%U(I,J,K) + M%US(I,J,K) - DT*(M%FVX(I,J,K) + M%RDXN(I)*(M%HS(I+1,J,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%V(I,J,K) = 0.5_EB*( M%V(I,J,K) + M%VS(I,J,K) - DT*(M%FVY(I,J,K) + M%RDYN(J)*(M%HS(I,J+1,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%W(I,J,K) = 0.5_EB*( M%W(I,J,K) + M%WS(I,J,K) - DT*(M%FVZ(I,J,K) + M%RDZN(K)*(M%HS(I,J,K+1)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO


END SUBROUTINE VELOCITY_CORRECTOR_KERNEL


!> \brief Block-decomposed velocity corrector: updates U, V, W for K-range [K1, K2].
!> \details Processes only a sub-range of the K dimension, enabling intra-mesh parallelism.
!> Each block writes to non-overlapping regions of the velocity arrays.
!> For U and V (staggered in I and J): loop K=K1:K2.
!> For W (staggered in K): loop K=K1-1:K2-1, plus K=KBAR for last block.
!> \param M Mesh data structure
!> \param DT Time step (s)
!> \param K1 Start of K cell range (1-based inclusive)
!> \param K2 End of K cell range (1-based inclusive)

RECURSIVE SUBROUTINE VELOCITY_CORRECTOR_BLOCK_KERNEL(M,DT,K1,K2)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: K1,K2
INTEGER :: I,J,K,K1_W,K2_W

IF (FREEZE_VELOCITY) THEN
   M%U(0:M%IBAR,1:M%JBAR,K1:K2) = M%US(0:M%IBAR,1:M%JBAR,K1:K2)
   M%V(1:M%IBAR,0:M%JBAR,K1:K2) = M%VS(1:M%IBAR,0:M%JBAR,K1:K2)
   K1_W = K1 - 1
   K2_W = K2 - 1
   IF (K2==M%KBAR) K2_W = M%KBAR
   M%W(1:M%IBAR,1:M%JBAR,K1_W:K2_W) = M%WS(1:M%IBAR,1:M%JBAR,K1_W:K2_W)
   RETURN
ENDIF

DO K=K1,K2
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%U(I,J,K) = 0.5_EB*( M%U(I,J,K) + M%US(I,J,K) - DT*(M%FVX(I,J,K) + M%RDXN(I)*(M%HS(I+1,J,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

DO K=K1,K2
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%V(I,J,K) = 0.5_EB*( M%V(I,J,K) + M%VS(I,J,K) - DT*(M%FVY(I,J,K) + M%RDYN(J)*(M%HS(I,J+1,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

! W-faces: block owns K1-1:K2-1, last block extends to KBAR
K1_W = K1 - 1
K2_W = K2 - 1
IF (K2==M%KBAR) K2_W = M%KBAR

DO K=K1_W,K2_W
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%W(I,J,K) = 0.5_EB*( M%W(I,J,K) + M%WS(I,J,K) - DT*(M%FVZ(I,J,K) + M%RDZN(K)*(M%HS(I,J,K+1)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE VELOCITY_CORRECTOR_BLOCK_KERNEL


!> \brief Compute the velocity flux terms (vorticity, stress tensor, momentum RHS).
!> \param M Mesh data structure
!> \param T Current simulation time (s)
!> \param DT Time step (s)
!> \param NM Mesh number (needed for external calls)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables
!> \param GX Gravity component arrays (output, for use by CC_IBM in wrapper)
!> \param GY Gravity component arrays (output, for use by CC_IBM in wrapper)
!> \param GZ Gravity component arrays (output, for use by CC_IBM in wrapper)

RECURSIVE SUBROUTINE VELOCITY_FLUX_KERNEL(M,T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES,GX,GY,GZ)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY: COMPUTE_WIND_COMPONENTS
USE CC_VERIFICATION, ONLY : ROTATED_CUBE_VELOCITY_FLUX

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB), INTENT(OUT) :: GX(0:IBAR_MAX),GY(0:IBAR_MAX),GZ(0:IBAR_MAX)
REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
            RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ
INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,IC1,IC2
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ,UU,VV,WW,RHOP,DP

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
   DP => M%DS
   RHOP => M%RHOS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
   DP => M%D
   RHOP => M%RHO
ENDIF

TXY => M%WORK1
TXZ => M%WORK2
TYZ => M%WORK3
OMX => M%WORK4
OMY => M%WORK5
OMZ => M%WORK6

! Compute vorticity and stress tensor components

DO K=0,M%KBAR
   DO J=0,M%JBAR
      DO I=0,M%IBAR
         DUDY = M%RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = M%RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = M%RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = M%RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = M%RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = M%RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
         MUX = 0.25_EB*(M%MU(I,J+1,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I,J+1,K+1))
         MUY = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I+1,J,K+1))
         MUZ = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J+1,K)+M%MU(I+1,J+1,K))
         TXY(I,J,K) = MUZ*(DVDX + DUDY)
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
   ENDDO
ENDDO

! Compute gravity components

IF (.NOT.SPATIAL_GRAVITY_VARIATION) THEN
   GX(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GX)*GVEC(1)
   GY(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GY)*GVEC(2)
   GZ(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GZ)*GVEC(3)
ELSE
   DO I=0,M%IBAR
      GX(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GX)*GVEC(1)
      GY(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GY)*GVEC(2)
      GZ(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GZ)*GVEC(3)
   ENDDO
ENDIF

! Compute x-direction flux term FVX


DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I,J-1,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEYP  = M%CELL(IC)%EDGE_INDEX(8)
         IEYM  = M%CELL(IC)%EDGE_INDEX(6)
         IEZP  = M%CELL(IC)%EDGE_INDEX(12)
         IEZM  = M%CELL(IC)%EDGE_INDEX(10)
         IF (M%EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
            OMYP = M%EDGE(IEYP)%OMEGA(-1)
            TXZP = M%EDGE(IEYP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
            OMYM = M%EDGE(IEYM)%OMEGA( 1)
            TXZM = M%EDGE(IEYM)%TAU( 1)
         ENDIF
         IF (M%EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
            OMZP = M%EDGE(IEZP)%OMEGA(-2)
            TXYP = M%EDGE(IEZP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
            OMZM = M%EDGE(IEZM)%OMEGA( 2)
            TXYM = M%EDGE(IEZM)%TAU( 2)
         ENDIF
         WOMY  = WP*OMYP + WM*OMYM
         VOMZ  = VP*OMZP + VM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*M%RDY(J)
         DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*M%RDZ(K)
         TXXP  = M%MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*M%RDY(J)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*M%RDZ(K)
         TXXM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
         DTXXDX= M%RDXN(I)*(TXXP-TXXM)
         DTXYDY= M%RDY(J) *(TXYP-TXYM)
         DTXZDZ= M%RDZ(K) *(TXZP-TXZM)
         VTRM  = DTXXDX + DTXYDY + DTXZDZ
         M%FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*M%RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO

! Compute y-direction flux term FVY

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I-1,J,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEXP  = M%CELL(IC)%EDGE_INDEX(4)
         IEXM  = M%CELL(IC)%EDGE_INDEX(2)
         IEZP  = M%CELL(IC)%EDGE_INDEX(12)
         IEZM  = M%CELL(IC)%EDGE_INDEX(11)
         IF (M%EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
            OMXP = M%EDGE(IEXP)%OMEGA(-2)
            TYZP = M%EDGE(IEXP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
            OMXM = M%EDGE(IEXM)%OMEGA( 2)
            TYZM = M%EDGE(IEXM)%TAU( 2)
         ENDIF
         IF (M%EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
            OMZP = M%EDGE(IEZP)%OMEGA(-1)
            TXYP = M%EDGE(IEZP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
            OMZM = M%EDGE(IEZM)%OMEGA( 1)
            TXYM = M%EDGE(IEZM)%TAU( 1)
         ENDIF
         WOMX  = WP*OMXP + WM*OMXM
         UOMZ  = UP*OMZP + UM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*M%RDX(I)
         DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*M%RDZ(K)
         TYYP  = M%MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*M%RDX(I)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*M%RDZ(K)
         TYYM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
         DTXYDX= M%RDX(I) *(TXYP-TXYM)
         DTYYDY= M%RDYN(J)*(TYYP-TYYM)
         DTYZDZ= M%RDZ(K) *(TYZP-TYZM)
         VTRM  = DTXYDX + DTYYDY + DTYZDZ
         M%FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*M%RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO

! Compute z-direction flux term FVZ

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I-1,J,K)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J-1,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEXP  = M%CELL(IC)%EDGE_INDEX(4)
         IEXM  = M%CELL(IC)%EDGE_INDEX(3)
         IEYP  = M%CELL(IC)%EDGE_INDEX(8)
         IEYM  = M%CELL(IC)%EDGE_INDEX(7)
         IF (M%EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
            OMXP = M%EDGE(IEXP)%OMEGA(-1)
            TYZP = M%EDGE(IEXP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
            OMXM = M%EDGE(IEXM)%OMEGA( 1)
            TYZM = M%EDGE(IEXM)%TAU( 1)
         ENDIF
         IF (M%EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
            OMYP = M%EDGE(IEYP)%OMEGA(-2)
            TXZP = M%EDGE(IEYP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
            OMYM = M%EDGE(IEYM)%OMEGA( 2)
            TXZM = M%EDGE(IEYM)%TAU( 2)
         ENDIF
         UOMY  = UP*OMYP + UM*OMYM
         VOMX  = VP*OMXP + VM*OMXM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*M%RDX(I)
         DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*M%RDY(J)
         TZZP  = M%MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*M%RDX(I)
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*M%RDY(J)
         TZZM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
         DTXZDX= M%RDX(I) *(TXZP-TXZM)
         DTYZDY= M%RDY(J) *(TYZP-TYZM)
         DTZZDZ= M%RDZN(K)*(TZZP-TZZM)
         VTRM  = DTXZDX + DTYZDY + DTZZDZ
         M%FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(M%RHO_0(K)+M%RHO_0(K+1)) - VTRM)
      ENDDO
   ENDDO
ENDDO


! Additional force terms

IF (OPEN_WIND_BOUNDARY) CALL COMPUTE_WIND_COMPONENTS(T,NM)

IF (ANY(ABS(FVEC)>TWENTY_EPSILON_EB) .OR. CTRL_DIRECT_FORCE) CALL DIRECT_FORCE        ! Direct force
IF (ANY(ABS(OVEC)>TWENTY_EPSILON_EB))                        CALL CORIOLIS_FORCE      ! Coriolis force
IF (PATCH_VELOCITY)                                       CALL PATCH_VELOCITY_FLUX ! Specified patch velocity
IF (PERIODIC_TEST==7)                                     CALL MMS_VELOCITY_FLUX   ! Source term in manufactured solution
IF (PERIODIC_TEST==21 .OR. PERIODIC_TEST==22 .OR. PERIODIC_TEST==23) CALL ROTATED_CUBE_VELOCITY_FLUX(NM,T)

CONTAINS

SUBROUTINE DIRECT_FORCE()

USE CONTROL_VARIABLES, ONLY: CONTROL,N_CTRL

REAL(EB) :: TIME_RAMP_FACTOR,SIN_THETA,COS_THETA,THETA
INTEGER :: N

! CTRL_DIRECT_FORCE overrides FORCE_VECTOR

IF (CTRL_DIRECT_FORCE) THEN
   DO N=1,N_CTRL
      IF (CONTROL(N)%CONTROL_FORCE(1)) FVEC(1) = FVEC(1) - CONTROL(N)%INSTANT_VALUE
      IF (CONTROL(N)%CONTROL_FORCE(2)) FVEC(2) = FVEC(2) - CONTROL(N)%INSTANT_VALUE
      IF (CONTROL(N)%CONTROL_FORCE(3)) FVEC(3) = FVEC(3) - CONTROL(N)%INSTANT_VALUE
   ENDDO
ENDIF

IF (I_RAMP_DIRECTION_T/=0) THEN
   THETA = EVALUATE_RAMP(T,I_RAMP_DIRECTION_T)*DEG2RAD
   SIN_THETA = -SIN(THETA)
   COS_THETA = -COS(THETA)
ELSE
   SIN_THETA = 1._EB
   COS_THETA = 1._EB
ENDIF

IF (ABS(FVEC(1))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVX_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVX_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=0,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
            M%FVX(I,J,K) = M%FVX(I,J,K) - RRHO*FVEC(1)*TIME_RAMP_FACTOR*SIN_THETA
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (ABS(FVEC(2))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVY_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVY_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
            M%FVY(I,J,K) = M%FVY(I,J,K) - RRHO*FVEC(2)*TIME_RAMP_FACTOR*COS_THETA
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (ABS(FVEC(3))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVZ_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVZ_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=0,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
            M%FVZ(I,J,K) = M%FVZ(I,J,K) - RRHO*FVEC(3)*TIME_RAMP_FACTOR
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE DIRECT_FORCE


SUBROUTINE CORIOLIS_FORCE()

REAL(EB), POINTER, DIMENSION(:,:,:) :: UP,VP,WP
REAL(EB) :: UBAR,VBAR,WBAR
INTEGER :: IW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! Velocities relative to the p-cell center (same work done in Deardorff eddy viscosity)

UP => M%WORK7
VP => M%WORK8
WP => M%WORK9
UP=0._EB
VP=0._EB
WP=0._EB

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
      ENDDO
   ENDDO
ENDDO

DO IW=1,M%N_EXTERNAL_WALL_CELLS
   WC=>M%WALL(IW)
   BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
   UP(BC%II,BC%JJ,BC%KK) = M%U_GHOST(IW)
   VP(BC%II,BC%JJ,BC%KK) = M%V_GHOST(IW)
   WP(BC%II,BC%JJ,BC%KK) = M%W_GHOST(IW)
ENDDO

! x momentum

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I+1,J,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I+1,J,K))
         M%FVX(I,J,K) = M%FVX(I,J,K) + 2._EB*(OVEC(2)*WBAR-OVEC(3)*VBAR)
      ENDDO
   ENDDO
ENDDO

! y momentum

DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J+1,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I,J+1,K))
         M%FVY(I,J,K) = M%FVY(I,J,K) + 2._EB*(OVEC(3)*UBAR - OVEC(1)*WBAR)
      ENDDO
   ENDDO
ENDDO

! z momentum

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J,K+1))
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I,J,K+1))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + 2._EB*(OVEC(1)*VBAR - OVEC(2)*UBAR)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CORIOLIS_FORCE


SUBROUTINE MMS_VELOCITY_FLUX

! Shunn et al., JCP (2012) prob 3

USE MANUFACTURED_SOLUTIONS, ONLY: VD2D_MMS_U_SRC_3,VD2D_MMS_V_SRC_3

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%FVX(I,J,K) = M%FVX(I,J,K) - VD2D_MMS_U_SRC_3(M%X(I),M%ZC(K),T)
      ENDDO
   ENDDO
ENDDO

DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%FVZ(I,J,K) = M%FVZ(I,J,K) - VD2D_MMS_V_SRC_3(M%XC(I),M%Z(K),T)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE MMS_VELOCITY_FLUX


!> \brief Compute the velocity flux at a user-specified patch
!> \details The user may specify a polynomial profile using the PROP and DEVC lines. This routine
!> specifies the source term in the momentum equation to drive the local velocity toward
!> this user-specified value, in much the same way as the immersed boundary method
!> (see CC_VELOCITY_FLUX).

SUBROUTINE PATCH_VELOCITY_FLUX

USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE,PROPERTY_TYPE,N_DEVC,DEVICE,PROPERTY
USE TRAN, ONLY: GINV
TYPE(DEVICE_TYPE), POINTER :: DV
TYPE(PROPERTY_TYPE), POINTER :: PY
INTEGER :: N,I1,I2,J1,J2,K1,K2
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
REAL(EB) :: VELP,DX0,DY0,DZ0

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   HP => M%HS
ELSE
   HP => M%H
ENDIF

DEVC_LOOP: DO N=1,N_DEVC

   DV=>DEVICE(N)
   IF (DV%QUANTITY(1)/='VELOCITY PATCH') CYCLE DEVC_LOOP
   IF (DV%PROP_INDEX<1)               CYCLE DEVC_LOOP
   IF (.NOT.DEVICE(DV%DEVC_INDEX(1))%CURRENT_STATE) CYCLE DEVC_LOOP

   IF (DV%X1 > M%XF .OR. DV%X2 < M%XS .OR. &
       DV%Y1 > M%YF .OR. DV%Y2 < M%YS .OR. &
       DV%Z1 > M%ZF .OR. DV%Z2 < M%ZS) CYCLE DEVC_LOOP

   PY=>PROPERTY(DV%PROP_INDEX)

   I_VEL_SELECT: SELECT CASE(PY%I_VEL)

      CASE(1) I_VEL_SELECT

         I1 = MAX(0,      NINT( GINV(DV%X1-M%XS,1,NM)*M%RDXI   )-1)
         I2 = MIN(M%IBAR, NINT( GINV(DV%X2-M%XS,1,NM)*M%RDXI   )+1)
         J1 = MAX(0,      NINT( GINV(DV%Y1-M%YS,2,NM)*M%RDETA  )-1)
         J2 = MIN(M%JBAR, NINT( GINV(DV%Y2-M%YS,2,NM)*M%RDETA  )+1)
         K1 = MAX(0,      NINT( GINV(DV%Z1-M%ZS,3,NM)*M%RDZETA )-1)
         K2 = MIN(M%KBAR, NINT( GINV(DV%Z2-M%ZS,3,NM)*M%RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = M%CELL_INDEX(I,J,K)
                  IC2 = M%CELL_INDEX(I+1,J,K)
                  IF (M%CELL(IC1)%SOLID .OR. M%CELL(IC2)%SOLID) CYCLE

                  IF ( M%X(I)<DV%X1 .OR.  M%X(I)>DV%X2) CYCLE ! Inefficient but simple
                  IF (M%YC(J)<DV%Y1 .OR. M%YC(J)>DV%Y2) CYCLE
                  IF (M%ZC(K)<DV%Z1 .OR. M%ZC(K)>DV%Z2) CYCLE

                  DX0 =  M%X(I)-DV%X
                  DY0 = M%YC(J)-DV%Y
                  DZ0 = M%ZC(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  M%FVX(I,J,K) = -M%RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - (VELP-UU(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

      CASE(2) I_VEL_SELECT

         I1 = MAX(0,      NINT( GINV(DV%X1-M%XS,1,NM)*M%RDXI   )-1)
         I2 = MIN(M%IBAR, NINT( GINV(DV%X2-M%XS,1,NM)*M%RDXI   )+1)
         J1 = MAX(0,      NINT( GINV(DV%Y1-M%YS,2,NM)*M%RDETA  )-1)
         J2 = MIN(M%JBAR, NINT( GINV(DV%Y2-M%YS,2,NM)*M%RDETA  )+1)
         K1 = MAX(0,      NINT( GINV(DV%Z1-M%ZS,3,NM)*M%RDZETA )-1)
         K2 = MIN(M%KBAR, NINT( GINV(DV%Z2-M%ZS,3,NM)*M%RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = M%CELL_INDEX(I,J,K)
                  IC2 = M%CELL_INDEX(I,J+1,K)

                  IF (M%CELL(IC1)%SOLID .OR. M%CELL(IC2)%SOLID) CYCLE

                  IF (M%XC(I)<DV%X1 .OR. M%XC(I)>DV%X2) CYCLE
                  IF ( M%Y(J)<DV%Y1 .OR.  M%Y(J)>DV%Y2) CYCLE
                  IF (M%ZC(K)<DV%Z1 .OR. M%ZC(K)>DV%Z2) CYCLE

                  DX0 = M%XC(I)-DV%X
                  DY0 =  M%Y(J)-DV%Y
                  DZ0 = M%ZC(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  M%FVY(I,J,K) = -M%RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - (VELP-VV(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

      CASE(3) I_VEL_SELECT

         I1 = MAX(0,      NINT( GINV(DV%X1-M%XS,1,NM)*M%RDXI   )-1)
         I2 = MIN(M%IBAR, NINT( GINV(DV%X2-M%XS,1,NM)*M%RDXI   )+1)
         J1 = MAX(0,      NINT( GINV(DV%Y1-M%YS,2,NM)*M%RDETA  )-1)
         J2 = MIN(M%JBAR, NINT( GINV(DV%Y2-M%YS,2,NM)*M%RDETA  )+1)
         K1 = MAX(0,      NINT( GINV(DV%Z1-M%ZS,3,NM)*M%RDZETA )-1)
         K2 = MIN(M%KBAR, NINT( GINV(DV%Z2-M%ZS,3,NM)*M%RDZETA )+1)

         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IC1 = M%CELL_INDEX(I,J,K)
                  IC2 = M%CELL_INDEX(I,J,K+1)
                  IF (M%CELL(IC1)%SOLID .OR. M%CELL(IC2)%SOLID) CYCLE

                  IF (M%XC(I)<DV%X1 .OR. M%XC(I)>DV%X2) CYCLE
                  IF (M%YC(J)<DV%Y1 .OR. M%YC(J)>DV%Y2) CYCLE
                  IF ( M%Z(K)<DV%Z1 .OR.  M%Z(K)>DV%Z2) CYCLE

                  DX0 = M%XC(I)-DV%X
                  DY0 = M%YC(J)-DV%Y
                  DZ0 =  M%Z(K)-DV%Z
                  VELP = PY%P0 + DX0*PY%PX(1) + 0.5_EB*(DX0*DX0*PY%PXX(1,1)+DX0*DY0*PY%PXX(1,2)+DX0*DZ0*PY%PXX(1,3)) &
                               + DY0*PY%PX(2) + 0.5_EB*(DY0*DX0*PY%PXX(2,1)+DY0*DY0*PY%PXX(2,2)+DY0*DZ0*PY%PXX(2,3)) &
                               + DZ0*PY%PX(3) + 0.5_EB*(DZ0*DX0*PY%PXX(3,1)+DZ0*DY0*PY%PXX(3,2)+DZ0*DZ0*PY%PXX(3,3))

                  M%FVZ(I,J,K) = -M%RDZN(K)*(HP(I,J,K)-HP(I,J,K+1)) - (VELP-WW(I,J,K))/DT
               ENDDO
            ENDDO
         ENDDO

   END SELECT I_VEL_SELECT

ENDDO DEVC_LOOP

END SUBROUTINE PATCH_VELOCITY_FLUX

END SUBROUTINE VELOCITY_FLUX_KERNEL


!> \brief Block-decomposed velocity flux: computes vorticity, FVX, FVY, FVZ for K-range [K1, K2].
!> \details Only includes main cell loops and DIRECT_FORCE (without CTRL_DIRECT_FORCE).
!> Features requiring mesh-level execution (Coriolis, patch velocity, open wind, periodic tests)
!> are excluded — caller must fall back to VELOCITY_FLUX_KERNEL when these are active.
!> \param M Mesh data structure
!> \param T Current time (s)
!> \param DT Time step (s)
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables
!> \param K1 Start of K cell range (1-based inclusive)
!> \param K2 End of K cell range (1-based inclusive)

RECURSIVE SUBROUTINE VELOCITY_FLUX_BLOCK_KERNEL(M,T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES,K1,K2)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM,K1,K2
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
            RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ
INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC
INTEGER :: K1_VORT,K1_FVZ,K2_FVZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ,UU,VV,WW,RHOP,DP
REAL(EB) :: GX(0:M%IBAR),GY(0:M%IBAR),GZ(0:M%IBAR)

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
   DP => M%DS
   RHOP => M%RHOS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
   DP => M%D
   RHOP => M%RHO
ENDIF

TXY => M%WORK1
TXZ => M%WORK2
TYZ => M%WORK3
OMX => M%WORK4
OMY => M%WORK5
OMZ => M%WORK6

! K ranges: vorticity extended 1 cell below for FVX/FVY dependency at K-1
K1_VORT = MAX(0, K1-1)

! FVZ staggered K range (same pattern as W in velocity predictor/corrector)
K1_FVZ = K1 - 1
K2_FVZ = K2 - 1
IF (K2==M%KBAR) K2_FVZ = M%KBAR

! Compute vorticity and stress tensor components

DO K=K1_VORT,K2
   DO J=0,M%JBAR
      DO I=0,M%IBAR
         DUDY = M%RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = M%RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = M%RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = M%RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = M%RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = M%RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
         MUX = 0.25_EB*(M%MU(I,J+1,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I,J+1,K+1))
         MUY = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I+1,J,K+1))
         MUZ = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J+1,K)+M%MU(I+1,J+1,K))
         TXY(I,J,K) = MUZ*(DVDX + DUDY)
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
   ENDDO
ENDDO

! Compute gravity components (1D, cheap — computed redundantly per block)

IF (.NOT.SPATIAL_GRAVITY_VARIATION) THEN
   GX(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GX)*GVEC(1)
   GY(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GY)*GVEC(2)
   GZ(0:M%IBAR) = EVALUATE_RAMP(T,I_RAMP_GZ)*GVEC(3)
ELSE
   DO I=0,M%IBAR
      GX(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GX)*GVEC(1)
      GY(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GY)*GVEC(2)
      GZ(I) = EVALUATE_RAMP(M%X(I),I_RAMP_GZ)*GVEC(3)
   ENDDO
ENDIF

! Compute x-direction flux term FVX

DO K=K1,K2
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I,J-1,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEYP  = M%CELL(IC)%EDGE_INDEX(8)
         IEYM  = M%CELL(IC)%EDGE_INDEX(6)
         IEZP  = M%CELL(IC)%EDGE_INDEX(12)
         IEZM  = M%CELL(IC)%EDGE_INDEX(10)
         IF (M%EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
            OMYP = M%EDGE(IEYP)%OMEGA(-1)
            TXZP = M%EDGE(IEYP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
            OMYM = M%EDGE(IEYM)%OMEGA( 1)
            TXZM = M%EDGE(IEYM)%TAU( 1)
         ENDIF
         IF (M%EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
            OMZP = M%EDGE(IEZP)%OMEGA(-2)
            TXYP = M%EDGE(IEZP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
            OMZM = M%EDGE(IEZM)%OMEGA( 2)
            TXYM = M%EDGE(IEZM)%TAU( 2)
         ENDIF
         WOMY  = WP*OMYP + WM*OMYM
         VOMZ  = VP*OMZP + VM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*M%RDY(J)
         DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*M%RDZ(K)
         TXXP  = M%MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*M%RDY(J)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*M%RDZ(K)
         TXXM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
         DTXXDX= M%RDXN(I)*(TXXP-TXXM)
         DTXYDY= M%RDY(J) *(TXYP-TXYM)
         DTXZDZ= M%RDZ(K) *(TXZP-TXZM)
         VTRM  = DTXXDX + DTXYDY + DTXZDZ
         M%FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*M%RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO

! Compute y-direction flux term FVY

DO K=K1,K2
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I-1,J,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEXP  = M%CELL(IC)%EDGE_INDEX(4)
         IEXM  = M%CELL(IC)%EDGE_INDEX(2)
         IEZP  = M%CELL(IC)%EDGE_INDEX(12)
         IEZM  = M%CELL(IC)%EDGE_INDEX(11)
         IF (M%EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
            OMXP = M%EDGE(IEXP)%OMEGA(-2)
            TYZP = M%EDGE(IEXP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
            OMXM = M%EDGE(IEXM)%OMEGA( 2)
            TYZM = M%EDGE(IEXM)%TAU( 2)
         ENDIF
         IF (M%EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
            OMZP = M%EDGE(IEZP)%OMEGA(-1)
            TXYP = M%EDGE(IEZP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
            OMZM = M%EDGE(IEZM)%OMEGA( 1)
            TXYM = M%EDGE(IEZM)%TAU( 1)
         ENDIF
         WOMX  = WP*OMXP + WM*OMXM
         UOMZ  = UP*OMZP + UM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*M%RDX(I)
         DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*M%RDZ(K)
         TYYP  = M%MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*M%RDX(I)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*M%RDZ(K)
         TYYM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
         DTXYDX= M%RDX(I) *(TXYP-TXYM)
         DTYYDY= M%RDYN(J)*(TYYP-TYYM)
         DTYZDZ= M%RDZ(K) *(TYZP-TYZM)
         VTRM  = DTXYDX + DTYYDY + DTYZDZ
         M%FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*M%RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO

! Compute z-direction flux term FVZ

DO K=K1_FVZ,K2_FVZ
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I-1,J,K)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J-1,K)
         IC    = M%CELL_INDEX(I,J,K)
         IEXP  = M%CELL(IC)%EDGE_INDEX(4)
         IEXM  = M%CELL(IC)%EDGE_INDEX(3)
         IEYP  = M%CELL(IC)%EDGE_INDEX(8)
         IEYM  = M%CELL(IC)%EDGE_INDEX(7)
         IF (M%EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
            OMXP = M%EDGE(IEXP)%OMEGA(-1)
            TYZP = M%EDGE(IEXP)%TAU(-1)
         ENDIF
         IF (M%EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
            OMXM = M%EDGE(IEXM)%OMEGA( 1)
            TYZM = M%EDGE(IEXM)%TAU( 1)
         ENDIF
         IF (M%EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
            OMYP = M%EDGE(IEYP)%OMEGA(-2)
            TXZP = M%EDGE(IEYP)%TAU(-2)
         ENDIF
         IF (M%EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
            OMYM = M%EDGE(IEYM)%OMEGA( 2)
            TXZM = M%EDGE(IEYM)%TAU( 2)
         ENDIF
         UOMY  = UP*OMYP + UM*OMYM
         VOMX  = VP*OMXP + VM*OMXM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*M%RDX(I)
         DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*M%RDY(J)
         TZZP  = M%MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*M%RDX(I)
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*M%RDY(J)
         TZZM  = M%MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
         DTXZDX= M%RDX(I) *(TXZP-TXZM)
         DTYZDY= M%RDY(J) *(TYZP-TYZM)
         DTZZDZ= M%RDZN(K)*(TZZP-TZZM)
         VTRM  = DTXZDX + DTYZDY + DTZZDZ
         M%FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(M%RHO_0(K)+M%RHO_0(K+1)) - VTRM)
      ENDDO
   ENDDO
ENDDO

! DIRECT_FORCE (K-restricted cell loops, excludes CTRL_DIRECT_FORCE)

IF (ANY(ABS(FVEC)>TWENTY_EPSILON_EB)) CALL DIRECT_FORCE_BLOCK

CONTAINS

SUBROUTINE DIRECT_FORCE_BLOCK()

REAL(EB) :: TIME_RAMP_FACTOR,SIN_THETA,COS_THETA,THETA

IF (I_RAMP_DIRECTION_T/=0) THEN
   THETA = EVALUATE_RAMP(T,I_RAMP_DIRECTION_T)*DEG2RAD
   SIN_THETA = -SIN(THETA)
   COS_THETA = -COS(THETA)
ELSE
   SIN_THETA = 1._EB
   COS_THETA = 1._EB
ENDIF

IF (ABS(FVEC(1))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVX_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVX_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=K1,K2
      DO J=1,M%JBAR
         DO I=0,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
            M%FVX(I,J,K) = M%FVX(I,J,K) - RRHO*FVEC(1)*TIME_RAMP_FACTOR*SIN_THETA
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (ABS(FVEC(2))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVY_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVY_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=K1,K2
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
            M%FVY(I,J,K) = M%FVY(I,J,K) - RRHO*FVEC(2)*TIME_RAMP_FACTOR*COS_THETA
         ENDDO
      ENDDO
   ENDDO
ENDIF

IF (ABS(FVEC(3))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVZ_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVZ_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   DO K=K1_FVZ,K2_FVZ
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
            M%FVZ(I,J,K) = M%FVZ(I,J,K) - RRHO*FVEC(3)*TIME_RAMP_FACTOR
         ENDDO
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE DIRECT_FORCE_BLOCK

END SUBROUTINE VELOCITY_FLUX_BLOCK_KERNEL


!> \brief Compute the turbulent viscosity.
!> \param M Mesh data structure
!> \param NM Mesh number (needed for external calls)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables

RECURSIVE SUBROUTINE COMPUTE_VISCOSITY_KERNEL(M,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY,GET_POTENTIAL_TEMPERATURE,GET_CONDUCTIVITY,GET_SPECIFIC_HEAT
USE TURB_KERNELS, ONLY: WALE_VISCOSITY,FILL_EDGES_KERNEL,TEST_FILTER_KERNEL,VARDEN_DYNSMAG_KERNEL
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE CC_VELOCITY_KERNELS, ONLY : CC_COMPUTE_KRES,CC_COMPUTE_VISCOSITY,CUTFACE_VELOCITIES

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
REAL(EB) :: NU_EDDY,DELTA,KSGS,U2,V2,W2,AA,A_IJ(3,3),BB,B_IJ(3,3),&
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,VDF,WGT
REAL(EB), PARAMETER :: RAPLUS=1._EB/26._EB
INTEGER :: I,J,K,IIG,JJG,KKG,II,JJ,KK,IW,IOR,IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UP,VP,WP, &
                                       UP_HAT,VP_HAT,WP_HAT, &
                                       UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_COUNTER
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(SURFACE_TYPE), POINTER :: SF

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   RHOP => M%RHOS
   UU   => M%US
   VV   => M%VS
   WW   => M%WS
   ZZP  => M%ZZS
ELSE
   RHOP => M%RHO
   UU   => M%U
   VV   => M%V
   WW   => M%W
   ZZP  => M%ZZ
ENDIF

! Compute viscosity for DNS using primitive species

IF (SIM_MODE==SVLES_MODE) THEN

   M%MU_DNS = MU_AIR_0

ELSE

   ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,M%MU_DNS(I,J,K),M%TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(ZZ_GET)

ENDIF

IF (CC_IBM) CALL CUTFACE_VELOCITIES(M,UU,VV,WW, &
   CUTFACES=.TRUE.)

CALL COMPUTE_STRAIN_RATE

SELECT_TURB: SELECT CASE (TURB_MODEL)

   CASE (NO_TURB_MODEL)

      M%MU = M%MU_DNS

   CASE (CONSMAG,DYNSMAG) SELECT_TURB ! Smagorinsky (1963) eddy viscosity

      IF (PREDICTOR .AND. TURB_MODEL==DYNSMAG) CALL VARDEN_DYNSMAG_KERNEL(M) ! dynamic procedure, Moin et al. (1991)

      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*M%CSD2(I,J,K)*M%STRAIN_RATE(I,J,K)
            ENDDO
         ENDDO
      ENDDO

   CASE (DEARDORFF) SELECT_TURB ! Deardorff (1980) eddy viscosity model (current default)

      ! Velocities relative to the p-cell center

      UP => M%WORK1
      VP => M%WORK2
      WP => M%WORK3
      UP=0._EB
      VP=0._EB
      WP=0._EB


      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
               VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
               WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
            ENDDO
         ENDDO
      ENDDO

      ! fill mesh boundary ghost cells

      DO IW=1,M%N_EXTERNAL_WALL_CELLS
         WC=>M%WALL(IW)
         BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
         SELECT CASE(WC%BOUNDARY_TYPE)
            CASE(INTERPOLATED_BOUNDARY)
               UP(BC%II,BC%JJ,BC%KK) = M%U_GHOST(IW)
               VP(BC%II,BC%JJ,BC%KK) = M%V_GHOST(IW)
               WP(BC%II,BC%JJ,BC%KK) = M%W_GHOST(IW)
            CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)
               UP(BC%II,BC%JJ,BC%KK) = UP(BC%IIG,BC%JJG,BC%KKG)
               VP(BC%II,BC%JJ,BC%KK) = VP(BC%IIG,BC%JJG,BC%KKG)
               WP(BC%II,BC%JJ,BC%KK) = WP(BC%IIG,BC%JJG,BC%KKG)
         END SELECT
      ENDDO


      ! fill edge and corner ghost cells

      CALL FILL_EDGES_KERNEL(M,UP)
      CALL FILL_EDGES_KERNEL(M,VP)
      CALL FILL_EDGES_KERNEL(M,WP)

      UP_HAT => M%WORK4
      VP_HAT => M%WORK5
      WP_HAT => M%WORK6
      UP_HAT=0._EB
      VP_HAT=0._EB
      WP_HAT=0._EB

      CALL TEST_FILTER_KERNEL(M,UP_HAT,UP)
      CALL TEST_FILTER_KERNEL(M,VP_HAT,VP)
      CALL TEST_FILTER_KERNEL(M,WP_HAT,WP)

      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DELTA = M%LES_FILTER_WIDTH(I,J,K)
               KSGS = 0.5_EB*( (UP(I,J,K)-UP_HAT(I,J,K))**2 + (VP(I,J,K)-VP_HAT(I,J,K))**2 + (WP(I,J,K)-WP_HAT(I,J,K))**2 )
               NU_EDDY = C_DEARDORFF*DELTA*SQRT(KSGS)
               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
            ENDDO
         ENDDO
      ENDDO

   CASE (VREMAN) SELECT_TURB ! Vreman (2004) eddy viscosity model (experimental)

      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))

               ! Vreman, Eq. (6)
               A_IJ(1,1)=DUDX; A_IJ(2,1)=DUDY; A_IJ(3,1)=DUDZ
               A_IJ(1,2)=DVDX; A_IJ(2,2)=DVDY; A_IJ(3,2)=DVDZ
               A_IJ(1,3)=DWDX; A_IJ(2,3)=DWDY; A_IJ(3,3)=DWDZ

               AA=0._EB
               DO JJ=1,3
                  DO II=1,3
                     AA = AA + A_IJ(II,JJ)*A_IJ(II,JJ)
                  ENDDO
               ENDDO

               ! Vreman, Eq. (7)
               B_IJ(1,1)=(M%DX(I)*A_IJ(1,1))**2 + (M%DY(J)*A_IJ(2,1))**2 + (M%DZ(K)*A_IJ(3,1))**2
               B_IJ(2,2)=(M%DX(I)*A_IJ(1,2))**2 + (M%DY(J)*A_IJ(2,2))**2 + (M%DZ(K)*A_IJ(3,2))**2
               B_IJ(3,3)=(M%DX(I)*A_IJ(1,3))**2 + (M%DY(J)*A_IJ(2,3))**2 + (M%DZ(K)*A_IJ(3,3))**2

               B_IJ(1,2)=M%DX(I)**2*A_IJ(1,1)*A_IJ(1,2) + M%DY(J)**2*A_IJ(2,1)*A_IJ(2,2) + M%DZ(K)**2*A_IJ(3,1)*A_IJ(3,2)
               B_IJ(1,3)=M%DX(I)**2*A_IJ(1,1)*A_IJ(1,3) + M%DY(J)**2*A_IJ(2,1)*A_IJ(2,3) + M%DZ(K)**2*A_IJ(3,1)*A_IJ(3,3)
               B_IJ(2,3)=M%DX(I)**2*A_IJ(1,2)*A_IJ(1,3) + M%DY(J)**2*A_IJ(2,2)*A_IJ(2,3) + M%DZ(K)**2*A_IJ(3,2)*A_IJ(3,3)

               BB = B_IJ(1,1)*B_IJ(2,2) - B_IJ(1,2)**2 &
                  + B_IJ(1,1)*B_IJ(3,3) - B_IJ(1,3)**2 &
                  + B_IJ(2,2)*B_IJ(3,3) - B_IJ(2,3)**2    ! Vreman, Eq. (8)

               IF (ABS(AA)>TWENTY_EPSILON_EB .AND. BB>TWENTY_EPSILON_EB) THEN
                  NU_EDDY = C_VREMAN*SQRT(BB/AA)  ! Vreman, Eq. (5)
               ELSE
                  NU_EDDY=0._EB
               ENDIF

               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY

            ENDDO
         ENDDO
      ENDDO

   CASE (WALE) SELECT_TURB

      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DELTA = M%LES_FILTER_WIDTH(I,J,K)
               ! compute velocity gradient tensor
               DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
               A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
               A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ

               CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)

               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
            ENDDO
         ENDDO
      ENDDO

END SELECT SELECT_TURB

! Compute resolved kinetic energy per unit mass

DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         U2 = 0.25_EB*(UU(I-1,J,K)+UU(I,J,K))**2
         V2 = 0.25_EB*(VV(I,J-1,K)+VV(I,J,K))**2
         W2 = 0.25_EB*(WW(I,J,K-1)+WW(I,J,K))**2
         M%KRES(I,J,K) = 0.5_EB*(U2+V2+W2)
      ENDDO
   ENDDO
ENDDO

IF (CC_IBM) CALL CC_COMPUTE_KRES(M, &
   APPLY_TO_ESTIMATED_VARIABLES)

! Mirror viscosity into solids and exterior boundary cells

CELL_COUNTER => M%IWORK1 ; CELL_COUNTER = 0

WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS

   WC=>M%WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
   B1=>M%BOUNDARY_PROP1(WC%B1_INDEX)
   B2=>M%BOUNDARY_PROP2(WC%B2_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IC  = M%CELL_INDEX(II,JJ,KK)
   IOR = BC%IOR
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   SF=>SURFACE(WC%SURF_INDEX)

   IF (M%CELL(IC)%SOLID .OR. M%CELL(IC)%EXTERIOR) M%KRES(II,JJ,KK) = M%KRES(IIG,JJG,KKG)

   SELECT CASE(WC%BOUNDARY_TYPE)

      CASE(SOLID_BOUNDARY)

         IF (SIM_MODE/=DNS_MODE) THEN
            DELTA = M%LES_FILTER_WIDTH(IIG,JJG,KKG)
            SELECT CASE(SF%NEAR_WALL_TURB_MODEL)
               CASE DEFAULT
                  NU_EDDY = 0._EB
               CASE(CONSTANT_EDDY_VISCOSITY)
                  NU_EDDY = SF%NEAR_WALL_EDDY_VISCOSITY
               CASE(CONSMAG) ! Constant Smagorinsky with Van Driest damping
                  VDF = 1._EB-EXP(-B2%Y_PLUS*RAPLUS)
                  NU_EDDY = (VDF*C_SMAGORINSKY*DELTA)**2*M%STRAIN_RATE(IIG,JJG,KKG)
               CASE(WALE)
                  ! compute velocity gradient tensor
                  DUDX = M%RDX(IIG)*(UU(IIG,JJG,KKG)-UU(IIG-1,JJG,KKG))
                  DVDY = M%RDY(JJG)*(VV(IIG,JJG,KKG)-VV(IIG,JJG-1,KKG))
                  DWDZ = M%RDZ(KKG)*(WW(IIG,JJG,KKG)-WW(IIG,JJG,KKG-1))
                  DUDY = 0.25_EB*M%RDY(JJG)*(UU(IIG,JJG+1,KKG)-UU(IIG,JJG-1,KKG)+UU(IIG-1,JJG+1,KKG)-UU(IIG-1,JJG-1,KKG))
                  DUDZ = 0.25_EB*M%RDZ(KKG)*(UU(IIG,JJG,KKG+1)-UU(IIG,JJG,KKG-1)+UU(IIG-1,JJG,KKG+1)-UU(IIG-1,JJG,KKG-1))
                  DVDX = 0.25_EB*M%RDX(IIG)*(VV(IIG+1,JJG,KKG)-VV(IIG-1,JJG,KKG)+VV(IIG+1,JJG-1,KKG)-VV(IIG-1,JJG-1,KKG))
                  DVDZ = 0.25_EB*M%RDZ(KKG)*(VV(IIG,JJG,KKG+1)-VV(IIG,JJG,KKG-1)+VV(IIG,JJG-1,KKG+1)-VV(IIG,JJG-1,KKG-1))
                  DWDX = 0.25_EB*M%RDX(IIG)*(WW(IIG+1,JJG,KKG)-WW(IIG-1,JJG,KKG)+WW(IIG+1,JJG,KKG-1)-WW(IIG-1,JJG,KKG-1))
                  DWDY = 0.25_EB*M%RDY(JJG)*(WW(IIG,JJG+1,KKG)-WW(IIG,JJG-1,KKG)+WW(IIG,JJG+1,KKG-1)-WW(IIG,JJG-1,KKG-1))
                  A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
                  A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
                  A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ
                  CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)
            END SELECT
            IF (CELL_COUNTER(IIG,JJG,KKG)==0) M%MU(IIG,JJG,KKG) = 0._EB
            CELL_COUNTER(IIG,JJG,KKG) = CELL_COUNTER(IIG,JJG,KKG) + 1
            WGT = 1._EB/REAL(CELL_COUNTER(IIG,JJG,KKG),EB)
            M%MU(IIG,JJG,KKG) = (1._EB-WGT)*M%MU(IIG,JJG,KKG) + WGT*(M%MU_DNS(IIG,JJG,KKG) + RHOP(IIG,JJG,KKG)*NU_EDDY)
         ELSE
            M%MU(IIG,JJG,KKG) = M%MU_DNS(IIG,JJG,KKG)
         ENDIF

         IF (M%CELL(M%CELL_INDEX(II,JJ,KK))%SOLID) M%MU(II,JJ,KK) = M%MU(IIG,JJG,KKG)

      CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)

         M%MU(II,JJ,KK) = M%MU(IIG,JJG,KKG)

   END SELECT

ENDDO WALL_LOOP

IF(CC_IBM) THEN
   CALL CC_COMPUTE_VISCOSITY(M,0._EB)
   CALL CUTFACE_VELOCITIES(M,UU,VV,WW, &
      CUTFACES=.FALSE.)
ENDIF

M%MU(   0,0:M%JBP1,   0) = M%MU(   1,0:M%JBP1,1)
M%MU(M%IBP1,0:M%JBP1,   0) = M%MU(M%IBAR,0:M%JBP1,1)
M%MU(M%IBP1,0:M%JBP1,M%KBP1) = M%MU(M%IBAR,0:M%JBP1,M%KBAR)
M%MU(   0,0:M%JBP1,M%KBP1) = M%MU(   1,0:M%JBP1,M%KBAR)
M%MU(0:M%IBP1,   0,   0) = M%MU(0:M%IBP1,   1,1)
M%MU(0:M%IBP1,M%JBP1,0)    = M%MU(0:M%IBP1,M%JBAR,1)
M%MU(0:M%IBP1,M%JBP1,M%KBP1) = M%MU(0:M%IBP1,M%JBAR,M%KBAR)
M%MU(0:M%IBP1,0,M%KBP1)    = M%MU(0:M%IBP1,   1,M%KBAR)
M%MU(0,   0,0:M%KBP1)    = M%MU(   1,   1,0:M%KBP1)
M%MU(M%IBP1,0,0:M%KBP1)    = M%MU(M%IBAR,   1,0:M%KBP1)
M%MU(M%IBP1,M%JBP1,0:M%KBP1) = M%MU(M%IBAR,M%JBAR,0:M%KBP1)
M%MU(0,M%JBP1,0:M%KBP1)    = M%MU(   1,M%JBAR,0:M%KBP1)

M%KRES(   0,0:M%JBP1,   0) = M%KRES(   1,0:M%JBP1,1)
M%KRES(M%IBP1,0:M%JBP1,   0) = M%KRES(M%IBAR,0:M%JBP1,1)
M%KRES(M%IBP1,0:M%JBP1,M%KBP1) = M%KRES(M%IBAR,0:M%JBP1,M%KBAR)
M%KRES(   0,0:M%JBP1,M%KBP1) = M%KRES(   1,0:M%JBP1,M%KBAR)
M%KRES(0:M%IBP1,   0,   0) = M%KRES(0:M%IBP1,   1,1)
M%KRES(0:M%IBP1,M%JBP1,0)    = M%KRES(0:M%IBP1,M%JBAR,1)
M%KRES(0:M%IBP1,M%JBP1,M%KBP1) = M%KRES(0:M%IBP1,M%JBAR,M%KBAR)
M%KRES(0:M%IBP1,0,M%KBP1)    = M%KRES(0:M%IBP1,   1,M%KBAR)
M%KRES(0,   0,0:M%KBP1)    = M%KRES(   1,   1,0:M%KBP1)
M%KRES(M%IBP1,0,0:M%KBP1)    = M%KRES(M%IBAR,   1,0:M%KBP1)
M%KRES(M%IBP1,M%JBP1,0:M%KBP1) = M%KRES(M%IBAR,M%JBAR,0:M%KBP1)
M%KRES(0,M%JBP1,0:M%KBP1)    = M%KRES(   1,M%JBAR,0:M%KBP1)

CONTAINS

SUBROUTINE COMPUTE_STRAIN_RATE

REAL(EB) :: S11,S22,S33,S12,S13,S23,ONTHDIV
INTEGER :: SURF_INDEX

SELECT CASE (TURB_MODEL)
   CASE DEFAULT
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
               S11 = DUDX - ONTHDIV
               S22 = DVDY - ONTHDIV
               S33 = DWDZ - ONTHDIV
               S12 = 0.5_EB*(DUDY+DVDX)
               S13 = 0.5_EB*(DUDZ+DWDX)
               S23 = 0.5_EB*(DVDZ+DWDY)
               M%STRAIN_RATE(I,J,K) = SQRT(2._EB*(S11**2 + S22**2 + S33**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
            ENDDO
         ENDDO
      ENDDO
   CASE (DEARDORFF)
      ! Here we omit the 3D loop, we only need the wall cell values of STRAIN_RATE
END SELECT

WALL_LOOP_SR: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
   WC=>M%WALL(IW)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP_SR

   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   SURF_INDEX = WC%SURF_INDEX
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   IOR = BC%IOR

   ! Handle the case where OBST lives on an external boundary
   IF (IW>M%N_EXTERNAL_WALL_CELLS) THEN
      SELECT CASE(IOR)
         CASE( 1); IF (IIG>M%IBAR) CYCLE WALL_LOOP_SR
         CASE(-1); IF (IIG<1)    CYCLE WALL_LOOP_SR
         CASE( 2); IF (JJG>M%JBAR) CYCLE WALL_LOOP_SR
         CASE(-2); IF (JJG<1)    CYCLE WALL_LOOP_SR
         CASE( 3); IF (KKG>M%KBAR) CYCLE WALL_LOOP_SR
         CASE(-3); IF (KKG<1)    CYCLE WALL_LOOP_SR
      END SELECT
   ENDIF

   DUDX = M%RDX(IIG)*(UU(IIG,JJG,KKG)-UU(IIG-1,JJG,KKG))
   DVDY = M%RDY(JJG)*(VV(IIG,JJG,KKG)-VV(IIG,JJG-1,KKG))
   DWDZ = M%RDZ(KKG)*(WW(IIG,JJG,KKG)-WW(IIG,JJG,KKG-1))
   ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
   S11 = DUDX - ONTHDIV
   S22 = DVDY - ONTHDIV
   S33 = DWDZ - ONTHDIV

   DUDY = 0.25_EB*M%RDY(JJG)*(UU(IIG,JJG+1,KKG)-UU(IIG,JJG-1,KKG)+UU(IIG-1,JJG+1,KKG)-UU(IIG-1,JJG-1,KKG))
   DUDZ = 0.25_EB*M%RDZ(KKG)*(UU(IIG,JJG,KKG+1)-UU(IIG,JJG,KKG-1)+UU(IIG-1,JJG,KKG+1)-UU(IIG-1,JJG,KKG-1))
   DVDX = 0.25_EB*M%RDX(IIG)*(VV(IIG+1,JJG,KKG)-VV(IIG-1,JJG,KKG)+VV(IIG+1,JJG-1,KKG)-VV(IIG-1,JJG-1,KKG))
   DVDZ = 0.25_EB*M%RDZ(KKG)*(VV(IIG,JJG,KKG+1)-VV(IIG,JJG,KKG-1)+VV(IIG,JJG-1,KKG+1)-VV(IIG,JJG-1,KKG-1))
   DWDX = 0.25_EB*M%RDX(IIG)*(WW(IIG+1,JJG,KKG)-WW(IIG-1,JJG,KKG)+WW(IIG+1,JJG,KKG-1)-WW(IIG-1,JJG,KKG-1))
   DWDY = 0.25_EB*M%RDY(JJG)*(WW(IIG,JJG+1,KKG)-WW(IIG,JJG-1,KKG)+WW(IIG,JJG+1,KKG-1)-WW(IIG,JJG-1,KKG-1))

   S12 = 0.5_EB*(DUDY+DVDX)
   S13 = 0.5_EB*(DUDZ+DWDX)
   S23 = 0.5_EB*(DVDZ+DWDY)

   M%STRAIN_RATE(IIG,JJG,KKG) = SQRT(2._EB*(S11**2 + S22**2 + S33**2 + 2._EB*(S12**2 + S13**2 + S23**2)))
ENDDO WALL_LOOP_SR

END SUBROUTINE COMPUTE_STRAIN_RATE

END SUBROUTINE COMPUTE_VISCOSITY_KERNEL


!> \brief Block-decomposed viscosity: computes MU_DNS, STRAIN_RATE, turb MU, KRES for K-range [K1, K2].
!> \details Only includes main cell loops for non-DEARDORFF/DYNSMAG turb models.
!> Wall loops and corner mirroring are excluded — caller must run POST_BLOCK after all blocks complete.
!> Features requiring mesh-level execution (DEARDORFF, DYNSMAG, CC_IBM) are excluded.
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables
!> \param K1 Start of K cell range (1-based inclusive)
!> \param K2 End of K cell range (1-based inclusive)

RECURSIVE SUBROUTINE COMPUTE_VISCOSITY_BLOCK_KERNEL(M,NM,APPLY_TO_ESTIMATED_VARIABLES,K1,K2)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY
USE TURB_KERNELS, ONLY: WALE_VISCOSITY

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM,K1,K2
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
REAL(EB) :: NU_EDDY,DELTA,U2,V2,W2,AA,A_IJ(3,3),BB,B_IJ(3,3), &
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ, &
            S11,S22,S33,S12,S13,S23,ONTHDIV
INTEGER :: I,J,K,II,JJ
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   RHOP => M%RHOS
   UU   => M%US
   VV   => M%VS
   WW   => M%WS
   ZZP  => M%ZZS
ELSE
   RHOP => M%RHO
   UU   => M%U
   VV   => M%V
   WW   => M%W
   ZZP  => M%ZZ
ENDIF

! Compute viscosity for DNS using primitive species (K=K1:K2)

IF (SIM_MODE==SVLES_MODE) THEN
   DO K=K1,K2
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            M%MU_DNS(I,J,K) = MU_AIR_0
         ENDDO
      ENDDO
   ENDDO
ELSE
   ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
   DO K=K1,K2
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,M%MU_DNS(I,J,K),M%TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   DEALLOCATE(ZZ_GET)
ENDIF

! Compute strain rate (K=K1:K2) — inlined from COMPUTE_STRAIN_RATE
! Reads UU/VV/WW stencil at K-1:K+1 (read-only, safe for parallel blocks)

DO K=K1,K2
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
         DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
         DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
         DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K) &
                +UU(I-1,J+1,K)-UU(I-1,J-1,K))
         DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1) &
                +UU(I-1,J,K+1)-UU(I-1,J,K-1))
         DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K) &
                +VV(I+1,J-1,K)-VV(I-1,J-1,K))
         DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1) &
                +VV(I,J-1,K+1)-VV(I,J-1,K-1))
         DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K) &
                +WW(I+1,J,K-1)-WW(I-1,J,K-1))
         DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K) &
                +WW(I,J+1,K-1)-WW(I,J-1,K-1))
         ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
         S11 = DUDX - ONTHDIV
         S22 = DVDY - ONTHDIV
         S33 = DWDZ - ONTHDIV
         S12 = 0.5_EB*(DUDY+DVDX)
         S13 = 0.5_EB*(DUDZ+DWDX)
         S23 = 0.5_EB*(DVDZ+DWDY)
         M%STRAIN_RATE(I,J,K) = SQRT(2._EB*(S11**2+S22**2+S33**2 &
                                +2._EB*(S12**2+S13**2+S23**2)))
      ENDDO
   ENDDO
ENDDO

! Compute turbulent viscosity (K=K1:K2)

SELECT CASE (TURB_MODEL)

   CASE (NO_TURB_MODEL)

      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               M%MU(I,J,K) = M%MU_DNS(I,J,K)
            ENDDO
         ENDDO
      ENDDO

   CASE (CONSMAG)

      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               M%MU(I,J,K) = M%MU_DNS(I,J,K) &
                  + RHOP(I,J,K)*M%CSD2(I,J,K)*M%STRAIN_RATE(I,J,K)
            ENDDO
         ENDDO
      ENDDO

   CASE (VREMAN)

      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K) &
                      +UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1) &
                      +UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K) &
                      +VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1) &
                      +VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K) &
                      +WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K) &
                      +WW(I,J+1,K-1)-WW(I,J-1,K-1))

               A_IJ(1,1)=DUDX; A_IJ(2,1)=DUDY; A_IJ(3,1)=DUDZ
               A_IJ(1,2)=DVDX; A_IJ(2,2)=DVDY; A_IJ(3,2)=DVDZ
               A_IJ(1,3)=DWDX; A_IJ(2,3)=DWDY; A_IJ(3,3)=DWDZ

               AA=0._EB
               DO JJ=1,3
                  DO II=1,3
                     AA = AA + A_IJ(II,JJ)*A_IJ(II,JJ)
                  ENDDO
               ENDDO

               B_IJ(1,1)=(M%DX(I)*A_IJ(1,1))**2 &
                  + (M%DY(J)*A_IJ(2,1))**2 + (M%DZ(K)*A_IJ(3,1))**2
               B_IJ(2,2)=(M%DX(I)*A_IJ(1,2))**2 &
                  + (M%DY(J)*A_IJ(2,2))**2 + (M%DZ(K)*A_IJ(3,2))**2
               B_IJ(3,3)=(M%DX(I)*A_IJ(1,3))**2 &
                  + (M%DY(J)*A_IJ(2,3))**2 + (M%DZ(K)*A_IJ(3,3))**2

               B_IJ(1,2)=M%DX(I)**2*A_IJ(1,1)*A_IJ(1,2) &
                  + M%DY(J)**2*A_IJ(2,1)*A_IJ(2,2) &
                  + M%DZ(K)**2*A_IJ(3,1)*A_IJ(3,2)
               B_IJ(1,3)=M%DX(I)**2*A_IJ(1,1)*A_IJ(1,3) &
                  + M%DY(J)**2*A_IJ(2,1)*A_IJ(2,3) &
                  + M%DZ(K)**2*A_IJ(3,1)*A_IJ(3,3)
               B_IJ(2,3)=M%DX(I)**2*A_IJ(1,2)*A_IJ(1,3) &
                  + M%DY(J)**2*A_IJ(2,2)*A_IJ(2,3) &
                  + M%DZ(K)**2*A_IJ(3,2)*A_IJ(3,3)

               BB = B_IJ(1,1)*B_IJ(2,2) - B_IJ(1,2)**2 &
                  + B_IJ(1,1)*B_IJ(3,3) - B_IJ(1,3)**2 &
                  + B_IJ(2,2)*B_IJ(3,3) - B_IJ(2,3)**2

               IF (ABS(AA)>TWENTY_EPSILON_EB &
                   .AND. BB>TWENTY_EPSILON_EB) THEN
                  NU_EDDY = C_VREMAN*SQRT(BB/AA)
               ELSE
                  NU_EDDY=0._EB
               ENDIF

               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY

            ENDDO
         ENDDO
      ENDDO

   CASE (WALE)

      DO K=K1,K2
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               DELTA = M%LES_FILTER_WIDTH(I,J,K)
               DUDX = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K) &
                      +UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1) &
                      +UU(I-1,J,K+1)-UU(I-1,J,K-1))
               DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K) &
                      +VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1) &
                      +VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K) &
                      +WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K) &
                      +WW(I,J+1,K-1)-WW(I,J-1,K-1))
               A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
               A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
               A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ
               CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)
               M%MU(I,J,K) = M%MU_DNS(I,J,K) + RHOP(I,J,K)*NU_EDDY
            ENDDO
         ENDDO
      ENDDO

END SELECT

! Compute resolved kinetic energy per unit mass (K=K1:K2)

DO K=K1,K2
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         U2 = 0.25_EB*(UU(I-1,J,K)+UU(I,J,K))**2
         V2 = 0.25_EB*(VV(I,J-1,K)+VV(I,J,K))**2
         W2 = 0.25_EB*(WW(I,J,K-1)+WW(I,J,K))**2
         M%KRES(I,J,K) = 0.5_EB*(U2+V2+W2)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE COMPUTE_VISCOSITY_BLOCK_KERNEL


!> \brief Post-processing for block-decomposed viscosity: wall loops and corner mirroring.
!> \details Must be called sequentially after all blocks have completed.
!> Runs STRAIN_RATE wall corrections, MU wall loop, and MU/KRES corner mirroring.
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables

RECURSIVE SUBROUTINE COMPUTE_VISCOSITY_POST_BLOCK(M,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY
USE TURB_KERNELS, ONLY: WALE_VISCOSITY
USE CC_VELOCITY_KERNELS, ONLY: CUTFACE_VELOCITIES, CC_COMPUTE_KRES, CC_COMPUTE_VISCOSITY

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: NU_EDDY,DELTA,VDF,WGT,A_IJ(3,3), &
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ, &
            S11,S22,S33,S12,S13,S23,ONTHDIV
REAL(EB), PARAMETER :: RAPLUS=1._EB/26._EB
INTEGER :: IIG,JJG,KKG,II,JJ,KK,IW,IOR,IC,SURF_INDEX
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_COUNTER
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(SURFACE_TYPE), POINTER :: SF

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   RHOP => M%RHOS
   UU   => M%US
   VV   => M%VS
   WW   => M%WS
ELSE
   RHOP => M%RHO
   UU   => M%U
   VV   => M%V
   WW   => M%W
ENDIF

! CC_IBM: override KRES at cut cells (must run before wall loops)
IF (CC_IBM) CALL CC_COMPUTE_KRES(M, APPLY_TO_ESTIMATED_VARIABLES)

! Strain rate wall loop (overwrites values near solid walls)

WALL_LOOP_SR: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
   WC=>M%WALL(IW)
   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP_SR

   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   SURF_INDEX = WC%SURF_INDEX
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   IOR = BC%IOR

   IF (IW>M%N_EXTERNAL_WALL_CELLS) THEN
      SELECT CASE(IOR)
         CASE( 1); IF (IIG>M%IBAR) CYCLE WALL_LOOP_SR
         CASE(-1); IF (IIG<1)      CYCLE WALL_LOOP_SR
         CASE( 2); IF (JJG>M%JBAR) CYCLE WALL_LOOP_SR
         CASE(-2); IF (JJG<1)      CYCLE WALL_LOOP_SR
         CASE( 3); IF (KKG>M%KBAR) CYCLE WALL_LOOP_SR
         CASE(-3); IF (KKG<1)      CYCLE WALL_LOOP_SR
      END SELECT
   ENDIF

   DUDX = M%RDX(IIG)*(UU(IIG,JJG,KKG)-UU(IIG-1,JJG,KKG))
   DVDY = M%RDY(JJG)*(VV(IIG,JJG,KKG)-VV(IIG,JJG-1,KKG))
   DWDZ = M%RDZ(KKG)*(WW(IIG,JJG,KKG)-WW(IIG,JJG,KKG-1))
   ONTHDIV = ONTH*(DUDX+DVDY+DWDZ)
   S11 = DUDX - ONTHDIV
   S22 = DVDY - ONTHDIV
   S33 = DWDZ - ONTHDIV

   DUDY = 0.25_EB*M%RDY(JJG)*(UU(IIG,JJG+1,KKG) &
          -UU(IIG,JJG-1,KKG)+UU(IIG-1,JJG+1,KKG) &
          -UU(IIG-1,JJG-1,KKG))
   DUDZ = 0.25_EB*M%RDZ(KKG)*(UU(IIG,JJG,KKG+1) &
          -UU(IIG,JJG,KKG-1)+UU(IIG-1,JJG,KKG+1) &
          -UU(IIG-1,JJG,KKG-1))
   DVDX = 0.25_EB*M%RDX(IIG)*(VV(IIG+1,JJG,KKG) &
          -VV(IIG-1,JJG,KKG)+VV(IIG+1,JJG-1,KKG) &
          -VV(IIG-1,JJG-1,KKG))
   DVDZ = 0.25_EB*M%RDZ(KKG)*(VV(IIG,JJG,KKG+1) &
          -VV(IIG,JJG,KKG-1)+VV(IIG,JJG-1,KKG+1) &
          -VV(IIG,JJG-1,KKG-1))
   DWDX = 0.25_EB*M%RDX(IIG)*(WW(IIG+1,JJG,KKG) &
          -WW(IIG-1,JJG,KKG)+WW(IIG+1,JJG,KKG-1) &
          -WW(IIG-1,JJG,KKG-1))
   DWDY = 0.25_EB*M%RDY(JJG)*(WW(IIG,JJG+1,KKG) &
          -WW(IIG,JJG-1,KKG)+WW(IIG,JJG+1,KKG-1) &
          -WW(IIG,JJG-1,KKG-1))

   S12 = 0.5_EB*(DUDY+DVDX)
   S13 = 0.5_EB*(DUDZ+DWDX)
   S23 = 0.5_EB*(DVDZ+DWDY)

   M%STRAIN_RATE(IIG,JJG,KKG) = SQRT(2._EB*(S11**2+S22**2 &
      +S33**2+2._EB*(S12**2+S13**2+S23**2)))
ENDDO WALL_LOOP_SR

! Mirror viscosity into solids and exterior boundary cells

CELL_COUNTER => M%IWORK1 ; CELL_COUNTER = 0

WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS

   WC=>M%WALL(IW)
   IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP
   BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
   B1=>M%BOUNDARY_PROP1(WC%B1_INDEX)
   B2=>M%BOUNDARY_PROP2(WC%B2_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IC  = M%CELL_INDEX(II,JJ,KK)
   IOR = BC%IOR
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   SF=>SURFACE(WC%SURF_INDEX)

   IF (M%CELL(IC)%SOLID .OR. M%CELL(IC)%EXTERIOR) &
      M%KRES(II,JJ,KK) = M%KRES(IIG,JJG,KKG)

   SELECT CASE(WC%BOUNDARY_TYPE)

      CASE(SOLID_BOUNDARY)

         IF (SIM_MODE/=DNS_MODE) THEN
            DELTA = M%LES_FILTER_WIDTH(IIG,JJG,KKG)
            SELECT CASE(SF%NEAR_WALL_TURB_MODEL)
               CASE DEFAULT
                  NU_EDDY = 0._EB
               CASE(CONSTANT_EDDY_VISCOSITY)
                  NU_EDDY = SF%NEAR_WALL_EDDY_VISCOSITY
               CASE(CONSMAG)
                  VDF = 1._EB-EXP(-B2%Y_PLUS*RAPLUS)
                  NU_EDDY = (VDF*C_SMAGORINSKY*DELTA)**2 &
                     *M%STRAIN_RATE(IIG,JJG,KKG)
               CASE(WALE)
                  DUDX = M%RDX(IIG)*(UU(IIG,JJG,KKG) &
                         -UU(IIG-1,JJG,KKG))
                  DVDY = M%RDY(JJG)*(VV(IIG,JJG,KKG) &
                         -VV(IIG,JJG-1,KKG))
                  DWDZ = M%RDZ(KKG)*(WW(IIG,JJG,KKG) &
                         -WW(IIG,JJG,KKG-1))
                  DUDY = 0.25_EB*M%RDY(JJG) &
                     *(UU(IIG,JJG+1,KKG)-UU(IIG,JJG-1,KKG) &
                      +UU(IIG-1,JJG+1,KKG)-UU(IIG-1,JJG-1,KKG))
                  DUDZ = 0.25_EB*M%RDZ(KKG) &
                     *(UU(IIG,JJG,KKG+1)-UU(IIG,JJG,KKG-1) &
                      +UU(IIG-1,JJG,KKG+1)-UU(IIG-1,JJG,KKG-1))
                  DVDX = 0.25_EB*M%RDX(IIG) &
                     *(VV(IIG+1,JJG,KKG)-VV(IIG-1,JJG,KKG) &
                      +VV(IIG+1,JJG-1,KKG)-VV(IIG-1,JJG-1,KKG))
                  DVDZ = 0.25_EB*M%RDZ(KKG) &
                     *(VV(IIG,JJG,KKG+1)-VV(IIG,JJG,KKG-1) &
                      +VV(IIG,JJG-1,KKG+1)-VV(IIG,JJG-1,KKG-1))
                  DWDX = 0.25_EB*M%RDX(IIG) &
                     *(WW(IIG+1,JJG,KKG)-WW(IIG-1,JJG,KKG) &
                      +WW(IIG+1,JJG,KKG-1)-WW(IIG-1,JJG,KKG-1))
                  DWDY = 0.25_EB*M%RDY(JJG) &
                     *(WW(IIG,JJG+1,KKG)-WW(IIG,JJG-1,KKG) &
                      +WW(IIG,JJG+1,KKG-1)-WW(IIG,JJG-1,KKG-1))
                  A_IJ(1,1)=DUDX; A_IJ(1,2)=DUDY; A_IJ(1,3)=DUDZ
                  A_IJ(2,1)=DVDX; A_IJ(2,2)=DVDY; A_IJ(2,3)=DVDZ
                  A_IJ(3,1)=DWDX; A_IJ(3,2)=DWDY; A_IJ(3,3)=DWDZ
                  CALL WALE_VISCOSITY(NU_EDDY,A_IJ,DELTA)
            END SELECT
            IF (CELL_COUNTER(IIG,JJG,KKG)==0) M%MU(IIG,JJG,KKG) = 0._EB
            CELL_COUNTER(IIG,JJG,KKG) = CELL_COUNTER(IIG,JJG,KKG) + 1
            WGT = 1._EB/REAL(CELL_COUNTER(IIG,JJG,KKG),EB)
            M%MU(IIG,JJG,KKG) = (1._EB-WGT)*M%MU(IIG,JJG,KKG) &
               + WGT*(M%MU_DNS(IIG,JJG,KKG) &
               + RHOP(IIG,JJG,KKG)*NU_EDDY)
         ELSE
            M%MU(IIG,JJG,KKG) = M%MU_DNS(IIG,JJG,KKG)
         ENDIF

         IF (M%CELL(M%CELL_INDEX(II,JJ,KK))%SOLID) &
            M%MU(II,JJ,KK) = M%MU(IIG,JJG,KKG)

      CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY)

         M%MU(II,JJ,KK) = M%MU(IIG,JJG,KKG)

   END SELECT

ENDDO WALL_LOOP

! CC_IBM: compute viscosity on cut-cell region + reset cut-face velocities
IF (CC_IBM) THEN
   CALL CC_COMPUTE_VISCOSITY(M,0._EB)
   CALL CUTFACE_VELOCITIES(M,UU,VV,WW,CUTFACES=.FALSE.)
ENDIF

! Corner mirroring for MU

M%MU(   0,0:M%JBP1,   0) = M%MU(   1,0:M%JBP1,1)
M%MU(M%IBP1,0:M%JBP1,   0) = M%MU(M%IBAR,0:M%JBP1,1)
M%MU(M%IBP1,0:M%JBP1,M%KBP1) = M%MU(M%IBAR,0:M%JBP1,M%KBAR)
M%MU(   0,0:M%JBP1,M%KBP1) = M%MU(   1,0:M%JBP1,M%KBAR)
M%MU(0:M%IBP1,   0,   0) = M%MU(0:M%IBP1,   1,1)
M%MU(0:M%IBP1,M%JBP1,0)    = M%MU(0:M%IBP1,M%JBAR,1)
M%MU(0:M%IBP1,M%JBP1,M%KBP1) = M%MU(0:M%IBP1,M%JBAR,M%KBAR)
M%MU(0:M%IBP1,0,M%KBP1)    = M%MU(0:M%IBP1,   1,M%KBAR)
M%MU(0,   0,0:M%KBP1)    = M%MU(   1,   1,0:M%KBP1)
M%MU(M%IBP1,0,0:M%KBP1)    = M%MU(M%IBAR,   1,0:M%KBP1)
M%MU(M%IBP1,M%JBP1,0:M%KBP1) = M%MU(M%IBAR,M%JBAR,0:M%KBP1)
M%MU(0,M%JBP1,0:M%KBP1)    = M%MU(   1,M%JBAR,0:M%KBP1)

! Corner mirroring for KRES

M%KRES(   0,0:M%JBP1,   0) = M%KRES(   1,0:M%JBP1,1)
M%KRES(M%IBP1,0:M%JBP1,   0) = M%KRES(M%IBAR,0:M%JBP1,1)
M%KRES(M%IBP1,0:M%JBP1,M%KBP1) = M%KRES(M%IBAR,0:M%JBP1,M%KBAR)
M%KRES(   0,0:M%JBP1,M%KBP1) = M%KRES(   1,0:M%JBP1,M%KBAR)
M%KRES(0:M%IBP1,   0,   0) = M%KRES(0:M%IBP1,   1,1)
M%KRES(0:M%IBP1,M%JBP1,0)    = M%KRES(0:M%IBP1,M%JBAR,1)
M%KRES(0:M%IBP1,M%JBP1,M%KBP1) = M%KRES(0:M%IBP1,M%JBAR,M%KBAR)
M%KRES(0:M%IBP1,0,M%KBP1)    = M%KRES(0:M%IBP1,   1,M%KBAR)
M%KRES(0,   0,0:M%KBP1)    = M%KRES(   1,   1,0:M%KBP1)
M%KRES(M%IBP1,0,0:M%KBP1)    = M%KRES(M%IBAR,   1,0:M%KBP1)
M%KRES(M%IBP1,M%JBP1,0:M%KBP1) = M%KRES(M%IBAR,M%JBAR,0:M%KBP1)
M%KRES(0,M%JBP1,0:M%KBP1)    = M%KRES(   1,M%JBAR,0:M%KBP1)

END SUBROUTINE COMPUTE_VISCOSITY_POST_BLOCK


!> \brief Check the CFL and Von Neumann stability criteria.
!> \param M Mesh data structure
!> \param DT Current time step (s)
!> \param DT_NEW_MESH New time step for this mesh (output)
!> \param T Current simulation time (s)
!> \param NM Mesh number (needed for CHANGE_TIME_STEP_INDEX and diagnostic writes)

RECURSIVE SUBROUTINE CHECK_STABILITY_KERNEL(M,DT,DT_NEW_MESH,T,NM)

USE CC_VELOCITY, ONLY : CHECK_CFLVN_LINKED_CELLS
USE OUTPUT_CLOCKS, ONLY: RAMP_TIME_INDEX,RAMP_DT_INDEX
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT,T
REAL(EB), INTENT(OUT) :: DT_NEW_MESH
INTEGER, INTENT(IN) :: NM
REAL(EB) :: UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM,PART_CFL,MU_TMP, UVWMAX_TMP, DT_CLIP
INTEGER  :: I,J,K,IW,IIG,JJG,KKG, ICFL_TMP, JCFL_TMP, KCFL_TMP
REAL(EB), PARAMETER :: DT_EPS = 1.E-10_EB
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(RAMPS_TYPE), POINTER :: RP

UVWMAX = 0._EB
M%VN     = 0._EB
MUTRM  = 1.E-9_EB
R_DX2  = 1.E-9_EB
M%ICFL   = 0; M%JCFL   = 0; M%KCFL   = 0
M%I_VN   = 0; M%J_VN   = 0; M%K_VN   = 0

! Determine max CFL number from all grid cells

UVWMAX_TMP = 0._EB
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
         UODX = MAXVAL(ABS(M%US(I-1:I,J,K)))*M%RDX(I)
         VODY = MAXVAL(ABS(M%VS(I,J-1:J,K)))*M%RDY(J)
         WODZ = MAXVAL(ABS(M%WS(I,J,K-1:K)))*M%RDZ(K)
         SELECT CASE (CFL_VELOCITY_NORM)
            CASE(0) ; UVW = MAX(UODX,VODY,WODZ) + ABS(M%DS(I,J,K))
            CASE(1) ; UVW = UODX + VODY + WODZ  + ABS(M%DS(I,J,K))
            CASE(2) ; UVW = SQRT(UODX**2+VODY**2+WODZ**2) + ABS(M%DS(I,J,K))
            CASE(3) ; UVW = MAX(UODX,VODY,WODZ)
         END SELECT
         IF (UVW>=UVWMAX_TMP) THEN
            UVWMAX_TMP = UVW
            ICFL_TMP = I
            JCFL_TMP = J
            KCFL_TMP = K
         ENDIF
      ENDDO
   ENDDO
ENDDO
IF(UVWMAX_TMP>UVWMAX) THEN
   UVWMAX = UVWMAX_TMP
   M%ICFL = ICFL_TMP
   M%JCFL = JCFL_TMP
   M%KCFL = KCFL_TMP
ENDIF

HEAT_TRANSFER_IF: IF (CHECK_HT) THEN
   WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
      WC=>M%WALL(IW)
      IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE WALL_LOOP
      BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
      B1=>M%BOUNDARY_PROP1(WC%B1_INDEX)
      IIG = BC%IIG
      JJG = BC%JJG
      KKG = BC%KKG
      UVW = (ABS(B1%Q_CON_F)/B1%RHO_F)**ONTH * 2._EB*B1%RDN
      IF (UVW>=UVWMAX) THEN
         UVWMAX = UVW
         M%ICFL=IIG
         M%JCFL=JJG
         M%KCFL=KKG
      ENDIF
   ENDDO WALL_LOOP
ENDIF HEAT_TRANSFER_IF

M%CFL = DT*UVWMAX
! Include surface vegetation drag if necessary
IF (M%DRAG_UVWMAX>0._EB) M%PART_UVWMAX = MAX(M%PART_UVWMAX,M%DRAG_UVWMAX)
PART_CFL = DT*M%PART_UVWMAX

! Determine max Von Neumann Number for fine grid calcs

PARABOLIC_IF: IF (CHECK_VN) THEN

   MU_MAX = 0._EB
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         I_LOOP: DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE I_LOOP
            MU_TMP = MAX(M%D_Z_MAX(I,J,K),M%MU(I,J,K)/M%RHOS(I,J,K))
            IF (MU_TMP>=MU_MAX) THEN
               MU_MAX = MU_TMP
               M%I_VN=I
               M%J_VN=J
               M%K_VN=K
            ENDIF
         ENDDO I_LOOP
      ENDDO
   ENDDO

   IF (TWO_D) THEN
      R_DX2 = M%RDX(M%I_VN)**2 + M%RDZ(M%K_VN)**2
   ELSE
      R_DX2 = M%RDX(M%I_VN)**2 + M%RDY(M%J_VN)**2 + M%RDZ(M%K_VN)**2
   ENDIF

   MUTRM = MU_MAX
   M%VN = DT*2._EB*R_DX2*MUTRM

ENDIF PARABOLIC_IF

IF (CC_IBM) CALL CHECK_CFLVN_LINKED_CELLS(NM,DT,UVWMAX,R_DX2,MUTRM)

! Attempt DT restriction to avoid clippings

DT_CLIP = HUGE(1._EB)
IF (M%CLIP_RHOMIN .OR. M%CLIP_RHOMAX) THEN
   IF (M%DT_RESTRICT_COUNT>=CLIP_DT_RESTRICTIONS_MAX) THEN
      IF (M%CLIP_RHOMIN) WRITE(LU_ERR,'(A,F8.3,A,I0)') 'WARNING: Minimum density, ',RHOMIN,' kg/m3, clipped in Mesh ',NM
      IF (M%CLIP_RHOMAX) WRITE(LU_ERR,'(A,F8.3,A,I0)') 'WARNING: Maximum density, ',RHOMAX,' kg/m3, clipped in Mesh ',NM
   ELSE
      M%CFL = HUGE(1._EB)
      DT_CLIP = DT
      M%DT_RESTRICT_COUNT = M%DT_RESTRICT_COUNT + 1
      M%DT_RESTRICT_STORE = MAX(M%DT_RESTRICT_STORE,M%DT_RESTRICT_COUNT)
   ENDIF
ENDIF

RAMP_TIME_IF: IF (RAMP_TIME_INDEX>0) THEN

   ! User-specified time increments

   RP=>RAMPS(RAMP_TIME_INDEX)
   IF (ICYC==RP%NUMBER_DATA_POINTS) THEN
      DT_NEW_MESH = T_END - RP%INDEPENDENT_DATA(ICYC)
   ELSEIF (ICYC<=RP%NUMBER_DATA_POINTS-1) THEN
      DT_NEW_MESH = RP%INDEPENDENT_DATA(ICYC+1) - RP%INDEPENDENT_DATA(ICYC)
   ELSE
      DT_NEW_MESH = MAX(0._EB,T_END - T)
   ENDIF
   CHANGE_TIME_STEP_INDEX(NM) = 1

ELSE RAMP_TIME_IF

   ! Adjust time step size if necessary

   IF ((M%CFL<CFL_MAX .AND. M%VN<VN_MAX .AND. PART_CFL<PARTICLE_CFL_MAX) .OR. LOCK_TIME_STEP) THEN
      DT_NEW_MESH = DT
      IF (M%CFL<=CFL_MIN .AND. M%VN<VN_MIN .AND. PART_CFL<PARTICLE_CFL_MIN .AND. .NOT.LOCK_TIME_STEP) THEN
         SELECT CASE (RESTRICT_TIME_STEP)
            CASE (.TRUE.);  DT_NEW_MESH = MIN(1.1_EB*DT,DT_INITIAL)
            CASE (.FALSE.); DT_NEW_MESH =     1.1_EB*DT
         END SELECT
         CHANGE_TIME_STEP_INDEX(NM) = 1
      ENDIF
   ELSE
      DT_NEW_MESH = 0.9_EB*MIN( CFL_MAX/MAX(UVWMAX,DT_EPS)               , &
                               VN_MAX/(2._EB*R_DX2*MAX(MUTRM,DT_EPS))   , &
                               PARTICLE_CFL_MAX/MAX(M%PART_UVWMAX,DT_EPS) , &
                               DT_CLIP)
      CHANGE_TIME_STEP_INDEX(NM) = -1
   ENDIF

   IF (RAMP_DT_INDEX > 0) DT_NEW_MESH = MIN(DT_NEW_MESH,EVALUATE_RAMP(T,RAMP_DT_INDEX))

ENDIF RAMP_TIME_IF

END SUBROUTINE CHECK_STABILITY_KERNEL


!> \brief Process edge boundary conditions (parallelizable per-mesh kernel)
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param T Current time (s)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating that estimated (starred) variables are to be used

RECURSIVE SUBROUTINE VELOCITY_BC_PROCESS_EDGES_KERNEL(M,NM,T,APPLY_TO_ESTIMATED_VARIABLES)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE TURB_KERNELS, ONLY: WALL_MODEL
USE CC_VELOCITY, ONLY : GET_OPENBC_TANGENTIAL_CUTFACE_VEL

REAL(EB), INTENT(IN) :: T
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: MUA,TSI,WGT,RAMP_T,OMW,MU_WALL,RHO_WALL,SLIP_COEF,VEL_T, &
            UUP(2),UUM(2),DXX(2),MU_DUIDXJ(-2:2),DUIDXJ(-2:2),PROFILE_FACTOR,VEL_GAS,VEL_GHOST, &
            MU_DUIDXJ_USE(2),DUIDXJ_USE(2),VEL_EDDY,U_TAU,Y_PLUS,U_NORM, &
            DRAG_FACTOR,HT_SCALE_FACTOR,VEG_HT,VEL_N
INTEGER :: NOM(2),IIO(2),JJO(2),KKO(2),IE,II,JJ,KK,IEC,IOR,IWM,IWP,ICMM,ICMP,ICPM,ICPP,ICD,ICDO,IVL,I_SGN, &
           VELOCITY_BC_INDEX,IIGM,JJGM,KKGM,IIGP,JJGP,KKGP,SURF_INDEXM,SURF_INDEXP,ITMP,ICD_SGN,ICDO_SGN, &
           BOUNDARY_TYPE_M,BOUNDARY_TYPE_P,IS,IS2,IWPI,IWMI,VENT_INDEX
LOGICAL :: ALTERED_GRADIENT(-2:2),SYNTHETIC_EDDY_METHOD,HVAC_TANGENTIAL,INTERPOLATED_EDGE,&
           UPWIND_BOUNDARY,INFLOW_BOUNDARY
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,VEL_OTHER
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WCM,WCP,WCX
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: WCM_B1,WCP_B1,WCX_B1
TYPE (EDGE_TYPE), POINTER :: ED
TYPE(SURFACE_TYPE), POINTER :: SF

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==12) RETURN
IF (PERIODIC_TEST==13) RETURN

! Point to the appropriate velocity field

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
   RHOP => M%RHOS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
   RHOP => M%RHO
ENDIF

! Loop over all cell edges and determine the appropriate velocity BCs

EDGE_LOOP: DO IE=1,EDGE_COUNT(NM)

   ED => M%EDGE(IE)

   ED%OMEGA    = -1.E6_EB
   ED%TAU      = -1.E6_EB
   ED%U_AVG    = -1.E6_EB
   ED%V_AVG    = -1.E6_EB
   ED%W_AVG    = -1.E6_EB
   INTERPOLATED_EDGE = .FALSE.

   ! Throw out edges that are completely surrounded by blockages or the exterior of the domain

   ICMM = ED%CELL_INDEX_MM
   ICPM = ED%CELL_INDEX_PM
   ICMP = ED%CELL_INDEX_MP
   ICPP = ED%CELL_INDEX_PP

   IF ((M%CELL(ICMM)%EXTERIOR .OR. M%CELL(ICMM)%SOLID) .AND. &
       (M%CELL(ICPM)%EXTERIOR .OR. M%CELL(ICPM)%SOLID) .AND. &
       (M%CELL(ICMP)%EXTERIOR .OR. M%CELL(ICMP)%SOLID) .AND. &
       (M%CELL(ICPP)%EXTERIOR .OR. M%CELL(ICPP)%SOLID)) CYCLE EDGE_LOOP

   ! Unpack indices for the edge

   II     = ED%I
   JJ     = ED%J
   KK     = ED%K
   IEC    = ED%AXIS
   NOM(1) = ED%NOM_1
   IIO(1) = ED%IIO_1
   JJO(1) = ED%JJO_1
   KKO(1) = ED%KKO_1
   NOM(2) = ED%NOM_2
   IIO(2) = ED%IIO_2
   JJO(2) = ED%JJO_2
   KKO(2) = ED%KKO_2

   ! Get the velocity components at the appropriate cell faces

   COMPONENT: SELECT CASE(IEC)
      CASE(1) COMPONENT
         UUP(1)  = VV(II,JJ,KK+1)
         UUM(1)  = VV(II,JJ,KK)
         UUP(2)  = WW(II,JJ+1,KK)
         UUM(2)  = WW(II,JJ,KK)
         DXX(1)  = M%DY(JJ)
         DXX(2)  = M%DZ(KK)
      CASE(2) COMPONENT
         UUP(1)  = WW(II+1,JJ,KK)
         UUM(1)  = WW(II,JJ,KK)
         UUP(2)  = UU(II,JJ,KK+1)
         UUM(2)  = UU(II,JJ,KK)
         DXX(1)  = M%DZ(KK)
         DXX(2)  = M%DX(II)
      CASE(3) COMPONENT
         UUP(1)  = UU(II,JJ+1,KK)
         UUM(1)  = UU(II,JJ,KK)
         UUP(2)  = VV(II+1,JJ,KK)
         UUM(2)  = VV(II,JJ,KK)
         DXX(1)  = M%DX(II)
         DXX(2)  = M%DY(JJ)
   END SELECT COMPONENT

   ! Indicate that the velocity gradients in the two orthogonal directions have not been changed yet

   ALTERED_GRADIENT = .FALSE.

   ! Loop over all possible orientations of edge and reassign velocity gradients if appropriate

   SIGN_LOOP: DO I_SGN=-1,1,2
      ORIENTATION_LOOP: DO IS=1,3

         IF (IS==IEC) CYCLE ORIENTATION_LOOP

         ! IOR is the orientation of the wall cells adjacent to the edge

         IOR = I_SGN*IS

         ! IS2 is the other coordinate direction besides IOR.

         SELECT CASE(IEC)
            CASE(1)
               IF (IS==2) IS2 = 3
               IF (IS==3) IS2 = 2
            CASE(2)
               IF (IS==1) IS2 = 3
               IF (IS==3) IS2 = 1
            CASE(3)
               IF (IS==1) IS2 = 2
               IF (IS==2) IS2 = 1
            END SELECT

         ! Determine Index_Coordinate_Direction
         ! IEC=1, ICD=1 refers to DWDY; ICD=2 refers to DVDZ
         ! IEC=2, ICD=1 refers to DUDZ; ICD=2 refers to DWDX
         ! IEC=3, ICD=1 refers to DVDX; ICD=2 refers to DUDY

         IF (IS>IEC) ICD = IS-IEC
         IF (IS<IEC) ICD = IS-IEC+3
         ICD_SGN = I_SGN * ICD

         ! IWM and IWP are the wall cell indices of the boundary on either side of the edge.

         IF (IOR<0) THEN
            IWM  = M%CELL(ICMM)%WALL_INDEX(-IOR)
            IWMI = M%CELL(ICMM)%WALL_INDEX( IS2)
            IF (ICD==1) THEN
               IWP  = M%CELL(ICMP)%WALL_INDEX(-IOR)
               IWPI = M%CELL(ICMP)%WALL_INDEX(-IS2)
            ELSE ! ICD==2
               IWP  = M%CELL(ICPM)%WALL_INDEX(-IOR)
               IWPI = M%CELL(ICPM)%WALL_INDEX(-IS2)
            ENDIF
         ELSE
            IF (ICD==1) THEN
               IWM  = M%CELL(ICPM)%WALL_INDEX(-IOR)
               IWMI = M%CELL(ICPM)%WALL_INDEX( IS2)
            ELSE ! ICD==2
               IWM  = M%CELL(ICMP)%WALL_INDEX(-IOR)
               IWMI = M%CELL(ICMP)%WALL_INDEX( IS2)
            ENDIF
            IWP  = M%CELL(ICPP)%WALL_INDEX(-IOR)
            IWPI = M%CELL(ICPP)%WALL_INDEX(-IS2)
         ENDIF

         ! If both adjacent wall cells are undefined, cycle out of the loop.

         IF (IWM==0 .AND. IWP==0) CYCLE ORIENTATION_LOOP

         ! If there is a solid wall separating the two adjacent wall cells, cycle out of the loop.

         IF ((M%WALL(IWMI)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. &
             SURFACE(M%WALL(IWM)%SURF_INDEX)%VELOCITY_BC_INDEX/=FREE_SLIP_BC) .OR. &
             (M%WALL(IWPI)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. &
             SURFACE(M%WALL(IWP)%SURF_INDEX)%VELOCITY_BC_INDEX/=FREE_SLIP_BC)) &
            CYCLE ORIENTATION_LOOP

         ! If only one adjacent wall cell is defined, use its properties.

         IF (IWM>0) THEN
            WCM => M%WALL(IWM)
         ELSE
            WCM => M%WALL(IWP)
         ENDIF

         IF (IWP>0) THEN
            WCP => M%WALL(IWP)
         ELSE
            WCP => M%WALL(IWM)
         ENDIF

         WCM_B1 => M%BOUNDARY_PROP1(WCM%B1_INDEX)
         WCP_B1 => M%BOUNDARY_PROP1(WCP%B1_INDEX)

         ! If both adjacent wall cells are NULL, cycle out.

         BOUNDARY_TYPE_M = WCM%BOUNDARY_TYPE
         BOUNDARY_TYPE_P = WCP%BOUNDARY_TYPE

         IF (BOUNDARY_TYPE_M==NULL_BOUNDARY .AND. BOUNDARY_TYPE_P==NULL_BOUNDARY) CYCLE ORIENTATION_LOOP

         ! Set up synthetic eddy method

         SYNTHETIC_EDDY_METHOD = .FALSE.
         IF (IWM>0 .AND. IWP>0) THEN
            IF (WCM%VENT_INDEX==WCP%VENT_INDEX) THEN
               IF (WCM%VENT_INDEX>0) THEN
                  VT=>M%VENTS(WCM%VENT_INDEX)
                  IF (VT%N_EDDY>0) SYNTHETIC_EDDY_METHOD=.TRUE.
               ENDIF
            ENDIF
         ENDIF

         VEL_EDDY = 0._EB
         SYNTHETIC_EDDY_IF_1: IF (SYNTHETIC_EDDY_METHOD) THEN
            IS_SELECT_1: SELECT CASE(IS) ! unsigned vent orientation
               CASE(1) ! yz plane
                  SELECT CASE(IEC) ! edge orientation
                     CASE(2)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%U_EDDY(JJ,KK)+VT%U_EDDY(JJ,KK+1))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%W_EDDY(JJ,KK)+VT%W_EDDY(JJ,KK+1))
                     CASE(3)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%V_EDDY(JJ,KK)+VT%V_EDDY(JJ+1,KK))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%U_EDDY(JJ,KK)+VT%U_EDDY(JJ+1,KK))
                  END SELECT
               CASE(2) ! zx plane
                  SELECT CASE(IEC)
                     CASE(3)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,KK)+VT%V_EDDY(II+1,KK))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%U_EDDY(II,KK)+VT%U_EDDY(II+1,KK))
                     CASE(1)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,KK)+VT%W_EDDY(II,KK+1))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,KK)+VT%V_EDDY(II,KK+1))
                  END SELECT
               CASE(3) ! xy plane
                  SELECT CASE(IEC)
                     CASE(1)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,JJ)+VT%W_EDDY(II,JJ+1))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%V_EDDY(II,JJ)+VT%V_EDDY(II,JJ+1))
                     CASE(2)
                        IF (ICD==1) VEL_EDDY = 0.5_EB*(VT%U_EDDY(II,JJ)+VT%U_EDDY(II+1,JJ))
                        IF (ICD==2) VEL_EDDY = 0.5_EB*(VT%W_EDDY(II,JJ)+VT%W_EDDY(II+1,JJ))
                  END SELECT
            END SELECT IS_SELECT_1
         ENDIF SYNTHETIC_EDDY_IF_1

         ! OPEN boundary conditions, both varieties, with and without a wind

         OPEN_AND_WIND_BC: IF ((IWM==0 .OR. M%WALL(IWM)%BOUNDARY_TYPE==OPEN_BOUNDARY) .AND. &
                               (IWP==0 .OR. M%WALL(IWP)%BOUNDARY_TYPE==OPEN_BOUNDARY)       ) THEN

            VENT_INDEX = MAX(WCM%VENT_INDEX,WCP%VENT_INDEX)
            VT => M%VENTS(VENT_INDEX)

            UPWIND_BOUNDARY = .FALSE.
            INFLOW_BOUNDARY = .FALSE.

            IF (OPEN_WIND_BOUNDARY) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,   0,KK) + VV(II,   0,KK+1))
                     IF (JJ==M%JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,M%JBAR,KK) + VV(II,M%JBAR,KK+1))
                     IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,0)    + WW(II,JJ+1,   0))
                     IF (KK==M%KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,M%KBAR) + WW(II,JJ+1,M%KBAR))
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(   0,JJ,KK) + UU(   0,JJ,KK+1))
                     IF (II==M%IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(M%IBAR,JJ,KK) + UU(M%IBAR,JJ,KK+1))
                     IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,   0) + WW(II+1,JJ,   0))
                     IF (KK==M%KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,M%KBAR) + WW(II+1,JJ,M%KBAR))
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(   0,JJ,KK) + UU(   0,JJ+1,KK))
                     IF (II==M%IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(M%IBAR,JJ,KK) + UU(M%IBAR,JJ+1,KK))
                     IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,   0,KK) + VV(II+1,   0,KK))
                     IF (JJ==M%JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,M%JBAR,KK) + VV(II+1,M%JBAR,KK))
               END SELECT
               IF ((IOR==1.AND.M%U_WIND(KK)>=0._EB) .OR. (IOR==-1.AND.M%U_WIND(KK)<=0._EB)) &
                  UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==2.AND.M%V_WIND(KK)>=0._EB) .OR. (IOR==-2.AND.M%V_WIND(KK)<=0._EB)) &
                  UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==3.AND.M%W_WIND(KK)>=0._EB) .OR. (IOR==-3.AND.M%W_WIND(KK)<=0._EB)) &
                  UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==1.AND.U_NORM>=0._EB) .OR. (IOR==-1.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
               IF ((IOR==2.AND.U_NORM>=0._EB) .OR. (IOR==-2.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
               IF ((IOR==3.AND.U_NORM>=0._EB) .OR. (IOR==-3.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
            ENDIF

            WIND_NO_WIND_IF: IF (.NOT.UPWIND_BOUNDARY .OR. .NOT.INFLOW_BOUNDARY) THEN

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = WW(II,1,KK)
                     IF (JJ==M%JBAR .AND. IOR==-2) WW(II,M%JBP1,KK) = WW(II,M%JBAR,KK)
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = VV(II,JJ,1)
                     IF (KK==M%KBAR .AND. IOR==-3) VV(II,JJ,M%KBP1) = VV(II,JJ,M%KBAR)
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = WW(1,JJ,KK)
                     IF (II==M%IBAR .AND. IOR==-1) WW(M%IBP1,JJ,KK) = WW(M%IBAR,JJ,KK)
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = UU(II,JJ,1)
                     IF (KK==M%KBAR .AND. IOR==-3) UU(II,JJ,M%KBP1) = UU(II,JJ,M%KBAR)
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = VV(1,JJ,KK)
                     IF (II==M%IBAR .AND. IOR==-1) VV(M%IBP1,JJ,KK) = VV(M%IBAR,JJ,KK)
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = UU(II,1,KK)
                     IF (JJ==M%JBAR .AND. IOR==-2) UU(II,M%JBP1,KK) = UU(II,M%JBAR,KK)
               END SELECT

            ELSE WIND_NO_WIND_IF

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = M%W_WIND(KK) + VEL_EDDY
                     IF (JJ==M%JBAR .AND. IOR==-2) WW(II,M%JBP1,KK) = M%W_WIND(KK) + VEL_EDDY
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = M%V_WIND(KK) + VEL_EDDY
                     IF (KK==M%KBAR .AND. IOR==-3) VV(II,JJ,M%KBP1) = M%V_WIND(KK) + VEL_EDDY
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = M%W_WIND(KK) + VEL_EDDY
                     IF (II==M%IBAR .AND. IOR==-1) WW(M%IBP1,JJ,KK) = M%W_WIND(KK) + VEL_EDDY
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = M%U_WIND(KK) + VEL_EDDY
                     IF (KK==M%KBAR .AND. IOR==-3) UU(II,JJ,M%KBP1) = M%U_WIND(KK) + VEL_EDDY
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = M%V_WIND(KK) + VEL_EDDY
                     IF (II==M%IBAR .AND. IOR==-1) VV(M%IBP1,JJ,KK) = M%V_WIND(KK) + VEL_EDDY
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = M%U_WIND(KK) + VEL_EDDY
                     IF (JJ==M%JBAR .AND. IOR==-2) UU(II,M%JBP1,KK) = M%U_WIND(KK) + VEL_EDDY
               END SELECT

            ENDIF WIND_NO_WIND_IF

            IF (CC_IBM) CALL GET_OPENBC_TANGENTIAL_CUTFACE_VEL(APPLY_TO_ESTIMATED_VARIABLES,UPWIND_BOUNDARY,&
                                                               INFLOW_BOUNDARY,IEC,II,JJ,KK,IOR,UU,VV,WW)

            IF (IWM/=0 .AND. IWP/=0) THEN
               CYCLE EDGE_LOOP  ! Do no further processing of this edge if both cell faces are OPEN
            ELSE
               CYCLE ORIENTATION_LOOP
            ENDIF

         ENDIF OPEN_AND_WIND_BC

         ! Define the appropriate gas and ghost velocity

         IF (ICD==1) THEN ! Used to pick the appropriate velocity component
            IVL=2
         ELSE !ICD==2
            IVL=1
         ENDIF

         IF (IOR<0) THEN
            VEL_GAS   = UUM(IVL)
            VEL_GHOST = UUP(IVL)
            IIGM = M%CELL(ICMM)%I
            JJGM = M%CELL(ICMM)%J
            KKGM = M%CELL(ICMM)%K
            IF (ICD==1) THEN
               IIGP = M%CELL(ICMP)%I
               JJGP = M%CELL(ICMP)%J
               KKGP = M%CELL(ICMP)%K
            ELSE ! ICD==2
               IIGP = M%CELL(ICPM)%I
               JJGP = M%CELL(ICPM)%J
               KKGP = M%CELL(ICPM)%K
            ENDIF
         ELSE
            VEL_GAS   = UUP(IVL)
            VEL_GHOST = UUM(IVL)
            IF (ICD==1) THEN
               IIGM = M%CELL(ICPM)%I
               JJGM = M%CELL(ICPM)%J
               KKGM = M%CELL(ICPM)%K
            ELSE ! ICD==2
               IIGM = M%CELL(ICMP)%I
               JJGM = M%CELL(ICMP)%J
               KKGM = M%CELL(ICMP)%K
            ENDIF
            IIGP = M%CELL(ICPP)%I
            JJGP = M%CELL(ICPP)%J
            KKGP = M%CELL(ICPP)%K
         ENDIF

         ! Decide whether or not to process edge using data interpolated from another mesh

         INTERPOLATION_IF: IF (NOM(ICD)==0 .OR. &
                   (BOUNDARY_TYPE_M==SOLID_BOUNDARY .OR. BOUNDARY_TYPE_P==SOLID_BOUNDARY) .OR. &
                   (BOUNDARY_TYPE_M/=INTERPOLATED_BOUNDARY .AND. BOUNDARY_TYPE_P/=INTERPOLATED_BOUNDARY) .OR. &
                   (SYNTHETIC_EDDY_METHOD .AND. &
                   (BOUNDARY_TYPE_M==OPEN_BOUNDARY .OR. BOUNDARY_TYPE_P==OPEN_BOUNDARY)) ) THEN

            ! Determine appropriate velocity BC by assessing each adjacent wall cell.
            ! If the BCs are different on each side of the edge, choose the one with the
            ! specified velocity or velocity gradient, if there is one. If not, choose the
            ! max value of boundary condition index, simply for consistency.

            SURF_INDEXM = WCM%SURF_INDEX
            SURF_INDEXP = WCP%SURF_INDEX
            IF (SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_VELOCITY .OR. &
                SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXM)
            ELSEIF (SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_VELOCITY .OR. &
                    SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXP)
            ELSE
               SF=>SURFACE(MAX(SURF_INDEXM,SURF_INDEXP))
            ENDIF
            VELOCITY_BC_INDEX = SF%VELOCITY_BC_INDEX
            IF (WCM%VENT_INDEX==WCP%VENT_INDEX .AND. WCP%VENT_INDEX > 0) THEN
               IF(M%VENTS(WCM%VENT_INDEX)%NODE_INDEX>0 .AND. WCM_B1%U_NORMAL >= 0._EB) &
                  VELOCITY_BC_INDEX=FREE_SLIP_BC
            ENDIF
            IF (SYNTHETIC_EDDY_METHOD)         VELOCITY_BC_INDEX=NO_SLIP_BC
            IF (SF%ROUGHNESS>2._EB/WCM_B1%RDN) VELOCITY_BC_INDEX=NO_SLIP_BC ! see Basu et al. BLM 2017

            ! Compute the viscosity by averaging the two adjacent gas cells

            MUA = 0.5_EB*(M%MU(IIGM,JJGM,KKGM) + M%MU(IIGP,JJGP,KKGP))

            ! Check for HVAC tangential velocity

            HVAC_TANGENTIAL = .FALSE.
            IF (WCM%VENT_INDEX>0 .OR. WCP%VENT_INDEX>0) THEN
               IF (WCM%VENT_INDEX>0) THEN
                  WCX => WCM
               ELSE
                  WCX => WCP
               ENDIF
               VT => M%VENTS(WCX%VENT_INDEX)
               WCX_B1 => M%BOUNDARY_PROP1(WCX%B1_INDEX)
               IF (VT%NODE_INDEX>0 .AND. WCX_B1%U_NORMAL_S<0._EB) THEN
                  VELOCITY_BC_INDEX = NO_SLIP_BC
                  IF (ALL(VT%UVW>-1.E12_EB)) HVAC_TANGENTIAL = .TRUE.
               ENDIF
            ENDIF

            ! Determine if there is a tangential velocity component

            IF (.NOT.SF%SPECIFIED_TANGENTIAL_VELOCITY .AND. .NOT.SYNTHETIC_EDDY_METHOD .AND. &
                .NOT.HVAC_TANGENTIAL) THEN

               VEL_T = 0._EB

            ELSEIF (HVAC_TANGENTIAL) THEN

               VEL_T = 0._EB
               SELECT CASE(IEC) ! edge orientation
                  CASE (1)
                     IF (ICD==1) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(3)
                     IF (ICD==2) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(2)
                  CASE (2)
                     IF (ICD==1) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(1)
                     IF (ICD==2) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(3)
                  CASE (3)
                     IF (ICD==1) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(2)
                     IF (ICD==2) VEL_T = ABS(WCX_B1%U_NORMAL_S/VT%UVW(ABS(VT%IOR)))*VT%UVW(1)
               END SELECT

            ELSE

               VEL_N = 0.5_EB*(WCM_B1%U_NORMAL_S+WCP_B1%U_NORMAL_S)

               IF (ABS(SF%VEL)>0._EB .OR. VEL_N==0._EB) THEN
                  IF (ABS(SF%T_IGN-T_BEGIN)<=SPACING(SF%T_IGN) .AND. SF%RAMP(TIME_VELO)%INDEX>=1) THEN
                     TSI = T
                  ELSE
                     TSI=T-SF%T_IGN
                  ENDIF
                  PROFILE_FACTOR = 1._EB
                  RAMP_T = EVALUATE_RAMP(TSI,SF%RAMP(TIME_VELO)%INDEX,TAU=SF%RAMP(TIME_VELO)%TAU)
                  IF (SF%VEL < 0._EB) THEN
                     IF (SF%RAMP(VELO_PROF_Z)%INDEX>0) &
                        PROFILE_FACTOR = EVALUATE_RAMP(M%ZC(KK),SF%RAMP(VELO_PROF_Z)%INDEX)
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) &
                        VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(2) + VEL_EDDY))
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) &
                        VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(1) + VEL_EDDY))
                  ELSEIF (SF%VEL > 0._EB) THEN
                     IF (SF%PROFILE/=0) &
                        PROFILE_FACTOR = ABS(0.5_EB*(WCM_B1%U_NORMAL_0+WCP_B1%U_NORMAL_0)/SF%VEL)
                     IF (SF%RAMP(VELO_PROF_Z)%INDEX>0) &
                        PROFILE_FACTOR = EVALUATE_RAMP(M%ZC(KK),SF%RAMP(VELO_PROF_Z)%INDEX)
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) &
                        VEL_T = RAMP_T*PROFILE_FACTOR*VEL_EDDY
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) &
                        VEL_T = RAMP_T*PROFILE_FACTOR*VEL_EDDY
                  ELSE  ! User-specified VEL_T but with VEL=0
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = RAMP_T*SF%VEL_T(2)
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = RAMP_T*SF%VEL_T(1)
                  ENDIF
               ELSE ! VEL_N is due to something else besides a user-specified VEL, like a MASS_FLUX BC
                  IF (VEL_N < 0._EB) THEN
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = -SF%VEL_T(2)*VEL_N
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = -SF%VEL_T(1)*VEL_N
                  ENDIF
               ENDIF

            ENDIF

            ! Choose the appropriate boundary condition to apply

            BOUNDARY_CONDITION: SELECT CASE(VELOCITY_BC_INDEX)

               CASE (FREE_SLIP_BC) BOUNDARY_CONDITION

                  VEL_GHOST = VEL_GAS
                  DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                  MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
                  ALTERED_GRADIENT(ICD_SGN) = .TRUE.

               CASE (NO_SLIP_BC) BOUNDARY_CONDITION

                  VEL_GHOST = 2._EB*VEL_T - VEL_GAS
                  DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                  MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
                  ALTERED_GRADIENT(ICD_SGN) = .TRUE.

               CASE (WALL_MODEL_BC) BOUNDARY_CONDITION

                  ! SLIP_COEF = -1, no slip,   VEL_GHOST = 2*VEL_T - VEL_GAS
                  ! SLIP_COEF =  0, half slip, VEL_GHOST = VEL_T
                  ! SLIP_COEF =  1, free slip, VEL_GHOST = VEL_GAS

                  IF ((IWM==0.OR.IWP==0) .AND. .NOT.ED%EXTERNAL) THEN  ! Special case for a corner
                     VEL_GHOST = 2._EB*VEL_T - VEL_GAS
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = MUA*DUIDXJ(ICD_SGN)
                  ELSE
                     ITMP = MIN(I_MAX_TEMP,NINT(0.5_EB*(M%TMP(IIGM,JJGM,KKGM)+M%TMP(IIGP,JJGP,KKGP))))
                     MU_WALL = MU_RSQMW_Z(ITMP,1)/RSQ_MW_Z(1)
                     RHO_WALL = 0.5_EB*( RHOP(IIGM,JJGM,KKGM) + RHOP(IIGP,JJGP,KKGP) )
                     CALL WALL_MODEL(SLIP_COEF,U_TAU,Y_PLUS,MU_WALL/RHO_WALL,SF%ROUGHNESS, &
                                     0.5_EB*DXX(ICD),VEL_GAS-VEL_T)
                     VEL_GHOST = VEL_T + SLIP_COEF*(VEL_GAS-VEL_T)
                     DUIDXJ(ICD_SGN) = I_SGN*(VEL_GAS-VEL_GHOST)/DXX(ICD)
                     MU_DUIDXJ(ICD_SGN) = RHO_WALL*U_TAU**2 * SIGN(1._EB,DUIDXJ(ICD_SGN))
                  ENDIF
                  ALTERED_GRADIENT(ICD_SGN) = .TRUE.

               CASE (BOUNDARY_FUEL_MODEL_BC) BOUNDARY_CONDITION

                  RHO_WALL = 0.5_EB*( RHOP(IIGM,JJGM,KKGM) + RHOP(IIGP,JJGP,KKGP) )
                  VEL_T = SQRT(UU(IIGM,JJGM,KKGM)**2 + VV(IIGM,JJGM,KKGM)**2)
                  VEL_GHOST = 2._EB*VEL_T - VEL_GAS
                  DUIDXJ(ICD_SGN) = 0._EB
                  IF (SF%VEG_LSET_SPREAD) THEN
                     VEG_HT = SF%VEG_LSET_HT
                     DRAG_FACTOR = 0.5_EB*SF%DRAG_COEFFICIENT*SF%SHAPE_FACTOR*SF%VEG_LSET_BETA* &
                                   (SF%VEG_LSET_SIGMA*100._EB)
                  ELSE
                     VEG_HT = SF%LAYER_THICKNESS(1)
                     DRAG_FACTOR = 0.5_EB*SF%DRAG_COEFFICIENT*SF%SHAPE_FACTOR*SF%PACKING_RATIO(1)* &
                                   SF%SURFACE_VOLUME_RATIO(1)
                  ENDIF
                  HT_SCALE_FACTOR = MIN(1._EB,0.5_EB*(WCM_B1%RDN+WCP_B1%RDN)*VEG_HT)
                  MU_DUIDXJ(ICD_SGN) = I_SGN*RHO_WALL*DRAG_FACTOR*VEG_HT*HT_SCALE_FACTOR**2* &
                                       VEL_GAS*VEL_T
                  M%DRAG_UVWMAX = MAX(M%DRAG_UVWMAX,DRAG_FACTOR*HT_SCALE_FACTOR**2*VEL_T)
                  ALTERED_GRADIENT(ICD_SGN) = .TRUE.

            END SELECT BOUNDARY_CONDITION

         ELSE INTERPOLATION_IF  ! Use data from another mesh

            INTERPOLATED_EDGE = .TRUE.
            OM => M%OMESH(ABS(NOM(ICD)))

            IF (PREDICTOR) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%WS
                     ELSE ! ICD=2
                        VEL_OTHER => OM%VS
                     ENDIF
                  CASE(2)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%US
                     ELSE ! ICD=2
                        VEL_OTHER => OM%WS
                     ENDIF
                  CASE(3)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%VS
                     ELSE ! ICD=2
                        VEL_OTHER => OM%US
                     ENDIF
               END SELECT
            ELSE
               SELECT CASE(IEC)
                  CASE(1)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%W
                     ELSE ! ICD=2
                        VEL_OTHER => OM%V
                     ENDIF
                  CASE(2)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%U
                     ELSE ! ICD=2
                        VEL_OTHER => OM%W
                     ENDIF
                  CASE(3)
                     IF (ICD==1) THEN
                        VEL_OTHER => OM%V
                     ELSE ! ICD=2
                        VEL_OTHER => OM%U
                     ENDIF
               END SELECT
            ENDIF

            WGT = ED%EDGE_INTERPOLATION_FACTOR(ICD)
            OMW = 1._EB-WGT

            SELECT CASE(IEC)
               CASE(1)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ENDIF
               CASE(2)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ENDIF
               CASE(3)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ELSE ! ICD==2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + &
                                 OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ENDIF
            END SELECT

            ! At the exterior edge of Mesh NM, which abuts Mesh NOM, assign the appropriate
            ! velocity component to the ghost cell.

            IF (CORRECTOR) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. KK==0    .AND. ABS(IOR)==2) &
                        UU(II,JJ  ,KK  ) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)-1)
                     IF (JJ==0    .AND. KK==0    .AND. ABS(IOR)==3) &
                        UU(II,JJ  ,KK  ) = OM%U(IIO(ICD),JJO(ICD)-1,KKO(ICD)  )
                     IF (JJ==0    .AND. KK==M%KBAR .AND. ABS(IOR)==2) &
                        UU(II,JJ  ,KK+1) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)+1)
                     IF (JJ==0    .AND. KK==M%KBAR .AND. ABS(IOR)==3) &
                        UU(II,JJ  ,KK+1) = OM%U(IIO(ICD),JJO(ICD)-1,KKO(ICD)  )
                     IF (JJ==M%JBAR .AND. KK==0    .AND. ABS(IOR)==2) &
                        UU(II,JJ+1,KK  ) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)-1)
                     IF (JJ==M%JBAR .AND. KK==0    .AND. ABS(IOR)==3) &
                        UU(II,JJ+1,KK  ) = OM%U(IIO(ICD),JJO(ICD)+1,KKO(ICD)  )
                     IF (JJ==M%JBAR .AND. KK==M%KBAR .AND. ABS(IOR)==2) &
                        UU(II,JJ+1,KK+1) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)+1)
                     IF (JJ==M%JBAR .AND. KK==M%KBAR .AND. ABS(IOR)==3) &
                        UU(II,JJ+1,KK+1) = OM%U(IIO(ICD),JJO(ICD)+1,KKO(ICD)  )
                  CASE(2)
                     IF (II==0    .AND. KK==0    .AND. ABS(IOR)==1) &
                        VV(II  ,JJ,KK  ) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)-1)
                     IF (II==0    .AND. KK==0    .AND. ABS(IOR)==3) &
                        VV(II  ,JJ,KK  ) = OM%V(IIO(ICD)-1,JJO(ICD),KKO(ICD)  )
                     IF (II==0    .AND. KK==M%KBAR .AND. ABS(IOR)==1) &
                        VV(II  ,JJ,KK+1) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)+1)
                     IF (II==0    .AND. KK==M%KBAR .AND. ABS(IOR)==3) &
                        VV(II  ,JJ,KK+1) = OM%V(IIO(ICD)-1,JJO(ICD),KKO(ICD)  )
                     IF (II==M%IBAR .AND. KK==0    .AND. ABS(IOR)==1) &
                        VV(II+1,JJ,KK  ) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)-1)
                     IF (II==M%IBAR .AND. KK==0    .AND. ABS(IOR)==3) &
                        VV(II+1,JJ,KK  ) = OM%V(IIO(ICD)+1,JJO(ICD),KKO(ICD)  )
                     IF (II==M%IBAR .AND. KK==M%KBAR .AND. ABS(IOR)==1) &
                        VV(II+1,JJ,KK+1) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)+1)
                     IF (II==M%IBAR .AND. KK==M%KBAR .AND. ABS(IOR)==3) &
                        VV(II+1,JJ,KK+1) = OM%V(IIO(ICD)+1,JJO(ICD),KKO(ICD)  )
                  CASE(3)
                     IF (II==0    .AND. JJ==0    .AND. ABS(IOR)==1) &
                        WW(II  ,JJ  ,KK) = OM%W(IIO(ICD)  ,JJO(ICD)-1,KKO(ICD))
                     IF (II==0    .AND. JJ==0    .AND. ABS(IOR)==2) &
                        WW(II  ,JJ  ,KK) = OM%W(IIO(ICD)-1,JJO(ICD)  ,KKO(ICD))
                     IF (II==0    .AND. JJ==M%JBAR .AND. ABS(IOR)==1) &
                        WW(II  ,JJ+1,KK) = OM%W(IIO(ICD)  ,JJO(ICD)+1,KKO(ICD))
                     IF (II==0    .AND. JJ==M%JBAR .AND. ABS(IOR)==2) &
                        WW(II  ,JJ+1,KK) = OM%W(IIO(ICD)-1,JJO(ICD)  ,KKO(ICD))
                     IF (II==M%IBAR .AND. JJ==0    .AND. ABS(IOR)==1) &
                        WW(II+1,JJ  ,KK) = OM%W(IIO(ICD)  ,JJO(ICD)-1,KKO(ICD))
                     IF (II==M%IBAR .AND. JJ==0    .AND. ABS(IOR)==2) &
                        WW(II+1,JJ  ,KK) = OM%W(IIO(ICD)+1,JJO(ICD)  ,KKO(ICD))
                     IF (II==M%IBAR .AND. JJ==M%JBAR .AND. ABS(IOR)==1) &
                        WW(II+1,JJ+1,KK) = OM%W(IIO(ICD)  ,JJO(ICD)+1,KKO(ICD))
                     IF (II==M%IBAR .AND. JJ==M%JBAR .AND. ABS(IOR)==2) &
                        WW(II+1,JJ+1,KK) = OM%W(IIO(ICD)+1,JJO(ICD)  ,KKO(ICD))
               END SELECT
            ENDIF

         ENDIF INTERPOLATION_IF

         ! Set ghost cell values at edge of computational domain

         SELECT CASE(IEC)
            CASE(1)
               IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
               IF (JJ==M%JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
               IF (KK==M%KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. .NOT.INTERPOLATED_EDGE) THEN
                 IF (ICD==1) THEN
                    ED%W_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    ED%V_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(2)
               IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
               IF (II==M%IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
               IF (KK==M%KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. .NOT.INTERPOLATED_EDGE) THEN
                 IF (ICD==1) THEN
                    ED%U_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    ED%W_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(3)
               IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
               IF (II==M%IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
               IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
               IF (JJ==M%JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
               IF (CORRECTOR .AND. .NOT.INTERPOLATED_EDGE) THEN
                 IF (ICD==1) THEN
                    ED%V_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    ED%U_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
         END SELECT

      ENDDO ORIENTATION_LOOP
   ENDDO SIGN_LOOP

   ! Cycle out of the EDGE_LOOP if no tangential gradients have been altered.

   IF (.NOT.ANY(ALTERED_GRADIENT)) CYCLE EDGE_LOOP

   ! If the edge is on an interpolated boundary, and all cells around it are not solid, cycle

   IF (INTERPOLATED_EDGE) THEN
      IF (.NOT.M%CELL(ICMM)%SOLID .AND. .NOT.M%CELL(ICPM)%SOLID .AND. &
          .NOT.M%CELL(ICMP)%SOLID .AND. .NOT.M%CELL(ICPP)%SOLID) CYCLE EDGE_LOOP
   ENDIF

   ! Loop over all 4 normal directions and compute vorticity and stress tensor components for each

   SIGN_LOOP_2: DO I_SGN=-1,1,2
      ORIENTATION_LOOP_2: DO ICD=1,2
         IF (ICD==1) THEN
            ICDO=2
         ELSE ! ICD=2
            ICDO=1
         ENDIF
         ICD_SGN = I_SGN*ICD
         IF (ALTERED_GRADIENT(ICD_SGN)) THEN
               DUIDXJ_USE(ICD) =    DUIDXJ(ICD_SGN)
            MU_DUIDXJ_USE(ICD) = MU_DUIDXJ(ICD_SGN)
         ELSEIF (ALTERED_GRADIENT(-ICD_SGN)) THEN
               DUIDXJ_USE(ICD) =    DUIDXJ(-ICD_SGN)
            MU_DUIDXJ_USE(ICD) = MU_DUIDXJ(-ICD_SGN)
         ELSE
            CYCLE ORIENTATION_LOOP_2
         ENDIF
         ICDO_SGN = I_SGN*ICDO
         IF (ALTERED_GRADIENT(ICDO_SGN)) THEN
               DUIDXJ_USE(ICDO) =    DUIDXJ(ICDO_SGN)
            MU_DUIDXJ_USE(ICDO) = MU_DUIDXJ(ICDO_SGN)
         ELSEIF (ALTERED_GRADIENT(-ICDO_SGN)) THEN
               DUIDXJ_USE(ICDO) =    DUIDXJ(-ICDO_SGN)
            MU_DUIDXJ_USE(ICDO) = MU_DUIDXJ(-ICDO_SGN)
         ELSE
               DUIDXJ_USE(ICDO) = 0._EB
            MU_DUIDXJ_USE(ICDO) = 0._EB
         ENDIF
         ED%OMEGA(ICD_SGN) =    DUIDXJ_USE(1) -    DUIDXJ_USE(2)
         ED%TAU(ICD_SGN)   = MU_DUIDXJ_USE(1) + MU_DUIDXJ_USE(2)
      ENDDO ORIENTATION_LOOP_2
   ENDDO SIGN_LOOP_2

ENDDO EDGE_LOOP
END SUBROUTINE VELOCITY_BC_PROCESS_EDGES_KERNEL


!> \brief Thread-safe kernel version of VISCOSITY_BC.
!> Fills ghost cells of MU, KRES, D/DS from neighboring mesh OMESH data.
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Use estimated (starred) variables

SUBROUTINE VISCOSITY_BC_KERNEL(M,NM,APPLY_TO_ESTIMATED_VARIABLES)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: MU_OTHER,DP_OTHER,KRES_OTHER
INTEGER :: II,JJ,KK,IW,IIO,JJO,KKO,NOM,N_INT_CELLS
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
   WC =>M%WALL(IW)
   EWC=>M%EXTERNAL_WALL(IW)
   IF (EWC%NOM==0) CYCLE WALL_LOOP
   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   NOM = EWC%NOM
   MU_OTHER   = 0._EB
   DP_OTHER   = 0._EB
   KRES_OTHER = 0._EB
   DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
      DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
         DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
            MU_OTHER = MU_OTHER + M%OMESH(NOM)%MU(IIO,JJO,KKO)
            KRES_OTHER = KRES_OTHER + M%OMESH(NOM)%KRES(IIO,JJO,KKO)
            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
               DP_OTHER = DP_OTHER + M%OMESH(NOM)%DS(IIO,JJO,KKO)
            ELSE
               DP_OTHER = DP_OTHER + M%OMESH(NOM)%D(IIO,JJO,KKO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
   MU_OTHER = MU_OTHER/REAL(N_INT_CELLS,EB)
   KRES_OTHER = KRES_OTHER/REAL(N_INT_CELLS,EB)
   DP_OTHER = DP_OTHER/REAL(N_INT_CELLS,EB)
   M%MU(II,JJ,KK) = MU_OTHER
   M%KRES(II,JJ,KK) = KRES_OTHER
   IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
      M%DS(II,JJ,KK) = DP_OTHER
   ELSE
      M%D(II,JJ,KK) = DP_OTHER
   ENDIF
ENDDO WALL_LOOP

END SUBROUTINE VISCOSITY_BC_KERNEL


!> \brief Thread-safe kernel version of MATCH_VELOCITY.
!> Forces normal component of velocity to match at interpolated boundaries.
!> Handles non-CC_IBM case only; CC_IBM dispatches to CC_MATCH_VELOCITY at wrapper level.
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param PREDICTOR_FLAG .TRUE. for predictor phase, .FALSE. for corrector

SUBROUTINE MATCH_VELOCITY_KERNEL(M,NM,PREDICTOR_FLAG)

USE COMPLEX_GEOMETRY, ONLY : CC_IDCF
USE MESH_VARIABLES, ONLY: MESHES

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: PREDICTOR_FLAG

INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
REAL(EB) :: DA_OTHER,UU_OTHER,VV_OTHER,WW_OTHER,NOM_CELLS
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,OM_UU,OM_VV,OM_WW
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

INTEGER  :: ICF
REAL(EB) :: AU,AU1,AV,AV1,AW,AW1

! Point to the appropriate velocity field

IF (PREDICTOR_FLAG) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
ENDIF

! Loop over all external wall cells and force adjacent normal
! components of velocity at interpolated boundaries to match.

EXTERNAL_WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS

   WC=>M%WALL(IW)
   EWC=>M%EXTERNAL_WALL(IW)
   EWC%BOUNDARY_TYPE_PREVIOUS = WC%BOUNDARY_TYPE

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   BC =>M%BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR
   NOM = EWC%NOM
   OM => M%OMESH(NOM)
   M2 => MESHES(NOM)

   ! Determine the area of the interpolated cell face

   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         IF (PREDICTOR_FLAG) OM_UU => OM%US
         IF (.NOT.PREDICTOR_FLAG) OM_UU => OM%U
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         IF (PREDICTOR_FLAG) OM_VV => OM%VS
         IF (.NOT.PREDICTOR_FLAG) OM_VV => OM%V
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         IF (PREDICTOR_FLAG) OM_WW => OM%WS
         IF (.NOT.PREDICTOR_FLAG) OM_WW => OM%W
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   ! Determine the normal component of velocity from the other mesh

   SELECT CASE(IOR)

      CASE( 1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + OM_UU(IIO,JJO,KKO) &
                     *M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_UU(IIO,JJO,KKO) = &
                     0.5_EB*(OM_UU(IIO,JJO,KKO)+UU(0,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = UU(0,JJ,KK)
         UU(0,JJ,KK) = 0.5_EB*(UU(0,JJ,KK) + UU_OTHER)

      CASE(-1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + &
                     OM_UU(IIO-1,JJO,KKO) &
                     *M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_UU(IIO-1,JJO,KKO) = &
                     0.5_EB*(OM_UU(IIO-1,JJO,KKO) &
                     +UU(M%IBAR,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = UU(M%IBAR,JJ,KK)
         UU(M%IBAR,JJ,KK) = &
            0.5_EB*(UU(M%IBAR,JJ,KK) + UU_OTHER)

      CASE( 2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO,KKO) &
                     *M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_VV(IIO,JJO,KKO) = &
                     0.5_EB*(OM_VV(IIO,JJO,KKO)+VV(II,0,KK))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = VV(II,0,KK)
         VV(II,0,KK) = 0.5_EB*(VV(II,0,KK) + VV_OTHER)

      CASE(-2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + &
                     OM_VV(IIO,JJO-1,KKO) &
                     *M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_VV(IIO,JJO-1,KKO) = &
                     0.5_EB*(OM_VV(IIO,JJO-1,KKO) &
                     +VV(II,M%JBAR,KK))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = VV(II,M%JBAR,KK)
         VV(II,M%JBAR,KK) = &
            0.5_EB*(VV(II,M%JBAR,KK) + VV_OTHER)

      CASE( 3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO) &
                     *M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_WW(IIO,JJO,KKO) = &
                     0.5_EB*(OM_WW(IIO,JJO,KKO)+WW(II,JJ,0))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = WW(II,JJ,0)
         WW(II,JJ,0) = 0.5_EB*(WW(II,JJ,0) + WW_OTHER)

      CASE(-3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + &
                     OM_WW(IIO,JJO,KKO-1) &
                     *M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) &
                     OM_WW(IIO,JJO,KKO-1) = &
                     0.5_EB*(OM_WW(IIO,JJO,KKO-1) &
                     +WW(II,JJ,M%KBAR))
               ENDDO
            ENDDO
         ENDDO
         M%UVW_SAVE(IW) = WW(II,JJ,M%KBAR)
         WW(II,JJ,M%KBAR) = &
            0.5_EB*(WW(II,JJ,M%KBAR) + WW_OTHER)

   END SELECT

   ! Save velocity components at the ghost cell midpoint

   M%U_GHOST(IW) = 0._EB
   M%V_GHOST(IW) = 0._EB
   M%W_GHOST(IW) = 0._EB

   IF (PREDICTOR_FLAG) OM_UU => OM%US
   IF (.NOT.PREDICTOR_FLAG) OM_UU => OM%U
   IF (PREDICTOR_FLAG) OM_VV => OM%VS
   IF (.NOT.PREDICTOR_FLAG) OM_VV => OM%V
   IF (PREDICTOR_FLAG) OM_WW => OM%WS
   IF (.NOT.PREDICTOR_FLAG) OM_WW => OM%W

   IF (CC_IBM) THEN
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               AU =1._EB
               ICF=M2%FCVAR(IIO  ,JJO,KKO,CC_IDCF,IAXIS)
               IF(ICF>0) AU =M2%CUT_FACE(ICF)%ALPHA_CF
               AU1=1._EB
               ICF=M2%FCVAR(IIO-1,JJO,KKO,CC_IDCF,IAXIS)
               IF(ICF>0) AU1=M2%CUT_FACE(ICF)%ALPHA_CF
               AV =1._EB
               ICF=M2%FCVAR(IIO,JJO  ,KKO,CC_IDCF,JAXIS)
               IF(ICF>0) AV =M2%CUT_FACE(ICF)%ALPHA_CF
               AV1=1._EB
               ICF=M2%FCVAR(IIO,JJO-1,KKO,CC_IDCF,JAXIS)
               IF(ICF>0) AV1=M2%CUT_FACE(ICF)%ALPHA_CF
               AW =1._EB
               ICF=M2%FCVAR(IIO,JJO,KKO  ,CC_IDCF,KAXIS)
               IF(ICF>0) AW =M2%CUT_FACE(ICF)%ALPHA_CF
               AW1=1._EB
               ICF=M2%FCVAR(IIO,JJO,KKO-1,CC_IDCF,KAXIS)
               IF(ICF>0) AW1=M2%CUT_FACE(ICF)%ALPHA_CF
               M%U_GHOST(IW) = M%U_GHOST(IW) + &
                  (OM_UU(IIO,JJO,KKO) &
                  +OM_UU(IIO-1,JJO,KKO))/(AU+AU1)
               M%V_GHOST(IW) = M%V_GHOST(IW) + &
                  (OM_VV(IIO,JJO,KKO) &
                  +OM_VV(IIO,JJO-1,KKO))/(AV+AV1)
               M%W_GHOST(IW) = M%W_GHOST(IW) + &
                  (OM_WW(IIO,JJO,KKO) &
                  +OM_WW(IIO,JJO,KKO-1))/(AW+AW1)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               M%U_GHOST(IW) = M%U_GHOST(IW) + &
                  0.5_EB*(OM_UU(IIO,JJO,KKO) &
                  +OM_UU(IIO-1,JJO,KKO))
               M%V_GHOST(IW) = M%V_GHOST(IW) + &
                  0.5_EB*(OM_VV(IIO,JJO,KKO) &
                  +OM_VV(IIO,JJO-1,KKO))
               M%W_GHOST(IW) = M%W_GHOST(IW) + &
                  0.5_EB*(OM_WW(IIO,JJO,KKO) &
                  +OM_WW(IIO,JJO,KKO-1))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   NOM_CELLS = REAL((EWC%IIO_MAX-EWC%IIO_MIN+1) &
      *(EWC%JJO_MAX-EWC%JJO_MIN+1) &
      *(EWC%KKO_MAX-EWC%KKO_MIN+1),EB)
   M%U_GHOST(IW) = M%U_GHOST(IW)/NOM_CELLS
   M%V_GHOST(IW) = M%V_GHOST(IW)/NOM_CELLS
   M%W_GHOST(IW) = M%W_GHOST(IW)/NOM_CELLS

ENDDO EXTERNAL_WALL_LOOP

END SUBROUTINE MATCH_VELOCITY_KERNEL


!> \brief Apply no-flux boundary conditions for pressure solver.
!> \details Sets velocity flux (FVX/FVY/FVZ) at solid boundaries and fills exterior ghost cells with H/HS from OMESH.
!> Thread-safe kernel version of NO_FLUX (velo.f90).
!> \param M Mesh data structure
!> \param DT Time step (s)

SUBROUTINE NO_FLUX_KERNEL(M,DT)

USE MESH_VARIABLES, ONLY: MESHES

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,OM_HP
REAL(EB) :: RFODT,H_OTHER,DUUDT,DVVDT,DWWDT,UN,DHFCT
INTEGER  :: IC2,IC1,N,I,J,K,IW,II,JJ,KK,IOR,N_INT_CELLS,IIO,JJO,KKO,NOM
TYPE(OBSTRUCTION_TYPE), POINTER :: OB
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

RFODT = RELAXATION_FACTOR/DT

IF (PREDICTOR) THEN
   HP => M%H
ELSE
   HP => M%HS
ENDIF

! Fill in exterior cells of mesh with values of HP from neighboring meshes

DO IW=1,M%N_EXTERNAL_WALL_CELLS
   EWC=>M%EXTERNAL_WALL(IW)
   NOM =EWC%NOM
   IF (NOM==0) CYCLE
   WC=>M%WALL(IW)
   IF (PREDICTOR) THEN
      OM_HP=>M%OMESH(NOM)%H
   ELSE
      OM_HP=>M%OMESH(NOM)%HS
   ENDIF
   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   II = BC%II
   JJ = BC%JJ
   KK = BC%KK
   H_OTHER = 0._EB
   DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
      DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
         DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
            H_OTHER = H_OTHER + OM_HP(IIO,JJO,KKO)
         ENDDO
      ENDDO
   ENDDO
   N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
   HP(II,JJ,KK) = H_OTHER/REAL(N_INT_CELLS,EB)
ENDDO

! Set FVX, FVY and FVZ to drive velocity components at solid boundaries within obstructions towards zero

OBST_LOOP: DO N=1,M%N_OBST

   OB=>M%OBSTRUCTION(N)

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1  ,OB%I2
            IC1 = M%CELL_INDEX(I,J,K)
            IC2 = M%CELL_INDEX(I+1,J,K)
            IF (M%CELL(IC1)%SOLID .AND. M%CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DUUDT = -RFODT*M%U(I,J,K)
               ELSE
                  DUUDT = -RFODT*(M%U(I,J,K)+M%US(I,J,K))
               ENDIF
               M%FVX(I,J,K) = -M%RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - DUUDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1  ,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = M%CELL_INDEX(I,J,K)
            IC2 = M%CELL_INDEX(I,J+1,K)
            IF (M%CELL(IC1)%SOLID .AND. M%CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DVVDT = -RFODT*M%V(I,J,K)
               ELSE
                  DVVDT = -RFODT*(M%V(I,J,K)+M%VS(I,J,K))
               ENDIF
               M%FVY(I,J,K) = -M%RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - DVVDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1  ,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = M%CELL_INDEX(I,J,K)
            IC2 = M%CELL_INDEX(I,J,K+1)
            IF (M%CELL(IC1)%SOLID .AND. M%CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DWWDT = -RFODT*M%W(I,J,K)
               ELSE
                  DWWDT = -RFODT*(M%W(I,J,K)+M%WS(I,J,K))
               ENDIF
               M%FVZ(I,J,K) = -M%RDZN(K)*(HP(I,J,K+1)-HP(I,J,K)) - DWWDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDDO OBST_LOOP

! Set FVX, FVY and FVZ to drive the normal velocity at solid boundaries towards the specified value

WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS

   WC => M%WALL(IW)

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==OPEN_BOUNDARY) CYCLE WALL_LOOP

   IF (IW<=M%N_EXTERNAL_WALL_CELLS) THEN
      NOM = M%EXTERNAL_WALL(IW)%NOM
   ELSE
      NOM = 0
   ENDIF

   IF (IW>M%N_EXTERNAL_WALL_CELLS .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. NOM==0) CYCLE WALL_LOOP

   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR

   DHFCT=1._EB
   SELECT CASE(PRES_FLAG)
      CASE(UGLMAT_FLAG,ULMAT_FLAG); DHFCT=0._EB
      CASE(GLMAT_FLAG); IF (IW<=M%N_EXTERNAL_WALL_CELLS) DHFCT=0._EB
   END SELECT

   IF (NOM/=0 .OR. WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY) THEN
      B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
      IF (PREDICTOR) THEN
         UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
      ELSE
         UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      ENDIF
      SELECT CASE(IOR)
         CASE( 1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-M%U(II,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(M%U(II,JJ,KK)+M%US(II,JJ,KK)) )
            ENDIF
            M%FVX(II,JJ,KK) = -M%RDXN(II)*(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT - DUUDT
         CASE(-1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-M%U(II-1,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(M%U(II-1,JJ,KK)+M%US(II-1,JJ,KK)) )
            ENDIF
            M%FVX(II-1,JJ,KK) = -M%RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT - DUUDT
         CASE( 2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-M%V(II,JJ,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(M%V(II,JJ,KK)+M%VS(II,JJ,KK)) )
            ENDIF
            M%FVY(II,JJ,KK) = -M%RDYN(JJ)*(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT - DVVDT
         CASE(-2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-M%V(II,JJ-1,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(M%V(II,JJ-1,KK)+M%VS(II,JJ-1,KK)) )
            ENDIF
            M%FVY(II,JJ-1,KK) = -M%RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT - DVVDT
         CASE( 3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-M%W(II,JJ,KK))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(M%W(II,JJ,KK)+M%WS(II,JJ,KK)) )
            ENDIF
            M%FVZ(II,JJ,KK) = -M%RDZN(KK)*(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT - DWWDT
         CASE(-3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-M%W(II,JJ,KK-1))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(M%W(II,JJ,KK-1)+M%WS(II,JJ,KK-1)) )
            ENDIF
            M%FVZ(II,JJ,KK-1) = -M%RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT - DWWDT
      END SELECT
   ENDIF

   IF (WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) THEN
      SELECT CASE(IOR)
         CASE( 1)
            M%FVX(II  ,JJ,KK) = 0._EB
         CASE(-1)
            M%FVX(II-1,JJ,KK) = 0._EB
         CASE( 2)
            M%FVY(II  ,JJ,KK) = 0._EB
         CASE(-2)
            M%FVY(II,JJ-1,KK) = 0._EB
         CASE( 3)
            M%FVZ(II  ,JJ,KK) = 0._EB
         CASE(-3)
            M%FVZ(II,JJ,KK-1) = 0._EB
      END SELECT
   ENDIF

ENDDO WALL_LOOP

END SUBROUTINE NO_FLUX_KERNEL


!> \brief Match velocity flux (FVX/FVY/FVZ) at interpolated mesh boundaries.
!> \details Thread-safe kernel version of MATCH_VELOCITY_FLUX (velo.f90).
!> Averages FVX/FVY/FVZ at exterior interpolated boundaries with values from neighboring meshes.
!> \param M Mesh data structure
!> \param NM Mesh index

SUBROUTINE MATCH_VELOCITY_FLUX_KERNEL(M,NM)

USE MESH_VARIABLES, ONLY: MESHES

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
REAL(EB) :: DA_OTHER,FVX_OTHER,FVY_OTHER,FVZ_OTHER
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

IF (NMESHES==1) RETURN
IF (SOLID_PHASE_ONLY) RETURN

! Loop over all external wall cells and match flux at interpolated boundaries

EXTERNAL_WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS

   WC=>M%WALL(IW)
   EWC=>M%EXTERNAL_WALL(IW)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR
   NOM = EWC%NOM
   OM => M%OMESH(NOM)
   M2 => MESHES(NOM)

   ! Determine the area of the interpolated cell face

   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   ! Determine the normal component of velocity flux from the other mesh and use it for average

   SELECT CASE(IOR)

      CASE( 1)
         FVX_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVX_OTHER = FVX_OTHER + OM%FVX(IIO,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVX(0,JJ,KK) = 0.5_EB*(M%FVX(0,JJ,KK) + FVX_OTHER)

      CASE(-1)
         FVX_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVX_OTHER = FVX_OTHER + OM%FVX(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVX(M%IBAR,JJ,KK) = 0.5_EB*(M%FVX(M%IBAR,JJ,KK) + FVX_OTHER)

      CASE( 2)
         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVY(II,0,KK) = 0.5_EB*(M%FVY(II,0,KK) + FVY_OTHER)

      CASE(-2)
         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVY(II,M%JBAR,KK) = 0.5_EB*(M%FVY(II,M%JBAR,KK) + FVY_OTHER)

      CASE( 3)
         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVZ(II,JJ,0) = 0.5_EB*(M%FVZ(II,JJ,0) + FVZ_OTHER)

      CASE(-3)
         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         M%FVZ(II,JJ,M%KBAR) = 0.5_EB*(M%FVZ(II,JJ,M%KBAR) + FVZ_OTHER)

   END SELECT

ENDDO EXTERNAL_WALL_LOOP

END SUBROUTINE MATCH_VELOCITY_FLUX_KERNEL


END MODULE VELO_KERNELS
