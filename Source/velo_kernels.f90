!> \brief Pure computation kernels extracted from VELO module.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE VELO_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE,BOUNDARY_PROP2_TYPE,SURFACE_TYPE,SURFACE, &
                 RAMPS_TYPE,RAMPS

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC BAROCLINIC_CORRECTION_KERNEL,VELOCITY_PREDICTOR_KERNEL,VELOCITY_CORRECTOR_KERNEL,VELOCITY_FLUX_KERNEL, &
       COMPUTE_VISCOSITY_KERNEL,CHECK_STABILITY_KERNEL

CONTAINS


!> \brief Compute the baroclinic torque correction terms.
!> \param M Mesh data structure
!> \param T Current simulation time (s)

SUBROUTINE BAROCLINIC_CORRECTION_KERNEL(M,T)

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

!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC)
DO K=0,M%KBP1
   DO J=0,M%JBP1
      DO I=0,M%IBP1
         P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-M%KRES(I,J,K))
         RRHO(I,J,K) = 1._EB/RHOP(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO

! Compute baroclinic term in the x momentum equation, p*d/dx(1/rho)

!$OMP DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%FVX_B(I,J,K) = -(P(I,J,K)*RHOP(I+1,J,K)+P(I+1,J,K)*RHOP(I,J,K))*(RRHO(I+1,J,K)-RRHO(I,J,K))*M%RDXN(I)/ &
                         (RHOP(I+1,J,K)+RHOP(I,J,K))
         M%FVX(I,J,K) = M%FVX(I,J,K) + M%FVX_B(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

! Compute baroclinic term in the y momentum equation, p*d/dy(1/rho)

IF (.NOT.TWO_D) THEN
!$OMP DO SCHEDULE(STATIC)
   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            M%FVY_B(I,J,K) = -(P(I,J,K)*RHOP(I,J+1,K)+P(I,J+1,K)*RHOP(I,J,K))*(RRHO(I,J+1,K)-RRHO(I,J,K))*M%RDYN(J)/ &
                            (RHOP(I,J+1,K)+RHOP(I,J,K))
            M%FVY(I,J,K) = M%FVY(I,J,K) + M%FVY_B(I,J,K)
         ENDDO
      ENDDO
   ENDDO
!$OMP END DO NOWAIT
ENDIF

! Compute baroclinic term in the z momentum equation, p*d/dz(1/rho)

!$OMP DO SCHEDULE(STATIC)
DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%FVZ_B(I,J,K) = -(P(I,J,K)*RHOP(I,J,K+1)+P(I,J,K+1)*RHOP(I,J,K))*(RRHO(I,J,K+1)-RRHO(I,J,K))*M%RDZN(K)/ &
                         (RHOP(I,J,K+1)+RHOP(I,J,K))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + M%FVZ_B(I,J,K)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

M%BAROCLINIC_TERMS_ATTACHED = .TRUE.

END SUBROUTINE BAROCLINIC_CORRECTION_KERNEL


!> \brief Predict the velocity components at the next time step.
!> \details Core OMP loops that compute US/VS/WS from U/V/W, FVX/FVY/FVZ, and H.
!> \param M Mesh data structure
!> \param DT Time step (s)

SUBROUTINE VELOCITY_PREDICTOR_KERNEL(M,DT)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER :: I,J,K

IF (FREEZE_VELOCITY) THEN
   M%US = M%U
   M%VS = M%V
   M%WS = M%W
   RETURN
ENDIF

!$OMP PARALLEL PRIVATE(I,J,K)

!$OMP DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%US(I,J,K) = M%U(I,J,K) - DT*( M%FVX(I,J,K) + M%RDXN(I)*(M%H(I+1,J,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%VS(I,J,K) = M%V(I,J,K) - DT*( M%FVY(I,J,K) + M%RDYN(J)*(M%H(I,J+1,K)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%WS(I,J,K) = M%W(I,J,K) - DT*( M%FVZ(I,J,K) + M%RDZN(K)*(M%H(I,J,K+1)-M%H(I,J,K)) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END SUBROUTINE VELOCITY_PREDICTOR_KERNEL


!> \brief Correct the velocity components at the next time step.
!> \details Core OMP loops that compute U/V/W from U/V/W, US/VS/WS, FVX/FVY/FVZ, and HS.
!> \param M Mesh data structure
!> \param DT Time step (s)

SUBROUTINE VELOCITY_CORRECTOR_KERNEL(M,DT)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER :: I,J,K

IF (FREEZE_VELOCITY) THEN
   M%U = M%US
   M%V = M%VS
   M%W = M%WS
   RETURN
ENDIF

!$OMP PARALLEL PRIVATE(I,J,K)

!$OMP DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         M%U(I,J,K) = 0.5_EB*( M%U(I,J,K) + M%US(I,J,K) - DT*(M%FVX(I,J,K) + M%RDXN(I)*(M%HS(I+1,J,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         M%V(I,J,K) = 0.5_EB*( M%V(I,J,K) + M%VS(I,J,K) - DT*(M%FVY(I,J,K) + M%RDYN(J)*(M%HS(I,J+1,K)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         M%W(I,J,K) = 0.5_EB*( M%W(I,J,K) + M%WS(I,J,K) - DT*(M%FVZ(I,J,K) + M%RDZN(K)*(M%HS(I,J,K+1)-M%HS(I,J,K))) )
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END SUBROUTINE VELOCITY_CORRECTOR_KERNEL


!> \brief Compute the velocity flux terms (vorticity, stress tensor, momentum RHS).
!> \param M Mesh data structure
!> \param T Current simulation time (s)
!> \param DT Time step (s)
!> \param NM Mesh number (needed for external calls)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables
!> \param GX Gravity component arrays (output, for use by CC_IBM in wrapper)
!> \param GY Gravity component arrays (output, for use by CC_IBM in wrapper)
!> \param GZ Gravity component arrays (output, for use by CC_IBM in wrapper)

SUBROUTINE VELOCITY_FLUX_KERNEL(M,T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES,GX,GY,GZ)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY: COMPUTE_WIND_COMPONENTS
USE CC_SCALARS, ONLY : ROTATED_CUBE_VELOCITY_FLUX

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

!$OMP PARALLEL DO PRIVATE(DUDY,DVDX,DUDZ,DWDX,DVDZ,DWDY,MUX,MUY,MUZ) SCHEDULE(STATIC)
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
!$OMP END PARALLEL DO

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

!$OMP PARALLEL PRIVATE(WP,WM,VP,VM,UP,UM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXZP,TXZM,TXYP,TXYM,TYZP,TYZM, &
!$OMP& IC,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,RRHO,DUDX,DVDY,DWDZ,VTRM)

!$OMP DO SCHEDULE(STATIC) PRIVATE(WOMY, VOMZ, TXXP, TXXM, DTXXDX, DTXYDY, DTXZDZ)
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
!$OMP END DO NOWAIT

! Compute y-direction flux term FVY

!$OMP DO SCHEDULE(STATIC) PRIVATE(WOMX, UOMZ, TYYP, TYYM, DTXYDX, DTYYDY, DTYZDZ)
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
!$OMP END DO NOWAIT

! Compute z-direction flux term FVZ

!$OMP DO SCHEDULE(STATIC) PRIVATE(UOMY, VOMX, TZZP, TZZM, DTXZDX, DTYZDY, DTZZDZ)
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
!$OMP END DO NOWAIT

!$OMP END PARALLEL

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

   !$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=0,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
            M%FVX(I,J,K) = M%FVX(I,J,K) - RRHO*FVEC(1)*TIME_RAMP_FACTOR*SIN_THETA
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

IF (ABS(FVEC(2))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVY_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVY_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   !$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
            M%FVY(I,J,K) = M%FVY(I,J,K) - RRHO*FVEC(2)*TIME_RAMP_FACTOR*COS_THETA
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
ENDIF

IF (ABS(FVEC(3))>TWENTY_EPSILON_EB) THEN
   IF (I_RAMP_FVZ_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_FVZ_T)
   ELSEIF (I_RAMP_PGF_T>0) THEN
      TIME_RAMP_FACTOR = EVALUATE_RAMP(T,I_RAMP_PGF_T)
   ELSE
      TIME_RAMP_FACTOR = 1._EB
   ENDIF

   !$OMP PARALLEL DO PRIVATE(RRHO) SCHEDULE(STATIC)
   DO K=0,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RRHO = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
            M%FVZ(I,J,K) = M%FVZ(I,J,K) - RRHO*FVEC(3)*TIME_RAMP_FACTOR
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
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

!$OMP PARALLEL DO SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

DO IW=1,M%N_EXTERNAL_WALL_CELLS
   WC=>M%WALL(IW)
   BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
   UP(BC%II,BC%JJ,BC%KK) = M%U_GHOST(IW)
   VP(BC%II,BC%JJ,BC%KK) = M%V_GHOST(IW)
   WP(BC%II,BC%JJ,BC%KK) = M%W_GHOST(IW)
ENDDO

! x momentum

!$OMP PARALLEL DO PRIVATE(VBAR,WBAR) SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=1,M%JBAR
      DO I=0,M%IBAR
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I+1,J,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I+1,J,K))
         M%FVX(I,J,K) = M%FVX(I,J,K) + 2._EB*(OVEC(2)*WBAR-OVEC(3)*VBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! y momentum

!$OMP PARALLEL DO PRIVATE(UBAR,WBAR) SCHEDULE(STATIC)
DO K=1,M%KBAR
   DO J=0,M%JBAR
      DO I=1,M%IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J+1,K))
         WBAR = 0.5_EB*(WP(I,J,K)+WP(I,J+1,K))
         M%FVY(I,J,K) = M%FVY(I,J,K) + 2._EB*(OVEC(3)*UBAR - OVEC(1)*WBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

! z momentum

!$OMP PARALLEL DO PRIVATE(UBAR,VBAR) SCHEDULE(STATIC)
DO K=0,M%KBAR
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         UBAR = 0.5_EB*(UP(I,J,K)+UP(I,J,K+1))
         VBAR = 0.5_EB*(VP(I,J,K)+VP(I,J,K+1))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + 2._EB*(OVEC(1)*VBAR - OVEC(2)*UBAR)
      ENDDO
   ENDDO
ENDDO
!$OMP END PARALLEL DO

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


!> \brief Compute the turbulent viscosity.
!> \param M Mesh data structure
!> \param NM Mesh number (needed for external calls)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag for estimated (starred) variables

SUBROUTINE COMPUTE_VISCOSITY_KERNEL(M,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY,GET_POTENTIAL_TEMPERATURE,GET_CONDUCTIVITY,GET_SPECIFIC_HEAT
USE TURB_KERNELS, ONLY: WALE_VISCOSITY,FILL_EDGES_KERNEL,TEST_FILTER_KERNEL,VARDEN_DYNSMAG_KERNEL
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE CC_VELOCITY, ONLY : CC_COMPUTE_KRES,CC_COMPUTE_VISCOSITY,CUTFACE_VELOCITIES

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

   !$OMP PARALLEL PRIVATE(ZZ_GET)
   ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
   !$OMP DO SCHEDULE(STATIC)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%CELL(M%CELL_INDEX(I,J,K))%SOLID) CYCLE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_VISCOSITY(ZZ_GET,M%MU_DNS(I,J,K),M%TMP(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO
   DEALLOCATE(ZZ_GET)
   !$OMP END PARALLEL

ENDIF

IF (CC_IBM) CALL CUTFACE_VELOCITIES(NM,UU,VV,WW,CUTFACES=.TRUE.)

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

      !$OMP PARALLEL

      !$OMP DO SCHEDULE(STATIC) PRIVATE(I,J,K)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
               VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
               WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

      ! fill mesh boundary ghost cells

      !$OMP DO SCHEDULE(STATIC) PRIVATE(IW,WC,BC)
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
      !$OMP END DO

      !$OMP END PARALLEL

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

      !$OMP PARALLEL DO PRIVATE(DELTA, KSGS, NU_EDDY) SCHEDULE(STATIC)
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
      !$OMP END PARALLEL DO

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

IF (CC_IBM) CALL CC_COMPUTE_KRES(APPLY_TO_ESTIMATED_VARIABLES,NM)

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
   CALL CC_COMPUTE_VISCOSITY(0._EB,NM)
   CALL CUTFACE_VELOCITIES(NM,UU,VV,WW,CUTFACES=.FALSE.)
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


!> \brief Check the CFL and Von Neumann stability criteria.
!> \param M Mesh data structure
!> \param DT Current time step (s)
!> \param DT_NEW_MESH New time step for this mesh (output)
!> \param T Current simulation time (s)
!> \param NM Mesh number (needed for CHANGE_TIME_STEP_INDEX and diagnostic writes)

SUBROUTINE CHECK_STABILITY_KERNEL(M,DT,DT_NEW_MESH,T,NM)

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

!$OMP PARALLEL PRIVATE(ICFL_TMP, JCFL_TMP, KCFL_TMP, UODX, VODY, WODZ, UVW, UVWMAX_TMP) SHARED(UVWMAX)
UVWMAX_TMP = 0._EB
!$OMP DO SCHEDULE(STATIC)
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
!$OMP END DO
!$OMP CRITICAL
IF(UVWMAX_TMP>UVWMAX) THEN
   UVWMAX = UVWMAX_TMP
   M%ICFL = ICFL_TMP
   M%JCFL = JCFL_TMP
   M%KCFL = KCFL_TMP
ENDIF
!$OMP END CRITICAL
!$OMP END PARALLEL

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


END MODULE VELO_KERNELS
