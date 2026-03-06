!> \brief Collection of velocity routines.
!> Computes the velocity flux terms, baroclinic torque correction terms, and performs the CFL check.

MODULE VELO

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC VELOCITY_PREDICTOR,VELOCITY_CORRECTOR,NO_FLUX,BAROCLINIC_CORRECTION,MATCH_VELOCITY,MATCH_VELOCITY_FLUX,&
       VELOCITY_BC,COMPUTE_VISCOSITY,VISCOSITY_BC,VELOCITY_FLUX,VELOCITY_FLUX_CYLINDRICAL,&
       CHECK_STABILITY


CONTAINS


!> \brief Compute the viscosity of the gas.
!> \callergraph
!> \callgraph
!> \param NM Mesh number.
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating \f$\mu(T,\mathbf{u})\f$ or \f$\mu(T^*,\mathbf{u}^*)\f$

SUBROUTINE COMPUTE_VISCOSITY(NM,APPLY_TO_ESTIMATED_VARIABLES)

USE VELO_KERNELS, ONLY: COMPUTE_VISCOSITY_KERNEL

INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: T_NOW

T_NOW = CURRENT_TIME()

CALL COMPUTE_VISCOSITY_KERNEL(MESHES(NM),NM,APPLY_TO_ESTIMATED_VARIABLES)

T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

END SUBROUTINE COMPUTE_VISCOSITY


!> \brief Compute boundary values of the viscosity, MU.
!> \param NM Mesh number.
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating \f$\mu(T,\mathbf{u})\f$ or \f$\mu(T^*,\mathbf{u}^*)\f$
!> \callergraph
!> \callgraph

SUBROUTINE VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES)

! Specify ghost cell values of the viscosity array MU

INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: MU_OTHER,DP_OTHER,KRES_OTHER,T_NOW
INTEGER :: II,JJ,KK,IW,IIO,JJO,KKO,NOM,N_INT_CELLS
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

T_NOW = CURRENT_TIME()

CALL POINT_TO_MESH(NM)

! Mirror viscosity into solids and exterior boundary cells

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
   WC =>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (EWC%NOM==0) CYCLE WALL_LOOP
   BC => BOUNDARY_COORD(WC%BC_INDEX)
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
            MU_OTHER = MU_OTHER + OMESH(NOM)%MU(IIO,JJO,KKO)
            KRES_OTHER = KRES_OTHER + OMESH(NOM)%KRES(IIO,JJO,KKO)
            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
               DP_OTHER = DP_OTHER + OMESH(NOM)%DS(IIO,JJO,KKO)
            ELSE
               DP_OTHER = DP_OTHER + OMESH(NOM)%D(IIO,JJO,KKO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
   MU_OTHER = MU_OTHER/REAL(N_INT_CELLS,EB)
   KRES_OTHER = KRES_OTHER/REAL(N_INT_CELLS,EB)
   DP_OTHER = DP_OTHER/REAL(N_INT_CELLS,EB)
   MU(II,JJ,KK) = MU_OTHER
   KRES(II,JJ,KK) = KRES_OTHER
   IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
      DS(II,JJ,KK) = DP_OTHER
   ELSE
      D(II,JJ,KK) = DP_OTHER
   ENDIF
ENDDO WALL_LOOP

T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

END SUBROUTINE VISCOSITY_BC


!> \brief Compute convective and diffusive terms of the momentum equations
!> \param T Current time (s)
!> \param DT Current time step (s)
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating whether to use estimated values of variables

SUBROUTINE VELOCITY_FLUX(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE CC_VELOCITY_KERNELS, ONLY : CC_VELOCITY_FLUX,CUTFACE_VELOCITIES
USE CC_VELOCITY, ONLY : CC_VELOCITY_BC
USE VELO_KERNELS, ONLY: VELOCITY_FLUX_KERNEL

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: T_NOW,GX(0:IBAR_MAX),GY(0:IBAR_MAX),GZ(0:IBAR_MAX)
INTEGER :: I
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP

T_NOW=CURRENT_TIME()

CALL POINT_TO_MESH(NM)

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
ELSE
   UU => U
   VV => V
   WW => W
   RHOP => RHO
ENDIF

! Define velocities on gas cut-faces underlaying Cartesian faces.

IF (CC_IBM) THEN
   T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW
   CALL CC_VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES,DO_IBEDGES=.FALSE.)
   T_NOW=CURRENT_TIME()
   CALL CUTFACE_VELOCITIES(MESHES(NM),UU,VV,WW, &
      CUTFACES=.TRUE.)
   T_USED(14) = T_USED(14) + CURRENT_TIME() - T_NOW
   T_NOW=CURRENT_TIME()
ENDIF

CALL VELOCITY_FLUX_KERNEL(MESHES(NM),T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES,GX,GY,GZ)

! Restore previous substep velocities to gas cut-faces underlaying Cartesian faces.

IF (CC_IBM) THEN
   T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW
   T_NOW=CURRENT_TIME()
   CALL CUTFACE_VELOCITIES(MESHES(NM),UU,VV,WW, &
      CUTFACES=.FALSE.)
   T_USED(14) = T_USED(14) + CURRENT_TIME() - T_NOW
   CALL CC_VELOCITY_FLUX(MESHES(NM),DT, &
      APPLY_TO_ESTIMATED_VARIABLES,RHOP, &
      CORRECT_GRAV=.TRUE.,GX=GX,GY=GY,GZ=GZ)
   T_NOW=CURRENT_TIME()
ENDIF

T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

END SUBROUTINE VELOCITY_FLUX


!> \brief Compute convective and diffusive terms of the momentum equations in 2-D cylindrical coordinates
!> \param T Current time (s)
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating whether to use estimated values of variables

SUBROUTINE VELOCITY_FLUX_CYLINDRICAL(T,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
REAL(EB) :: T,DMUDX
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
INTEGER :: I0
INTEGER, INTENT(IN) :: NM
REAL(EB) :: MUY,UP,UM,WP,WM,VTRM,DTXZDZ,DTXZDX,DUDX,DWDZ,DUDZ,DWDX,WOMY,UOMY,OMYP,OMYM,TXZP,TXZM, &
            AH,RRHO,GX,GZ,TXXP,TXXM,TZZP,TZZM,DTXXDX,DTZZDZ,T_NOW
INTEGER :: I,J,K,IEYP,IEYM,IC
REAL(EB), POINTER, DIMENSION(:,:,:) :: TXZ,OMY,UU,WW,RHOP,DP

T_NOW = CURRENT_TIME()

CALL POINT_TO_MESH(NM)

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => US
   WW => WS
   DP => DS
   RHOP => RHOS
ELSE
   UU => U
   WW => W
   DP => D
   RHOP => RHO
ENDIF

TXZ => WORK2
OMY => WORK5

! Compute vorticity and stress tensor components

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         OMY(I,J,K) = DUDZ - DWDX
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
      ENDDO
   ENDDO
ENDDO

! Compute gravity components

GX  = 0._EB
GZ  = EVALUATE_RAMP(T,I_RAMP_GZ)*GVEC(3)

! Compute r-direction flux term FVX

IF (ABS(XS)<=TWENTY_EPSILON_EB) THEN
   I0 = 1
ELSE
   I0 = 0
ENDIF

J = 1

DO K= 1,KBAR
   DO I=I0,IBAR
      WP    = WW(I,J,K)   + WW(I+1,J,K)
      WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
      OMYP  = OMY(I,J,K)
      OMYM  = OMY(I,J,K-1)
      TXZP  = TXZ(I,J,K)
      TXZM  = TXZ(I,J,K-1)
      IC    = CELL_INDEX(I,J,K)
      IEYP  = CELL(IC)%EDGE_INDEX(8)
      IEYM  = CELL(IC)%EDGE_INDEX(6)
      IF (EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
         OMYP = EDGE(IEYP)%OMEGA(-1)
         TXZP = EDGE(IEYP)%TAU(-1)
      ENDIF
      IF (EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
         OMYM = EDGE(IEYM)%OMEGA( 1)
         TXZM = EDGE(IEYM)%TAU( 1)
      ENDIF
      WOMY  = WP*OMYP + WM*OMYM
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
      AH    = RHO_0(K)*RRHO - 1._EB
      DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
      TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*DWDZ )
      DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
      TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K) -2._EB*DWDZ )
      DTXXDX= RDXN(I)*(TXXP-TXXM)
      DTXZDZ= RDZ(K) *(TXZP-TXZM)
      DMUDX = (MU(I+1,J,K)-MU(I,J,K))*RDXN(I)
      VTRM  = RRHO*( DTXXDX + DTXZDZ - 2._EB*UU(I,J,K)*DMUDX/R(I) )
      FVX(I,J,K) = 0.25_EB*WOMY + GX*AH - VTRM
   ENDDO
ENDDO

! Compute z-direction flux term FVZ

DO K=0,KBAR
   DO I=1,IBAR
      UP    = UU(I,J,K)   + UU(I,J,K+1)
      UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
      OMYP  = OMY(I,J,K)
      OMYM  = OMY(I-1,J,K)
      TXZP  = TXZ(I,J,K)
      TXZM  = TXZ(I-1,J,K)
      IC    = CELL_INDEX(I,J,K)
      IEYP  = CELL(IC)%EDGE_INDEX(8)
      IEYM  = CELL(IC)%EDGE_INDEX(7)
      IF (EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
         OMYP = EDGE(IEYP)%OMEGA(-2)
         TXZP = EDGE(IEYP)%TAU(-2)
      ENDIF
      IF (EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
         OMYM = EDGE(IEYM)%OMEGA( 2)
         TXZM = EDGE(IEYM)%TAU( 2)
      ENDIF
      UOMY  = UP*OMYP + UM*OMYM
      RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
      AH    = 0.5_EB*(RHO_0(K)+RHO_0(K+1))*RRHO - 1._EB
      DUDX  = (R(I)*UU(I,J,K+1)-R(I-1)*UU(I-1,J,K+1))*RDX(I)*RRN(I)
      TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*DUDX )
      DUDX  = (R(I)*UU(I,J,K)-R(I-1)*UU(I-1,J,K))*RDX(I)*RRN(I)
      TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*DUDX )
      DTXZDX= RDX(I) *(R(I)*TXZP-R(I-1)*TXZM)*RRN(I)
      DTZZDZ= RDZN(K)*(     TZZP       -TZZM)
      VTRM  = RRHO*(DTXZDX + DTZZDZ)
      FVZ(I,J,K) = -0.25_EB*UOMY + GZ*AH - VTRM
   ENDDO
ENDDO

T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

END SUBROUTINE VELOCITY_FLUX_CYLINDRICAL


!> \brief Set momentum fluxes inside and on the surface of solid obstructions to maintain user-specified flux
!> \param DT Time step (s)
!> \param NM Mesh number

SUBROUTINE NO_FLUX(DT,NM)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,OM_HP
REAL(EB) :: RFODT,H_OTHER,DUUDT,DVVDT,DWWDT,UN,T_NOW,DHFCT
INTEGER  :: IC2,IC1,N,I,J,K,IW,II,JJ,KK,IOR,N_INT_CELLS,IIO,JJO,KKO,NOM
TYPE(OBSTRUCTION_TYPE), POINTER :: OB
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

T_NOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

RFODT = RELAXATION_FACTOR/DT

IF (PREDICTOR) THEN
   HP => H
ELSE
   HP => HS
ENDIF

! Fill in exterior cells of mesh NM with values of HP from mesh NOM

DO IW=1,N_EXTERNAL_WALL_CELLS
   EWC=>EXTERNAL_WALL(IW)
   NOM =EWC%NOM
   IF (NOM==0) CYCLE
   WC=>WALL(IW)
   IF (PREDICTOR) THEN
      OM_HP=>OMESH(NOM)%H
   ELSE
      OM_HP=>OMESH(NOM)%HS
   ENDIF
   BC => BOUNDARY_COORD(WC%BC_INDEX)
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
   HP(II,JJ,KK)  = H_OTHER/REAL(N_INT_CELLS,EB)
ENDDO

! Set FVX, FVY and FVZ to drive velocity components at solid boundaries within obstructions towards zero

OBST_LOOP: DO N=1,N_OBST

   OB=>OBSTRUCTION(N)

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1  ,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I+1,J,K)
            IF (CELL(IC1)%SOLID .AND. CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DUUDT = -RFODT*U(I,J,K)
               ELSE
                  DUUDT = -RFODT*(U(I,J,K)+US(I,J,K))
               ENDIF
               FVX(I,J,K) = -RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - DUUDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1+1,OB%K2
      DO J=OB%J1  ,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J+1,K)
            IF (CELL(IC1)%SOLID .AND. CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DVVDT = -RFODT*V(I,J,K)
               ELSE
                  DVVDT = -RFODT*(V(I,J,K)+VS(I,J,K))
               ENDIF
               FVY(I,J,K) = -RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - DVVDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO K=OB%K1  ,OB%K2
      DO J=OB%J1+1,OB%J2
         DO I=OB%I1+1,OB%I2
            IC1 = CELL_INDEX(I,J,K)
            IC2 = CELL_INDEX(I,J,K+1)
            IF (CELL(IC1)%SOLID .AND. CELL(IC2)%SOLID) THEN
               IF (PREDICTOR) THEN
                  DWWDT = -RFODT*W(I,J,K)
               ELSE
                  DWWDT = -RFODT*(W(I,J,K)+WS(I,J,K))
               ENDIF
               FVZ(I,J,K) = -RDZN(K)*(HP(I,J,K+1)-HP(I,J,K)) - DWWDT
            ENDIF
         ENDDO
      ENDDO
   ENDDO

ENDDO OBST_LOOP

! Set FVX, FVY and FVZ to drive the normal velocity at solid boundaries towards the specified value (U_NORMAL or U_NORMAL_S)

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS

   WC => WALL(IW)

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. WC%BOUNDARY_TYPE==OPEN_BOUNDARY) CYCLE WALL_LOOP

   IF (IW<=N_EXTERNAL_WALL_CELLS) THEN
      NOM = EXTERNAL_WALL(IW)%NOM
   ELSE
      NOM = 0
   ENDIF

   IF (IW>N_EXTERNAL_WALL_CELLS .AND. WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. NOM==0) CYCLE WALL_LOOP

   BC => BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR

   DHFCT=1._EB
   SELECT CASE(PRES_FLAG)
      CASE(UGLMAT_FLAG,ULMAT_FLAG); DHFCT=0._EB
      CASE(GLMAT_FLAG); IF (IW<=N_EXTERNAL_WALL_CELLS) DHFCT=0._EB
   END SELECT

   IF (NOM/=0 .OR. WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY) THEN
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IF (PREDICTOR) THEN
         UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
      ELSE
         UN = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      ENDIF
      SELECT CASE(IOR)
         CASE( 1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-U(II,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)) )
            ENDIF
            FVX(II,JJ,KK) = -RDXN(II)*(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT - DUUDT
         CASE(-1)
            IF (PREDICTOR) THEN
               DUUDT = RFODT*(UN-U(II-1,JJ,KK))
            ELSE
               DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK)) )
            ENDIF
            FVX(II-1,JJ,KK) = -RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT - DUUDT
         CASE( 2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-V(II,JJ,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)) )
            ENDIF
            FVY(II,JJ,KK) = -RDYN(JJ)*(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT - DVVDT
         CASE(-2)
            IF (PREDICTOR) THEN
               DVVDT = RFODT*(UN-V(II,JJ-1,KK))
            ELSE
               DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK)) )
            ENDIF
            FVY(II,JJ-1,KK) = -RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT - DVVDT
         CASE( 3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-W(II,JJ,KK))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)) )
            ENDIF
            FVZ(II,JJ,KK) = -RDZN(KK)*(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT - DWWDT
         CASE(-3)
            IF (PREDICTOR) THEN
               DWWDT = RFODT*(UN-W(II,JJ,KK-1))
            ELSE
               DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1)) )
            ENDIF
            FVZ(II,JJ,KK-1) = -RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT - DWWDT
      END SELECT
   ENDIF

   IF (WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) THEN
      SELECT CASE(IOR)
         CASE( 1)
            FVX(II  ,JJ,KK) = 0._EB
         CASE(-1)
            FVX(II-1,JJ,KK) = 0._EB
         CASE( 2)
            FVY(II  ,JJ,KK) = 0._EB
         CASE(-2)
            FVY(II,JJ-1,KK) = 0._EB
         CASE( 3)
            FVZ(II  ,JJ,KK) = 0._EB
         CASE(-3)
            FVZ(II,JJ,KK-1) = 0._EB
      END SELECT
   ENDIF

ENDDO WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

END SUBROUTINE NO_FLUX


!> \brief Estimate the velocity components at the next time step
!> \param T Current time (s)
!> \param DT Time step (s)
!> \param DT_NEW New time step (if necessary)
!> \param NM Mesh number

SUBROUTINE VELOCITY_PREDICTOR(T,DT,DT_NEW,NM)

USE TURBULENCE, ONLY: COMPRESSION_WAVE
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
USE CC_VELOCITY, ONLY : CC_PROJECT_VELOCITY
USE VELO_KERNELS, ONLY: VELOCITY_PREDICTOR_KERNEL

REAL(EB) :: T_NOW,XHAT,ZHAT
INTEGER  :: I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: DT_NEW(NMESHES)

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==4) THEN
   CALL COMPRESSION_WAVE(NM,T,4)
   CALL CHECK_STABILITY(DT,DT_NEW,T,NM)
   RETURN
ENDIF

T_NOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

CALL VELOCITY_PREDICTOR_KERNEL(MESHES(NM),DT)

IF (.NOT.FREEZE_VELOCITY) THEN
   IF (CC_IBM) THEN
      T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW
      CALL CC_PROJECT_VELOCITY(NM,DT,STORE_FLG=.FALSE.)
      T_NOW=CURRENT_TIME()
   ENDIF
   SELECT CASE(PRES_FLAG)
      CASE(GLMAT_FLAG,UGLMAT_FLAG,ULMAT_FLAG); CALL WALL_VELOCITY_NO_GRADH(DT,.FALSE.)
   END SELECT
ENDIF

! Manufactured solution (debug)

IF (PERIODIC_TEST==7 .AND. .FALSE.) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            XHAT =  X(I) - UF_MMS*(T)
            ZHAT = ZC(K) - WF_MMS*(T)
            US(I,J,K) = VD2D_MMS_U(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            XHAT = XC(I) - UF_MMS*(T)
            ZHAT =  Z(K) - WF_MMS*(T)
            WS(I,J,K) = VD2D_MMS_V(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
ENDIF

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

! Check the stability criteria, and if the time step is too small, send back a signal to kill the job

CALL CHECK_STABILITY(DT,DT_NEW,T,NM)

IF (DT_NEW(NM)<DT_INITIAL*LIMITING_DT_RATIO .AND. (T+DT_NEW(NM)<(T_END-TWENTY_EPSILON_EB))) STOP_STATUS = INSTABILITY_STOP


END SUBROUTINE VELOCITY_PREDICTOR


!> \brief Correct the velocity components at the next time step
!> \param T Current time (s)
!> \param DT Time step (s)
!> \param NM Mesh number

SUBROUTINE VELOCITY_CORRECTOR(T,DT,NM)

USE TURBULENCE, ONLY: COMPRESSION_WAVE
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
USE CC_VELOCITY, ONLY : CC_PROJECT_VELOCITY
USE VELO_KERNELS, ONLY: VELOCITY_CORRECTOR_KERNEL

REAL(EB) :: T_NOW,XHAT,ZHAT
INTEGER  :: I,J,K
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==4) THEN
   CALL COMPRESSION_WAVE(NM,T,4)
   RETURN
ENDIF

T_NOW=CURRENT_TIME()
CALL POINT_TO_MESH(NM)

IF (.NOT.FREEZE_VELOCITY) THEN
   IF (CC_IBM) THEN
      T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW
      CALL CC_PROJECT_VELOCITY(NM,DT,.TRUE.)
      T_NOW=CURRENT_TIME()
   ENDIF
   SELECT CASE(PRES_FLAG)
      CASE(GLMAT_FLAG,UGLMAT_FLAG,ULMAT_FLAG)
         CALL WALL_VELOCITY_NO_GRADH(DT,.TRUE.)                    ! Store U velocities on OBST surfaces.
   END SELECT
ENDIF

CALL VELOCITY_CORRECTOR_KERNEL(MESHES(NM),DT)

IF (.NOT.FREEZE_VELOCITY) THEN
   IF (CC_IBM) THEN
      T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW
      CALL CC_PROJECT_VELOCITY(NM,DT,.FALSE.)
      T_NOW=CURRENT_TIME()
   ENDIF
   SELECT CASE(PRES_FLAG)
      CASE(GLMAT_FLAG,UGLMAT_FLAG,ULMAT_FLAG)
         CALL WALL_VELOCITY_NO_GRADH(DT,.FALSE.)
   END SELECT
ENDIF

! Manufactured solution (debug)

IF (PERIODIC_TEST==7 .AND. .FALSE.) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=0,IBAR
            XHAT =  X(I) - UF_MMS*T
            ZHAT = ZC(K) - WF_MMS*T
            U(I,J,K) = VD2D_MMS_U(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
   DO K=0,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            XHAT = XC(I) - UF_MMS*T
            ZHAT =  Z(K) - WF_MMS*T
            W(I,J,K) = VD2D_MMS_V(XHAT,ZHAT,T)
         ENDDO
      ENDDO
   ENDDO
ENDIF

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW
END SUBROUTINE VELOCITY_CORRECTOR


!> \brief Assert tangential velocity boundary conditions
!> \param T Current time (s)
!> \param NM Mesh number
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating that estimated (starred) variables are to be used

SUBROUTINE VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE TURB_KERNELS, ONLY: WALL_MODEL
USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT
USE CC_VELOCITY, ONLY : CC_VELOCITY_BC,GET_OPENBC_TANGENTIAL_CUTFACE_VEL

REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: MUA,TSI,WGT,T_NOW,RAMP_T,OMW,MU_WALL,RHO_WALL,SLIP_COEF,VEL_T, &
            UUP(2),UUM(2),DXX(2),MU_DUIDXJ(-2:2),DUIDXJ(-2:2),PROFILE_FACTOR,VEL_GAS,VEL_GHOST, &
            MU_DUIDXJ_USE(2),DUIDXJ_USE(2),VEL_EDDY,U_TAU,Y_PLUS,U_NORM, &
            DRAG_FACTOR,HT_SCALE_FACTOR,VEG_HT,VEL_N,UN_OTHER
INTEGER :: NOM(2),IIO(2),JJO(2),KKO(2),IE,IIOO,JJOO,KKOO,II,JJ,KK,IEC,IOR,IW,IWM,IWP,ICMM,ICMP,ICPM,ICPP,ICD,ICDO,IVL,I_SGN, &
           VELOCITY_BC_INDEX,IIGM,JJGM,KKGM,IIGP,JJGP,KKGP,SURF_INDEXM,SURF_INDEXP,ITMP,ICD_SGN,ICDO_SGN, &
           BOUNDARY_TYPE_M,BOUNDARY_TYPE_P,IS,IS2,IWPI,IWMI,VENT_INDEX,N_INT_CELLS
LOGICAL :: ALTERED_GRADIENT(-2:2),SYNTHETIC_EDDY_METHOD,HVAC_TANGENTIAL,INTERPOLATED_EDGE,&
           UPWIND_BOUNDARY,INFLOW_BOUNDARY
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP,VEL_OTHER
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC,WCM,WCP,WCX
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC
REAL(EB), POINTER, DIMENSION(:,:,:) :: OM_UU,OM_VV,OM_WW
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: WCM_B1,WCP_B1,WCX_B1
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (EDGE_TYPE), POINTER :: ED
TYPE(SURFACE_TYPE), POINTER :: SF

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==12) RETURN
IF (PERIODIC_TEST==13) RETURN

T_NOW = CURRENT_TIME()

! Assign local names to variables

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => US
   VV => VS
   WW => WS
   RHOP => RHOS
   ZZP => ZZS
ELSE
   UU => U
   VV => V
   WW => W
   RHOP => RHO
   ZZP => ZZ
ENDIF

! Transfer from neighboring mesh the normal component of velocity that is one grid cell beyond external boundary

WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS
   WC =>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (EWC%NOM==0) CYCLE WALL_LOOP
   IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
      OM_UU => OMESH(EWC%NOM)%US
      OM_VV => OMESH(EWC%NOM)%VS
      OM_WW => OMESH(EWC%NOM)%WS
   ELSE
      OM_UU => OMESH(EWC%NOM)%U
      OM_VV => OMESH(EWC%NOM)%V
      OM_WW => OMESH(EWC%NOM)%W
   ENDIF
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   UN_OTHER = 0._EB
   DO KKOO=EWC%KKO_MIN,EWC%KKO_MAX
      DO JJOO=EWC%JJO_MIN,EWC%JJO_MAX
         DO IIOO=EWC%IIO_MIN,EWC%IIO_MAX
            SELECT CASE(BC%IOR)
               CASE(-1) ; UN_OTHER = UN_OTHER + OM_UU(IIOO  ,JJOO  ,KKOO  )
               CASE( 1) ; UN_OTHER = UN_OTHER + OM_UU(IIOO-1,JJOO  ,KKOO  )
               CASE(-2) ; UN_OTHER = UN_OTHER + OM_VV(IIOO  ,JJOO  ,KKOO  )
               CASE( 2) ; UN_OTHER = UN_OTHER + OM_VV(IIOO  ,JJOO-1,KKOO  )
               CASE(-3) ; UN_OTHER = UN_OTHER + OM_WW(IIOO  ,JJOO  ,KKOO  )
               CASE( 3) ; UN_OTHER = UN_OTHER + OM_WW(IIOO  ,JJOO  ,KKOO-1)
            END SELECT
         ENDDO
      ENDDO
   ENDDO
   N_INT_CELLS = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
   UN_OTHER = UN_OTHER/REAL(N_INT_CELLS,EB)
   SELECT CASE(BC%IOR)
      CASE(-1) ; UU(BC%II  ,BC%JJ,BC%KK) = UN_OTHER
      CASE( 1) ; UU(BC%II-1,BC%JJ,BC%KK) = UN_OTHER
      CASE(-2) ; VV(BC%II,BC%JJ  ,BC%KK) = UN_OTHER
      CASE( 2) ; VV(BC%II,BC%JJ-1,BC%KK) = UN_OTHER
      CASE(-3) ; WW(BC%II,BC%JJ,BC%KK  ) = UN_OTHER
      CASE( 3) ; WW(BC%II,BC%JJ,BC%KK-1) = UN_OTHER
   END SELECT
ENDDO WALL_LOOP

DRAG_UVWMAX = 0._EB

! Loop over all cell edges and determine the appropriate velocity BCs

EDGE_LOOP: DO IE=1,EDGE_COUNT(NM)

   ED => EDGE(IE)

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

   IF ((CELL(ICMM)%EXTERIOR .OR. CELL(ICMM)%SOLID) .AND. &
       (CELL(ICPM)%EXTERIOR .OR. CELL(ICPM)%SOLID) .AND. &
       (CELL(ICMP)%EXTERIOR .OR. CELL(ICMP)%SOLID) .AND. &
       (CELL(ICPP)%EXTERIOR .OR. CELL(ICPP)%SOLID)) CYCLE EDGE_LOOP

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
         DXX(1)  = DY(JJ)
         DXX(2)  = DZ(KK)
      CASE(2) COMPONENT
         UUP(1)  = WW(II+1,JJ,KK)
         UUM(1)  = WW(II,JJ,KK)
         UUP(2)  = UU(II,JJ,KK+1)
         UUM(2)  = UU(II,JJ,KK)
         DXX(1)  = DZ(KK)
         DXX(2)  = DX(II)
      CASE(3) COMPONENT
         UUP(1)  = UU(II,JJ+1,KK)
         UUM(1)  = UU(II,JJ,KK)
         UUP(2)  = VV(II+1,JJ,KK)
         UUM(2)  = VV(II,JJ,KK)
         DXX(1)  = DX(II)
         DXX(2)  = DY(JJ)
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
            IWM  = CELL(ICMM)%WALL_INDEX(-IOR)
            IWMI = CELL(ICMM)%WALL_INDEX( IS2)
            IF (ICD==1) THEN
               IWP  = CELL(ICMP)%WALL_INDEX(-IOR)
               IWPI = CELL(ICMP)%WALL_INDEX(-IS2)
            ELSE ! ICD==2
               IWP  = CELL(ICPM)%WALL_INDEX(-IOR)
               IWPI = CELL(ICPM)%WALL_INDEX(-IS2)
            ENDIF
         ELSE
            IF (ICD==1) THEN
               IWM  = CELL(ICPM)%WALL_INDEX(-IOR)
               IWMI = CELL(ICPM)%WALL_INDEX( IS2)
            ELSE ! ICD==2
               IWM  = CELL(ICMP)%WALL_INDEX(-IOR)
               IWMI = CELL(ICMP)%WALL_INDEX( IS2)
            ENDIF
            IWP  = CELL(ICPP)%WALL_INDEX(-IOR)
            IWPI = CELL(ICPP)%WALL_INDEX(-IS2)
         ENDIF

         ! If both adjacent wall cells are undefined, cycle out of the loop.

         IF (IWM==0 .AND. IWP==0) CYCLE ORIENTATION_LOOP

         ! If there is a solid wall separating the two adjacent wall cells, cycle out of the loop.

         IF ((WALL(IWMI)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. SURFACE(WALL(IWM)%SURF_INDEX)%VELOCITY_BC_INDEX/=FREE_SLIP_BC) .OR. &
             (WALL(IWPI)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. SURFACE(WALL(IWP)%SURF_INDEX)%VELOCITY_BC_INDEX/=FREE_SLIP_BC)) &
            CYCLE ORIENTATION_LOOP

         ! If only one adjacent wall cell is defined, use its properties.

         IF (IWM>0) THEN
            WCM => WALL(IWM)
         ELSE
            WCM => WALL(IWP)
         ENDIF

         IF (IWP>0) THEN
            WCP => WALL(IWP)
         ELSE
            WCP => WALL(IWM)
         ENDIF

         WCM_B1 => BOUNDARY_PROP1(WCM%B1_INDEX)
         WCP_B1 => BOUNDARY_PROP1(WCP%B1_INDEX)

         ! If both adjacent wall cells are NULL, cycle out.

         BOUNDARY_TYPE_M = WCM%BOUNDARY_TYPE
         BOUNDARY_TYPE_P = WCP%BOUNDARY_TYPE

         IF (BOUNDARY_TYPE_M==NULL_BOUNDARY .AND. BOUNDARY_TYPE_P==NULL_BOUNDARY) CYCLE ORIENTATION_LOOP

         ! Set up synthetic eddy method

         SYNTHETIC_EDDY_METHOD = .FALSE.
         IF (IWM>0 .AND. IWP>0) THEN
            IF (WCM%VENT_INDEX==WCP%VENT_INDEX) THEN
               IF (WCM%VENT_INDEX>0) THEN
                  VT=>VENTS(WCM%VENT_INDEX)
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

         OPEN_AND_WIND_BC: IF ((IWM==0 .OR. WALL(IWM)%BOUNDARY_TYPE==OPEN_BOUNDARY) .AND. &
                               (IWP==0 .OR. WALL(IWP)%BOUNDARY_TYPE==OPEN_BOUNDARY)       ) THEN

            VENT_INDEX = MAX(WCM%VENT_INDEX,WCP%VENT_INDEX)
            VT => VENTS(VENT_INDEX)

            UPWIND_BOUNDARY = .FALSE.
            INFLOW_BOUNDARY = .FALSE.

            IF (OPEN_WIND_BOUNDARY) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,   0,KK) + VV(II,   0,KK+1))
                     IF (JJ==JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,JBAR,KK) + VV(II,JBAR,KK+1))
                     IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,0)    + WW(II,JJ+1,   0))
                     IF (KK==KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,KBAR) + WW(II,JJ+1,KBAR))
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(   0,JJ,KK) + UU(   0,JJ,KK+1))
                     IF (II==IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(IBAR,JJ,KK) + UU(IBAR,JJ,KK+1))
                     IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,   0) + WW(II+1,JJ,   0))
                     IF (KK==KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,KBAR) + WW(II+1,JJ,KBAR))
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(   0,JJ,KK) + UU(   0,JJ+1,KK))
                     IF (II==IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(IBAR,JJ,KK) + UU(IBAR,JJ+1,KK))
                     IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,   0,KK) + VV(II+1,   0,KK))
                     IF (JJ==JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,JBAR,KK) + VV(II+1,JBAR,KK))
               END SELECT
               IF ((IOR==1.AND.U_WIND(KK)>=0._EB) .OR. (IOR==-1.AND.U_WIND(KK)<=0._EB)) UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==2.AND.V_WIND(KK)>=0._EB) .OR. (IOR==-2.AND.V_WIND(KK)<=0._EB)) UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==3.AND.W_WIND(KK)>=0._EB) .OR. (IOR==-3.AND.W_WIND(KK)<=0._EB)) UPWIND_BOUNDARY = .TRUE.
               IF ((IOR==1.AND.U_NORM>=0._EB) .OR. (IOR==-1.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
               IF ((IOR==2.AND.U_NORM>=0._EB) .OR. (IOR==-2.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
               IF ((IOR==3.AND.U_NORM>=0._EB) .OR. (IOR==-3.AND.U_NORM<=0._EB)) INFLOW_BOUNDARY = .TRUE.
            ENDIF

            WIND_NO_WIND_IF: IF (.NOT.UPWIND_BOUNDARY .OR. .NOT.INFLOW_BOUNDARY) THEN  ! For regular OPEN boundary, (free-slip) BCs

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = WW(II,1,KK)
                     IF (JJ==JBAR .AND. IOR==-2) WW(II,JBP1,KK) = WW(II,JBAR,KK)
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = VV(II,JJ,1)
                     IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KBP1) = VV(II,JJ,KBAR)
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = WW(1,JJ,KK)
                     IF (II==IBAR .AND. IOR==-1) WW(IBP1,JJ,KK) = WW(IBAR,JJ,KK)
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = UU(II,JJ,1)
                     IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KBP1) = UU(II,JJ,KBAR)
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = VV(1,JJ,KK)
                     IF (II==IBAR .AND. IOR==-1) VV(IBP1,JJ,KK) = VV(IBAR,JJ,KK)
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = UU(II,1,KK)
                     IF (JJ==JBAR .AND. IOR==-2) UU(II,JBP1,KK) = UU(II,JBAR,KK)
               END SELECT

            ELSE WIND_NO_WIND_IF  ! For upwind, inflow boundaries, use the specified wind field for tangential velocity components

               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)    = W_WIND(KK) + VEL_EDDY
                     IF (JJ==JBAR .AND. IOR==-2) WW(II,JBP1,KK) = W_WIND(KK) + VEL_EDDY
                     IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)    = V_WIND(KK) + VEL_EDDY
                     IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KBP1) = V_WIND(KK) + VEL_EDDY
                  CASE(2)
                     IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)    = W_WIND(KK) + VEL_EDDY
                     IF (II==IBAR .AND. IOR==-1) WW(IBP1,JJ,KK) = W_WIND(KK) + VEL_EDDY
                     IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)    = U_WIND(KK) + VEL_EDDY
                     IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KBP1) = U_WIND(KK) + VEL_EDDY
                  CASE(3)
                     IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)    = V_WIND(KK) + VEL_EDDY
                     IF (II==IBAR .AND. IOR==-1) VV(IBP1,JJ,KK) = V_WIND(KK) + VEL_EDDY
                     IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)    = U_WIND(KK) + VEL_EDDY
                     IF (JJ==JBAR .AND. IOR==-2) UU(II,JBP1,KK) = U_WIND(KK) + VEL_EDDY
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
            IIGM = CELL(ICMM)%I
            JJGM = CELL(ICMM)%J
            KKGM = CELL(ICMM)%K
            IF (ICD==1) THEN
               IIGP = CELL(ICMP)%I
               JJGP = CELL(ICMP)%J
               KKGP = CELL(ICMP)%K
            ELSE ! ICD==2
               IIGP = CELL(ICPM)%I
               JJGP = CELL(ICPM)%J
               KKGP = CELL(ICPM)%K
            ENDIF
         ELSE
            VEL_GAS   = UUP(IVL)
            VEL_GHOST = UUM(IVL)
            IF (ICD==1) THEN
               IIGM = CELL(ICPM)%I
               JJGM = CELL(ICPM)%J
               KKGM = CELL(ICPM)%K
            ELSE ! ICD==2
               IIGM = CELL(ICMP)%I
               JJGM = CELL(ICMP)%J
               KKGM = CELL(ICMP)%K
            ENDIF
            IIGP = CELL(ICPP)%I
            JJGP = CELL(ICPP)%J
            KKGP = CELL(ICPP)%K
         ENDIF

         ! Decide whether or not to process edge using data interpolated from another mesh

         INTERPOLATION_IF: IF (NOM(ICD)==0 .OR. &
                   (BOUNDARY_TYPE_M==SOLID_BOUNDARY .OR. BOUNDARY_TYPE_P==SOLID_BOUNDARY) .OR. &
                   (BOUNDARY_TYPE_M/=INTERPOLATED_BOUNDARY .AND. BOUNDARY_TYPE_P/=INTERPOLATED_BOUNDARY) .OR. &
                   (SYNTHETIC_EDDY_METHOD .AND. (BOUNDARY_TYPE_M==OPEN_BOUNDARY .OR. BOUNDARY_TYPE_P==OPEN_BOUNDARY)) ) THEN

            ! Determine appropriate velocity BC by assessing each adjacent wall cell. If the BCs are different on each
            ! side of the edge, choose the one with the specified velocity or velocity gradient, if there is one.
            ! If not, choose the max value of boundary condition index, simply for consistency.

            SURF_INDEXM = WCM%SURF_INDEX
            SURF_INDEXP = WCP%SURF_INDEX
            IF (SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_VELOCITY .OR. SURFACE(SURF_INDEXM)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXM)
            ELSEIF (SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_VELOCITY .OR. SURFACE(SURF_INDEXP)%SPECIFIED_NORMAL_GRADIENT) THEN
               SF=>SURFACE(SURF_INDEXP)
            ELSE
               SF=>SURFACE(MAX(SURF_INDEXM,SURF_INDEXP))
            ENDIF
            VELOCITY_BC_INDEX = SF%VELOCITY_BC_INDEX
            IF (WCM%VENT_INDEX==WCP%VENT_INDEX .AND. WCP%VENT_INDEX > 0) THEN
               IF(VENTS(WCM%VENT_INDEX)%NODE_INDEX>0 .AND. WCM_B1%U_NORMAL >= 0._EB) VELOCITY_BC_INDEX=FREE_SLIP_BC
            ENDIF
            IF (SYNTHETIC_EDDY_METHOD)         VELOCITY_BC_INDEX=NO_SLIP_BC
            IF (SF%ROUGHNESS>2._EB/WCM_B1%RDN) VELOCITY_BC_INDEX=NO_SLIP_BC ! see Basu et al. BLM 2017

            ! Compute the viscosity by averaging the two adjacent gas cells

            MUA = 0.5_EB*(MU(IIGM,JJGM,KKGM) + MU(IIGP,JJGP,KKGP))

            ! Check for HVAC tangential velocity

            HVAC_TANGENTIAL = .FALSE.
            IF (WCM%VENT_INDEX>0 .OR. WCP%VENT_INDEX>0) THEN
               IF (WCM%VENT_INDEX>0) THEN
                  WCX => WCM
               ELSE
                  WCX => WCP
               ENDIF
               VT => VENTS(WCX%VENT_INDEX)
               WCX_B1 => BOUNDARY_PROP1(WCX%B1_INDEX)
               IF (VT%NODE_INDEX>0 .AND. WCX_B1%U_NORMAL_S<0._EB) THEN
                  VELOCITY_BC_INDEX = NO_SLIP_BC
                  IF (ALL(VT%UVW>-1.E12_EB)) HVAC_TANGENTIAL = .TRUE.  ! User-specified tangential components of velocity
               ENDIF
            ENDIF

            ! Determine if there is a tangential velocity component

            IF (.NOT.SF%SPECIFIED_TANGENTIAL_VELOCITY .AND. .NOT.SYNTHETIC_EDDY_METHOD .AND. .NOT.HVAC_TANGENTIAL) THEN

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

               IF (ABS(SF%VEL)>0._EB .OR. VEL_N==0._EB) THEN  ! User-specified normal velocity or no normal velocity at all
                  IF (ABS(SF%T_IGN-T_BEGIN)<=SPACING(SF%T_IGN) .AND. SF%RAMP(TIME_VELO)%INDEX>=1) THEN
                     TSI = T
                  ELSE
                     TSI=T-SF%T_IGN
                  ENDIF
                  PROFILE_FACTOR = 1._EB
                  RAMP_T = EVALUATE_RAMP(TSI,SF%RAMP(TIME_VELO)%INDEX,TAU=SF%RAMP(TIME_VELO)%TAU)
                  IF (SF%VEL < 0._EB) THEN
                     IF (SF%RAMP(VELO_PROF_Z)%INDEX>0) PROFILE_FACTOR = EVALUATE_RAMP(ZC(KK),SF%RAMP(VELO_PROF_Z)%INDEX)
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(2) + VEL_EDDY))
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = RAMP_T*(PROFILE_FACTOR*(SF%VEL_T(1) + VEL_EDDY))
                  ELSEIF (SF%VEL > 0._EB) THEN
                     IF (SF%PROFILE/=0) PROFILE_FACTOR = ABS(0.5_EB*(WCM_B1%U_NORMAL_0+WCP_B1%U_NORMAL_0)/SF%VEL)
                     IF (SF%RAMP(VELO_PROF_Z)%INDEX>0) PROFILE_FACTOR = EVALUATE_RAMP(ZC(KK),SF%RAMP(VELO_PROF_Z)%INDEX)
                     IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = RAMP_T*PROFILE_FACTOR*VEL_EDDY
                     IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = RAMP_T*PROFILE_FACTOR*VEL_EDDY
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
                     ITMP = MIN(I_MAX_TEMP,NINT(0.5_EB*(TMP(IIGM,JJGM,KKGM)+TMP(IIGP,JJGP,KKGP))))
                     MU_WALL = MU_RSQMW_Z(ITMP,1)/RSQ_MW_Z(1)
                     RHO_WALL = 0.5_EB*( RHOP(IIGM,JJGM,KKGM) + RHOP(IIGP,JJGP,KKGP) )
                     CALL WALL_MODEL(SLIP_COEF,U_TAU,Y_PLUS,MU_WALL/RHO_WALL,SF%ROUGHNESS,0.5_EB*DXX(ICD),VEL_GAS-VEL_T)
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
                     DRAG_FACTOR = 0.5_EB*SF%DRAG_COEFFICIENT*SF%SHAPE_FACTOR*SF%VEG_LSET_BETA*(SF%VEG_LSET_SIGMA*100._EB)
                  ELSE
                     VEG_HT = SF%LAYER_THICKNESS(1)
                     DRAG_FACTOR = 0.5_EB*SF%DRAG_COEFFICIENT*SF%SHAPE_FACTOR*SF%PACKING_RATIO(1)*SF%SURFACE_VOLUME_RATIO(1)
                  ENDIF
                  HT_SCALE_FACTOR = MIN(1._EB,0.5_EB*(WCM_B1%RDN+WCP_B1%RDN)*VEG_HT)
                  MU_DUIDXJ(ICD_SGN) = I_SGN*RHO_WALL*DRAG_FACTOR*VEG_HT*HT_SCALE_FACTOR**2*VEL_GAS*VEL_T
                  DRAG_UVWMAX = MAX(DRAG_UVWMAX,DRAG_FACTOR*HT_SCALE_FACTOR**2*VEL_T)
                  ALTERED_GRADIENT(ICD_SGN) = .TRUE.

            END SELECT BOUNDARY_CONDITION

         ELSE INTERPOLATION_IF  ! Use data from another mesh

            INTERPOLATED_EDGE = .TRUE.
            OM => OMESH(ABS(NOM(ICD)))

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
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ENDIF
               CASE(2)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ELSE ! ICD=2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)-1)
                  ENDIF
               CASE(3)
                  IF (ICD==1) THEN
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD),JJO(ICD)-1,KKO(ICD))
                  ELSE ! ICD==2
                     VEL_GHOST = WGT*VEL_OTHER(IIO(ICD),JJO(ICD),KKO(ICD)) + OMW*VEL_OTHER(IIO(ICD)-1,JJO(ICD),KKO(ICD))
                  ENDIF
            END SELECT

            ! At the exterior edge of Mesh NM, which abuts Mesh NOM, assign the appropriate velocity component to the
            ! ghost cell. For example, MESHES(NM)%VV(IBP1,JJ,KBP1) is stored in MESHES(NM)%OMESH(NOM)%V(IIO+1,JJO,KKO) if Mesh
            ! NOM is above Mesh NM (IOR=-3), or it is stored in MESHES(NM)%OMESH(NOM)%V(IIO,JJO,KKO+1) if Mesh NOM is to
            ! the right of NM (IOR=-1). The value of each is the same and it is transfered to Mesh NM with a one-step time lag.

            IF (CORRECTOR) THEN
               SELECT CASE(IEC)
                  CASE(1)
                     IF (JJ==0    .AND. KK==0    .AND. ABS(IOR)==2) UU(II,JJ  ,KK  ) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)-1)
                     IF (JJ==0    .AND. KK==0    .AND. ABS(IOR)==3) UU(II,JJ  ,KK  ) = OM%U(IIO(ICD),JJO(ICD)-1,KKO(ICD)  )
                     IF (JJ==0    .AND. KK==KBAR .AND. ABS(IOR)==2) UU(II,JJ  ,KK+1) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)+1)
                     IF (JJ==0    .AND. KK==KBAR .AND. ABS(IOR)==3) UU(II,JJ  ,KK+1) = OM%U(IIO(ICD),JJO(ICD)-1,KKO(ICD)  )
                     IF (JJ==JBAR .AND. KK==0    .AND. ABS(IOR)==2) UU(II,JJ+1,KK  ) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)-1)
                     IF (JJ==JBAR .AND. KK==0    .AND. ABS(IOR)==3) UU(II,JJ+1,KK  ) = OM%U(IIO(ICD),JJO(ICD)+1,KKO(ICD)  )
                     IF (JJ==JBAR .AND. KK==KBAR .AND. ABS(IOR)==2) UU(II,JJ+1,KK+1) = OM%U(IIO(ICD),JJO(ICD)  ,KKO(ICD)+1)
                     IF (JJ==JBAR .AND. KK==KBAR .AND. ABS(IOR)==3) UU(II,JJ+1,KK+1) = OM%U(IIO(ICD),JJO(ICD)+1,KKO(ICD)  )
                  CASE(2)
                     IF (II==0    .AND. KK==0    .AND. ABS(IOR)==1) VV(II  ,JJ,KK  ) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)-1)
                     IF (II==0    .AND. KK==0    .AND. ABS(IOR)==3) VV(II  ,JJ,KK  ) = OM%V(IIO(ICD)-1,JJO(ICD),KKO(ICD)  )
                     IF (II==0    .AND. KK==KBAR .AND. ABS(IOR)==1) VV(II  ,JJ,KK+1) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)+1)
                     IF (II==0    .AND. KK==KBAR .AND. ABS(IOR)==3) VV(II  ,JJ,KK+1) = OM%V(IIO(ICD)-1,JJO(ICD),KKO(ICD)  )
                     IF (II==IBAR .AND. KK==0    .AND. ABS(IOR)==1) VV(II+1,JJ,KK  ) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)-1)
                     IF (II==IBAR .AND. KK==0    .AND. ABS(IOR)==3) VV(II+1,JJ,KK  ) = OM%V(IIO(ICD)+1,JJO(ICD),KKO(ICD)  )
                     IF (II==IBAR .AND. KK==KBAR .AND. ABS(IOR)==1) VV(II+1,JJ,KK+1) = OM%V(IIO(ICD)  ,JJO(ICD),KKO(ICD)+1)
                     IF (II==IBAR .AND. KK==KBAR .AND. ABS(IOR)==3) VV(II+1,JJ,KK+1) = OM%V(IIO(ICD)+1,JJO(ICD),KKO(ICD)  )
                  CASE(3)
                     IF (II==0    .AND. JJ==0    .AND. ABS(IOR)==1) WW(II  ,JJ  ,KK) = OM%W(IIO(ICD)  ,JJO(ICD)-1,KKO(ICD))
                     IF (II==0    .AND. JJ==0    .AND. ABS(IOR)==2) WW(II  ,JJ  ,KK) = OM%W(IIO(ICD)-1,JJO(ICD)  ,KKO(ICD))
                     IF (II==0    .AND. JJ==JBAR .AND. ABS(IOR)==1) WW(II  ,JJ+1,KK) = OM%W(IIO(ICD)  ,JJO(ICD)+1,KKO(ICD))
                     IF (II==0    .AND. JJ==JBAR .AND. ABS(IOR)==2) WW(II  ,JJ+1,KK) = OM%W(IIO(ICD)-1,JJO(ICD)  ,KKO(ICD))
                     IF (II==IBAR .AND. JJ==0    .AND. ABS(IOR)==1) WW(II+1,JJ  ,KK) = OM%W(IIO(ICD)  ,JJO(ICD)-1,KKO(ICD))
                     IF (II==IBAR .AND. JJ==0    .AND. ABS(IOR)==2) WW(II+1,JJ  ,KK) = OM%W(IIO(ICD)+1,JJO(ICD)  ,KKO(ICD))
                     IF (II==IBAR .AND. JJ==JBAR .AND. ABS(IOR)==1) WW(II+1,JJ+1,KK) = OM%W(IIO(ICD)  ,JJO(ICD)+1,KKO(ICD))
                     IF (II==IBAR .AND. JJ==JBAR .AND. ABS(IOR)==2) WW(II+1,JJ+1,KK) = OM%W(IIO(ICD)+1,JJO(ICD)  ,KKO(ICD))
               END SELECT
            ENDIF

         ENDIF INTERPOLATION_IF

         ! Set ghost cell values at edge of computational domain

         SELECT CASE(IEC)
            CASE(1)
               IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
               IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
               IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. .NOT.INTERPOLATED_EDGE) THEN
                 IF (ICD==1) THEN
                    ED%W_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    ED%V_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(2)
               IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
               IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
               IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
               IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
               IF (CORRECTOR .AND. .NOT.INTERPOLATED_EDGE) THEN
                 IF (ICD==1) THEN
                    ED%U_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ELSE ! ICD=2
                    ED%W_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
                 ENDIF
               ENDIF
            CASE(3)
               IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
               IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
               IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
               IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
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
      IF (.NOT.CELL(ICMM)%SOLID .AND. .NOT.CELL(ICPM)%SOLID .AND. &
          .NOT.CELL(ICMP)%SOLID .AND. .NOT.CELL(ICPP)%SOLID) CYCLE EDGE_LOOP
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

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

IF(CC_IBM) CALL CC_VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES,DO_IBEDGES=.TRUE.)

END SUBROUTINE VELOCITY_BC


!> \brief Force normal component of velocity to match at interpolated boundaries
!> \param NM Mesh number

SUBROUTINE MATCH_VELOCITY(NM)

USE COMPLEX_GEOMETRY, ONLY : CC_IDCF
USE CC_VELOCITY, ONLY : CC_MATCH_VELOCITY
INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
INTEGER, INTENT(IN) :: NM
REAL(EB) :: T_NOW,DA_OTHER,UU_OTHER,VV_OTHER,WW_OTHER,NOM_CELLS
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,OM_UU,OM_VV,OM_WW
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

INTEGER  :: ICF
REAL(EB) :: AU,AU1,AV,AV1,AW,AW1

IF (SOLID_PHASE_ONLY) RETURN

IF(CC_IBM) THEN
   CALL CC_MATCH_VELOCITY(NM,PREDICTOR,.TRUE.)
   RETURN
ENDIF

T_NOW = CURRENT_TIME()

! Assign local variable names

CALL POINT_TO_MESH(NM)

! Point to the appropriate velocity field

IF (PREDICTOR) THEN
   UU => US
   VV => VS
   WW => WS
ELSE
   UU => U
   VV => V
   WW => W
ENDIF

! Loop over all external wall cells and force adjacent normal components of velocty at interpolated boundaries to match.
! BOUNDARY_TYPE_PREVIOUS will be used at the next phase of the time step to indicate if the velocity component has been changed.

EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   EWC%BOUNDARY_TYPE_PREVIOUS = WC%BOUNDARY_TYPE

   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   BC =>BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR
   NOM = EWC%NOM
   OM => OMESH(NOM)
   M2 => MESHES(NOM)

   ! Determine the area of the interpolated cell face

   DA_OTHER = 0._EB

   SELECT CASE(ABS(IOR))
      CASE(1)
         IF (PREDICTOR) OM_UU => OM%US
         IF (CORRECTOR) OM_UU => OM%U
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DY(JJO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         IF (PREDICTOR) OM_VV => OM%VS
         IF (CORRECTOR) OM_VV => OM%V
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DZ(KKO)
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         IF (PREDICTOR) OM_WW => OM%WS
         IF (CORRECTOR) OM_WW => OM%W
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  DA_OTHER = DA_OTHER + M2%DX(IIO)*M2%DY(JJO)
               ENDDO
            ENDDO
         ENDDO
   END SELECT

   ! Determine the normal component of velocity from the other mesh and use it for average

   SELECT CASE(IOR)

      CASE( 1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + OM_UU(IIO,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_UU(IIO,JJO,KKO) = 0.5_EB*(OM_UU(IIO,JJO,KKO)+UU(0,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW) = UU(0,JJ,KK)
         UU(0,JJ,KK)  = 0.5_EB*(UU(0,JJ,KK) + UU_OTHER)

      CASE(-1)

         UU_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  UU_OTHER = UU_OTHER + OM_UU(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_UU(IIO-1,JJO,KKO) = 0.5_EB*(OM_UU(IIO-1,JJO,KKO)+UU(IBAR,JJ,KK))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW) = UU(IBAR,JJ,KK)
         UU(IBAR,JJ,KK) = 0.5_EB*(UU(IBAR,JJ,KK) + UU_OTHER)

      CASE( 2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_VV(IIO,JJO,KKO) = 0.5_EB*(OM_VV(IIO,JJO,KKO)+VV(II,0,KK))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW) = VV(II,0,KK)
         VV(II,0,KK)  = 0.5_EB*(VV(II,0,KK) + VV_OTHER)

      CASE(-2)

         VV_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  VV_OTHER = VV_OTHER + OM_VV(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_VV(IIO,JJO-1,KKO) = 0.5_EB*(OM_VV(IIO,JJO-1,KKO)+VV(II,JBAR,KK))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW)   = VV(II,JBAR,KK)
         VV(II,JBAR,KK) = 0.5_EB*(VV(II,JBAR,KK) + VV_OTHER)

      CASE( 3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_WW(IIO,JJO,KKO) = 0.5_EB*(OM_WW(IIO,JJO,KKO)+WW(II,JJ,0))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW) = WW(II,JJ,0)
         WW(II,JJ,0)  = 0.5_EB*(WW(II,JJ,0) + WW_OTHER)

      CASE(-3)

         WW_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  WW_OTHER = WW_OTHER + OM_WW(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
                  IF (EWC%AREA_RATIO>0.9_EB) OM_WW(IIO,JJO,KKO-1) = 0.5_EB*(OM_WW(IIO,JJO,KKO-1)+WW(II,JJ,KBAR))
               ENDDO
            ENDDO
         ENDDO
         UVW_SAVE(IW)   = WW(II,JJ,KBAR)
         WW(II,JJ,KBAR) = 0.5_EB*(WW(II,JJ,KBAR) + WW_OTHER)

   END SELECT

   ! Save velocity components at the ghost cell midpoint

   U_GHOST(IW) = 0._EB
   V_GHOST(IW) = 0._EB
   W_GHOST(IW) = 0._EB

   IF (PREDICTOR) OM_UU => OM%US
   IF (CORRECTOR) OM_UU => OM%U
   IF (PREDICTOR) OM_VV => OM%VS
   IF (CORRECTOR) OM_VV => OM%V
   IF (PREDICTOR) OM_WW => OM%WS
   IF (CORRECTOR) OM_WW => OM%W

   IF (CC_IBM) THEN
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               AU =1._EB; ICF=M2%FCVAR(IIO  ,JJO,KKO,CC_IDCF,IAXIS); IF(ICF>0) AU =M2%CUT_FACE(ICF)%ALPHA_CF
               AU1=1._EB; ICF=M2%FCVAR(IIO-1,JJO,KKO,CC_IDCF,IAXIS); IF(ICF>0) AU1=M2%CUT_FACE(ICF)%ALPHA_CF
               AV =1._EB; ICF=M2%FCVAR(IIO,JJO  ,KKO,CC_IDCF,JAXIS); IF(ICF>0) AV =M2%CUT_FACE(ICF)%ALPHA_CF
               AV1=1._EB; ICF=M2%FCVAR(IIO,JJO-1,KKO,CC_IDCF,JAXIS); IF(ICF>0) AV1=M2%CUT_FACE(ICF)%ALPHA_CF
               AW =1._EB; ICF=M2%FCVAR(IIO,JJO,KKO  ,CC_IDCF,KAXIS); IF(ICF>0) AW =M2%CUT_FACE(ICF)%ALPHA_CF
               AW1=1._EB; ICF=M2%FCVAR(IIO,JJO,KKO-1,CC_IDCF,KAXIS); IF(ICF>0) AW1=M2%CUT_FACE(ICF)%ALPHA_CF
               U_GHOST(IW) = U_GHOST(IW) + (OM_UU(IIO,JJO,KKO)+OM_UU(IIO-1,JJO,KKO))/(AU+AU1)
               V_GHOST(IW) = V_GHOST(IW) + (OM_VV(IIO,JJO,KKO)+OM_VV(IIO,JJO-1,KKO))/(AV+AV1)
               W_GHOST(IW) = W_GHOST(IW) + (OM_WW(IIO,JJO,KKO)+OM_WW(IIO,JJO,KKO-1))/(AW+AW1)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
               U_GHOST(IW) = U_GHOST(IW) + 0.5_EB*(OM_UU(IIO,JJO,KKO)+OM_UU(IIO-1,JJO,KKO))
               V_GHOST(IW) = V_GHOST(IW) + 0.5_EB*(OM_VV(IIO,JJO,KKO)+OM_VV(IIO,JJO-1,KKO))
               W_GHOST(IW) = W_GHOST(IW) + 0.5_EB*(OM_WW(IIO,JJO,KKO)+OM_WW(IIO,JJO,KKO-1))
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   NOM_CELLS = REAL((EWC%IIO_MAX-EWC%IIO_MIN+1)*(EWC%JJO_MAX-EWC%JJO_MIN+1)*(EWC%KKO_MAX-EWC%KKO_MIN+1),EB)
   U_GHOST(IW) = U_GHOST(IW)/NOM_CELLS
   V_GHOST(IW) = V_GHOST(IW)/NOM_CELLS
   W_GHOST(IW) = W_GHOST(IW)/NOM_CELLS

ENDDO EXTERNAL_WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

END SUBROUTINE MATCH_VELOCITY


!> \brief Force normal component of velocity flux to match at interpolated boundaries
!> \param NM Mesh number

SUBROUTINE MATCH_VELOCITY_FLUX(NM)

USE CC_VELOCITY, ONLY : CC_MATCH_VELOCITY_FLUX
INTEGER  :: NOM,II,JJ,KK,IOR,IW,IIO,JJO,KKO
INTEGER, INTENT(IN) :: NM
REAL(EB) :: T_NOW,DA_OTHER,FVX_OTHER,FVY_OTHER,FVZ_OTHER
TYPE (OMESH_TYPE), POINTER :: OM
TYPE (MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

IF (NMESHES==1) RETURN
IF (SOLID_PHASE_ONLY) RETURN

IF(CC_IBM) THEN
   CALL CC_MATCH_VELOCITY_FLUX(NM)
   RETURN
ENDIF

T_NOW = CURRENT_TIME()

! Assign local variable names

CALL POINT_TO_MESH(NM)

! Loop over all cell edges and determine the appropriate velocity BCs

EXTERNAL_WALL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS

   WC=>WALL(IW)
   EWC=>EXTERNAL_WALL(IW)
   IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP

   BC => BOUNDARY_COORD(WC%BC_INDEX)
   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR
   NOM = EWC%NOM
   OM => OMESH(NOM)
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

   ! Determine the normal component of velocity from the other mesh and use it for average

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
         FVX(0,JJ,KK) = 0.5_EB*(FVX(0,JJ,KK) + FVX_OTHER)

      CASE(-1)

         FVX_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVX_OTHER = FVX_OTHER + OM%FVX(IIO-1,JJO,KKO)*M2%DY(JJO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVX(IBAR,JJ,KK) = 0.5_EB*(FVX(IBAR,JJ,KK) + FVX_OTHER)

      CASE( 2)

         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVY(II,0,KK) = 0.5_EB*(FVY(II,0,KK) + FVY_OTHER)

      CASE(-2)

         FVY_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVY_OTHER = FVY_OTHER + OM%FVY(IIO,JJO-1,KKO)*M2%DX(IIO)*M2%DZ(KKO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVY(II,JBAR,KK) = 0.5_EB*(FVY(II,JBAR,KK) + FVY_OTHER)

      CASE( 3)

         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVZ(II,JJ,0) = 0.5_EB*(FVZ(II,JJ,0) + FVZ_OTHER)

      CASE(-3)

         FVZ_OTHER = 0._EB
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                  FVZ_OTHER = FVZ_OTHER + OM%FVZ(IIO,JJO,KKO-1)*M2%DX(IIO)*M2%DY(JJO)/DA_OTHER
               ENDDO
            ENDDO
         ENDDO
         FVZ(II,JJ,KBAR) = 0.5_EB*(FVZ(II,JJ,KBAR) + FVZ_OTHER)

   END SELECT

ENDDO EXTERNAL_WALL_LOOP

T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

END SUBROUTINE MATCH_VELOCITY_FLUX


!> \brief Check the Courant and Von Neumann stability criteria, and if necessary, reduce or increase the time step
!> \param DT Time step (s)
!> \param DT_NEW New time step (s)
!> \param T Current time (s)
!> \param NM Mesh number

SUBROUTINE CHECK_STABILITY(DT,DT_NEW,T,NM)

USE VELO_KERNELS, ONLY: CHECK_STABILITY_KERNEL

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT,T
REAL(EB) :: DT_NEW(NMESHES)
REAL(EB) :: T_NOW

T_NOW = CURRENT_TIME()
CALL CHECK_STABILITY_KERNEL(MESHES(NM),DT,DT_NEW(NM),T,NM)
T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

END SUBROUTINE CHECK_STABILITY


!> \brief Add baroclinic term to the momentum equation
!> \param T Current time (s)
!> \param NM Mesh number

SUBROUTINE BAROCLINIC_CORRECTION(T,NM)

USE CC_VELOCITY, ONLY: CC_BAROCLINIC_CORRECTION
USE VELO_KERNELS, ONLY: BAROCLINIC_CORRECTION_KERNEL
REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
REAL(EB) :: T_NOW

IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

T_NOW = CURRENT_TIME()

CALL BAROCLINIC_CORRECTION_KERNEL(MESHES(NM),T)

T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

IF(CC_IBM) CALL CC_BAROCLINIC_CORRECTION(T,NM)

END SUBROUTINE BAROCLINIC_CORRECTION


!> \brief Recompute velocities on wall cells
!> \param DT Time step (s)
!> \param STORE_UN Flag indicating whether normal velocity component is to be saved
!> \details Ensure that the correct normal derivative of H is used on the projection. It is only used when the Poisson equation
!> for the pressure is solved .NOT. PRES_ON_WHOLE_DOMAIN (i.e. using the GLMAT solver).

SUBROUTINE WALL_VELOCITY_NO_GRADH(DT,STORE_UN)

REAL(EB), INTENT(IN) :: DT
LOGICAL, INTENT(IN) :: STORE_UN
INTEGER :: II,JJ,KK,IIG,JJG,KKG,IOR,IW,N_INTERNAL_WALL_CELLS_AUX,IC,ICG
REAL(EB) :: VEL_N
TYPE (WALL_TYPE), POINTER :: WC
REAL(EB), SAVE, ALLOCATABLE, DIMENSION(:) :: UN_WALLS
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

N_INTERNAL_WALL_CELLS_AUX=0
IF (.NOT.PRES_ON_WHOLE_DOMAIN) N_INTERNAL_WALL_CELLS_AUX=N_INTERNAL_WALL_CELLS

STORE_UN_COND : IF ( STORE_UN .AND. CORRECTOR) THEN

   ! These velocities from the beginning of step are needed for the velocity fix on wall cells at the corrector
   ! phase (i.e. the loops in VELOCITY_CORRECTOR will change U,V,W to wrong reults using (HP1-HP)/DX gradients,
   ! when the pressure solver in the GLMAT solver.

   IF (ALLOCATED(UN_WALLS)) DEALLOCATE(UN_WALLS)
   ALLOCATE( UN_WALLS(1:N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX) )
   UN_WALLS(:) = 0._EB

   STORE_LOOP : DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

      WC => WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      IIG= BC%IIG; JJG= BC%JJG; KKG= BC%KKG; IOR= BC%IOR

      SELECT CASE(IOR)
      CASE( IAXIS)
         UN_WALLS(IW) = U(IIG-1,JJG  ,KKG  )
      CASE(-IAXIS)
         UN_WALLS(IW) = U(IIG  ,JJG  ,KKG  )
      CASE( JAXIS)
         UN_WALLS(IW) = V(IIG  ,JJG-1,KKG  )
      CASE(-JAXIS)
         UN_WALLS(IW) = V(IIG  ,JJG  ,KKG  )
      CASE( KAXIS)
         UN_WALLS(IW) = W(IIG  ,JJG  ,KKG-1)
      CASE(-KAXIS)
         UN_WALLS(IW) = W(IIG  ,JJG  ,KKG  )
      END SELECT

   ENDDO STORE_LOOP

   RETURN

ENDIF STORE_UN_COND

! Case of not storing, recompute INTERNAL_WALL_CELL velocities, taking into acct that DHDN=0._EB:

PREDICTOR_COND : IF (PREDICTOR) THEN

   WALL_CELL_LOOP_1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

      WC => WALL(IW)
      IF ( .NOT. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR.  &
                  WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) ) CYCLE

      BC => BOUNDARY_COORD(WC%BC_INDEX)
      II  = BC%II; JJ  = BC%JJ; KK  = BC%KK; IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG
      IC  = CELL_INDEX(II,JJ,KK)
      ICG = CELL_INDEX(IIG,JJG,KKG)

      IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. .NOT.CELL(ICG)%SOLID .AND. .NOT.CELL(IC)%SOLID) CYCLE

      IOR = BC%IOR

      SELECT CASE(IOR)
      CASE( IAXIS)
         US(IIG-1,JJG  ,KKG  ) = (U(IIG-1,JJG  ,KKG  ) - DT*( FVX(IIG-1,JJG  ,KKG  ) ))
      CASE(-IAXIS)
         US(IIG  ,JJG  ,KKG  ) = (U(IIG  ,JJG  ,KKG  ) - DT*( FVX(IIG  ,JJG  ,KKG  ) ))
      CASE( JAXIS)
         VS(IIG  ,JJG-1,KKG  ) = (V(IIG  ,JJG-1,KKG  ) - DT*( FVY(IIG  ,JJG-1,KKG  ) ))
      CASE(-JAXIS)
         VS(IIG  ,JJG  ,KKG  ) = (V(IIG  ,JJG  ,KKG  ) - DT*( FVY(IIG  ,JJG  ,KKG  ) ))
      CASE( KAXIS)
         WS(IIG  ,JJG  ,KKG-1) = (W(IIG  ,JJG  ,KKG-1) - DT*( FVZ(IIG  ,JJG  ,KKG-1) ))
      CASE(-KAXIS)
         WS(IIG  ,JJG  ,KKG  ) = (W(IIG  ,JJG  ,KKG  ) - DT*( FVZ(IIG  ,JJG  ,KKG  ) ))
      END SELECT

   ENDDO WALL_CELL_LOOP_1

ELSE ! Corrector

  ! Loop internal wall cells -> on OBST surfaces:

  WALL_CELL_LOOP_2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX

     WC => WALL(IW)
     IF ( .NOT. (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .OR. WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR.  &
                 WC%BOUNDARY_TYPE==MIRROR_BOUNDARY) ) CYCLE

     BC => BOUNDARY_COORD(WC%BC_INDEX)
     II  = BC%II
     JJ  = BC%JJ
     KK  = BC%KK
     IIG = BC%IIG
     JJG = BC%JJG
     KKG = BC%KKG
     IC  = CELL_INDEX(II,JJ,KK)
     ICG = CELL_INDEX(IIG,JJG,KKG)

     IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .AND. .NOT.CELL(ICG)%SOLID .AND. .NOT.CELL(IC)%SOLID) CYCLE

     IOR = BC%IOR
     VEL_N = UN_WALLS(IW)

     SELECT CASE(IOR)
     CASE( IAXIS)                                 
         U(IIG-1,JJG  ,KKG  ) = 0.5_EB*(VEL_N + US(IIG-1,JJG  ,KKG  ) - DT*( FVX(IIG-1,JJG  ,KKG  )  ))
     CASE(-IAXIS)
         U(IIG  ,JJG  ,KKG  ) = 0.5_EB*(VEL_N + US(IIG  ,JJG  ,KKG  ) - DT*( FVX(IIG  ,JJG  ,KKG  )  ))
     CASE( JAXIS)
         V(IIG  ,JJG-1,KKG  ) = 0.5_EB*(VEL_N + VS(IIG  ,JJG-1,KKG  ) - DT*( FVY(IIG  ,JJG-1,KKG  )  ))
     CASE(-JAXIS)
         V(IIG  ,JJG  ,KKG  ) = 0.5_EB*(VEL_N + VS(IIG  ,JJG  ,KKG  ) - DT*( FVY(IIG  ,JJG  ,KKG  )  ))
     CASE( KAXIS)
         W(IIG  ,JJG  ,KKG-1) = 0.5_EB*(VEL_N + WS(IIG  ,JJG  ,KKG-1) - DT*( FVZ(IIG  ,JJG  ,KKG-1)  ))
     CASE(-KAXIS)
         W(IIG  ,JJG  ,KKG  ) = 0.5_EB*(VEL_N + WS(IIG  ,JJG  ,KKG  ) - DT*( FVZ(IIG  ,JJG  ,KKG  )  ))
     END SELECT

  ENDDO WALL_CELL_LOOP_2

  DEALLOCATE(UN_WALLS)

ENDIF PREDICTOR_COND

END SUBROUTINE WALL_VELOCITY_NO_GRADH

END MODULE VELO
