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
       VELOCITY_BC,VELOCITY_BC_PREPROCESSING,&
       COMPUTE_VISCOSITY,VISCOSITY_BC,VELOCITY_FLUX,VELOCITY_FLUX_CYLINDRICAL,&
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

!> \brief Preprocessing for VELOCITY_BC: Transfer velocities from neighboring meshes
!> \param M Mesh data structure
!> \param NM Mesh number
!> \param T Current time (s)
!> \param APPLY_TO_ESTIMATED_VARIABLES Flag indicating that estimated (starred) variables are to be used

SUBROUTINE VELOCITY_BC_PREPROCESSING(M,NM,T,APPLY_TO_ESTIMATED_VARIABLES)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
INTEGER :: IW,IIOO,JJOO,KKOO,N_INT_CELLS
REAL(EB) :: UN_OTHER
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,OM_UU,OM_VV,OM_WW
TYPE(WALL_TYPE), POINTER :: WC
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! Point to the appropriate velocity field

IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
   UU => M%US
   VV => M%VS
   WW => M%WS
ELSE
   UU => M%U
   VV => M%V
   WW => M%W
ENDIF

! Transfer from neighboring mesh the normal component of velocity that is one grid cell beyond external boundary

WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS
   WC =>M%WALL(IW)
   EWC=>M%EXTERNAL_WALL(IW)
   IF (EWC%NOM==0) CYCLE WALL_LOOP
   IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
      OM_UU => M%OMESH(EWC%NOM)%US
      OM_VV => M%OMESH(EWC%NOM)%VS
      OM_WW => M%OMESH(EWC%NOM)%WS
   ELSE
      OM_UU => M%OMESH(EWC%NOM)%U
      OM_VV => M%OMESH(EWC%NOM)%V
      OM_WW => M%OMESH(EWC%NOM)%W
   ENDIF
   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
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

M%DRAG_UVWMAX = 0._EB

END SUBROUTINE VELOCITY_BC_PREPROCESSING

SUBROUTINE VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES)

USE VELO_KERNELS, ONLY: VELOCITY_BC_PROCESS_EDGES_KERNEL
USE CC_VELOCITY, ONLY : CC_VELOCITY_BC

REAL(EB), INTENT(IN) :: T
INTEGER, INTENT(IN) :: NM
LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
REAL(EB) :: T_NOW

IF (SOLID_PHASE_ONLY) RETURN
IF (PERIODIC_TEST==12) RETURN
IF (PERIODIC_TEST==13) RETURN

T_NOW = CURRENT_TIME()

! Preprocessing: OMESH reads for wall boundary velocities
CALL VELOCITY_BC_PREPROCESSING(MESHES(NM),NM,T,APPLY_TO_ESTIMATED_VARIABLES)

! Edge processing kernel: process all edge boundary conditions
CALL VELOCITY_BC_PROCESS_EDGES_KERNEL(MESHES(NM),NM,T,APPLY_TO_ESTIMATED_VARIABLES)

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
