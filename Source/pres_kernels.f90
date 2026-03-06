!> \brief Pure computation kernels extracted from PRES module
!> These routines take TYPE(MESH_TYPE) as an explicit argument instead of relying on MESH_POINTERS.

MODULE PRES_KERNELS

USE PRECISION_PARAMETERS
USE TYPES
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC PRESSURE_SOLVER_COMPUTE_RHS, PRESSURE_SOLVER_FFT, PRESSURE_SOLVER_CHECK_RESIDUALS

CONTAINS


SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS(M,T,DT,NM)

USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: T,DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,HP,RHOP
INTEGER :: I,J,K,IW,IOR,NOM
REAL(EB) :: TRM1,TRM2,TRM3,TRM4, &
            TSI,TIME_RAMP_FACTOR,DX_OTHER,DY_OTHER,DZ_OTHER,P_EXTERNAL,VEL_EDDY,H0
TYPE (VENTS_TYPE), POINTER :: VT
TYPE (WALL_TYPE), POINTER :: WC
TYPE (BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE (BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE (EXTERNAL_WALL_TYPE), POINTER :: EWC

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
   HP => M%H
   RHOP => M%RHO
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
   HP => M%HS
   RHOP => M%RHOS
ENDIF

!$OMP PARALLEL

! Apply pressure boundary conditions at external cells.

!$OMP DO PRIVATE(IW,WC,EWC,BC,B1,I,J,K,IOR,NOM,DX_OTHER,DY_OTHER,DZ_OTHER,VT,TSI) &
!$OMP&   PRIVATE(TIME_RAMP_FACTOR,P_EXTERNAL,VEL_EDDY,H0)
WALL_CELL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS

   WC => M%WALL(IW)
   EWC => M%EXTERNAL_WALL(IW)
   BC => M%BOUNDARY_COORD(WC%BC_INDEX)
   I   = BC%II
   J   = BC%JJ
   K   = BC%KK
   IOR = BC%IOR

   ! Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt

   IF_NEUMANN: IF (EWC%PRESSURE_BC_TYPE==NEUMANN) THEN

      SELECT CASE(IOR)
         CASE( 1)
            M%BXS(J,K) = M%HX(0)      *(-M%FVX(0,J,K)       + EWC%DUNDT)
         CASE(-1)
            M%BXF(J,K) = M%HX(M%IBP1) *(-M%FVX(M%IBAR,J,K)  - EWC%DUNDT)
         CASE( 2)
            M%BYS(I,K) = M%HY(0)      *(-M%FVY(I,0,K)        + EWC%DUNDT)
         CASE(-2)
            M%BYF(I,K) = M%HY(M%JBP1) *(-M%FVY(I,M%JBAR,K)  - EWC%DUNDT)
         CASE( 3)
            M%BZS(I,J) = M%HZ(0)      *(-M%FVZ(I,J,0)        + EWC%DUNDT)
         CASE(-3)
            M%BZF(I,J) = M%HZ(M%KBP1) *(-M%FVZ(I,J,M%KBAR)  - EWC%DUNDT)
      END SELECT
   ENDIF IF_NEUMANN

   ! Apply pressures at DIRICHLET boundaries, depending on the specific type

   IF_DIRICHLET: IF (EWC%PRESSURE_BC_TYPE==DIRICHLET) THEN

      NOT_OPEN: IF (WC%BOUNDARY_TYPE/=OPEN_BOUNDARY .AND. &
                    WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) THEN

         SELECT CASE(IOR)
            CASE( 1) ; M%BXS(J,K) = 0.5_EB*(HP(0,J,K)      +HP(1,J,K))       + M%WALL_WORK1(IW)
            CASE(-1) ; M%BXF(J,K) = 0.5_EB*(HP(M%IBAR,J,K) +HP(M%IBP1,J,K))  + M%WALL_WORK1(IW)
            CASE( 2) ; M%BYS(I,K) = 0.5_EB*(HP(I,0,K)      +HP(I,1,K))       + M%WALL_WORK1(IW)
            CASE(-2) ; M%BYF(I,K) = 0.5_EB*(HP(I,M%JBAR,K) +HP(I,M%JBP1,K))  + M%WALL_WORK1(IW)
            CASE( 3) ; M%BZS(I,J) = 0.5_EB*(HP(I,J,0)      +HP(I,J,1))       + M%WALL_WORK1(IW)
            CASE(-3) ; M%BZF(I,J) = 0.5_EB*(HP(I,J,M%KBAR) +HP(I,J,M%KBP1))  + M%WALL_WORK1(IW)
         END SELECT

      ENDIF NOT_OPEN

      ! Interpolated boundary

      INTERPOLATED_ONLY: IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

         NOM = EWC%NOM

         SELECT CASE(IOR)
            CASE( 1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               M%BXS(J,K) = (DX_OTHER*HP(1,J,K) + M%DX(1)*HP(0,J,K))/(M%DX(1)+DX_OTHER) + M%WALL_WORK1(IW)
            CASE(-1)
               DX_OTHER = MESHES(NOM)%DX(EWC%IIO_MIN)
               M%BXF(J,K) = (DX_OTHER*HP(M%IBAR,J,K) + M%DX(M%IBAR)*HP(M%IBP1,J,K))/ &
                             (M%DX(M%IBAR)+DX_OTHER) + M%WALL_WORK1(IW)
            CASE( 2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               M%BYS(I,K) = (DY_OTHER*HP(I,1,K) + M%DY(1)*HP(I,0,K))/(M%DY(1)+DY_OTHER) + M%WALL_WORK1(IW)
            CASE(-2)
               DY_OTHER = MESHES(NOM)%DY(EWC%JJO_MIN)
               M%BYF(I,K) = (DY_OTHER*HP(I,M%JBAR,K) + M%DY(M%JBAR)*HP(I,M%JBP1,K))/ &
                             (M%DY(M%JBAR)+DY_OTHER) + M%WALL_WORK1(IW)
            CASE( 3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               M%BZS(I,J) = (DZ_OTHER*HP(I,J,1) + M%DZ(1)*HP(I,J,0))/(M%DZ(1)+DZ_OTHER) + M%WALL_WORK1(IW)
            CASE(-3)
               DZ_OTHER = MESHES(NOM)%DZ(EWC%KKO_MIN)
               M%BZF(I,J) = (DZ_OTHER*HP(I,J,M%KBAR) + M%DZ(M%KBAR)*HP(I,J,M%KBP1))/ &
                             (M%DZ(M%KBAR)+DZ_OTHER) + M%WALL_WORK1(IW)
         END SELECT

      ENDIF INTERPOLATED_ONLY

      ! OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.

      OPEN_IF: IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN

         B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
         VT => M%VENTS(WC%VENT_INDEX)
         IF (ABS(B1%T_IGN-T_BEGIN)<=TWENTY_EPSILON_EB .AND. &
             VT%PRESSURE_RAMP_INDEX >=1) THEN
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

         IF (OPEN_WIND_BOUNDARY) THEN
            SELECT CASE(IOR)
               CASE( 1)
                  H0 = HP(1,J,K)       + 0.5_EB/(DT*M%RDXN(0))      *(M%U_WIND(K) + VEL_EDDY - UU(0,J,K))
               CASE(-1)
                  H0 = HP(M%IBAR,J,K)  - 0.5_EB/(DT*M%RDXN(M%IBAR)) *(M%U_WIND(K) + VEL_EDDY - UU(M%IBAR,J,K))
               CASE( 2)
                  H0 = HP(I,1,K)       + 0.5_EB/(DT*M%RDYN(0))      *(M%V_WIND(K) + VEL_EDDY - VV(I,0,K))
               CASE(-2)
                  H0 = HP(I,M%JBAR,K)  - 0.5_EB/(DT*M%RDYN(M%JBAR)) *(M%V_WIND(K) + VEL_EDDY - VV(I,M%JBAR,K))
               CASE( 3)
                  H0 = HP(I,J,1)       + 0.5_EB/(DT*M%RDZN(0))      *(M%W_WIND(K) + VEL_EDDY - WW(I,J,0))
               CASE(-3)
                  H0 = HP(I,J,M%KBAR)  - 0.5_EB/(DT*M%RDZN(M%KBAR)) *(M%W_WIND(K) + VEL_EDDY - WW(I,J,M%KBAR))
            END SELECT
         ENDIF

         SELECT CASE(IOR)
            CASE( 1)
               IF (UU(0,J,K)<0._EB) THEN
                  M%BXS(J,K) = P_EXTERNAL/B1%RHO_F + M%KRES(1,J,K)
               ELSE
                  M%BXS(J,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-1)
               IF (UU(M%IBAR,J,K)>0._EB) THEN
                  M%BXF(J,K) = P_EXTERNAL/B1%RHO_F + M%KRES(M%IBAR,J,K)
               ELSE
                  M%BXF(J,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE( 2)
               IF (VV(I,0,K)<0._EB) THEN
                  M%BYS(I,K) = P_EXTERNAL/B1%RHO_F + M%KRES(I,1,K)
               ELSE
                  M%BYS(I,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-2)
               IF (VV(I,M%JBAR,K)>0._EB) THEN
                  M%BYF(I,K) = P_EXTERNAL/B1%RHO_F + M%KRES(I,M%JBAR,K)
               ELSE
                  M%BYF(I,K) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE( 3)
               IF (WW(I,J,0)<0._EB) THEN
                  M%BZS(I,J) = P_EXTERNAL/B1%RHO_F + M%KRES(I,J,1)
               ELSE
                  M%BZS(I,J) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
            CASE(-3)
               IF (WW(I,J,M%KBAR)>0._EB) THEN
                  M%BZF(I,J) = P_EXTERNAL/B1%RHO_F + M%KRES(I,J,M%KBAR)
               ELSE
                  M%BZF(I,J) = P_EXTERNAL/B1%RHO_F + H0
               ENDIF
         END SELECT

      ENDIF OPEN_IF

   ENDIF IF_DIRICHLET

ENDDO WALL_CELL_LOOP
!$OMP END DO

! Compute the RHS of the Poisson equation

SELECT CASE(M%IPS)

   CASE(:1,4,7)
      IF (CYLINDRICAL) THEN
         !$OMP DO PRIVATE(TRM1,TRM3,TRM4)
         DO K=1,M%KBAR
            DO I=1,M%IBAR
               TRM1 = (M%R(I-1)*M%FVX(I-1,1,K)-M%R(I)*M%FVX(I,1,K))*M%RDX(I)*M%RRN(I)
               TRM3 = (M%FVZ(I,1,K-1)-M%FVZ(I,1,K))*M%RDZ(K)
               TRM4 = -M%DDDT(I,1,K)
               M%PRHS(I,1,K) = TRM1 + TRM3 + TRM4
            ENDDO
         ENDDO
         !$OMP END DO
      ENDIF
      IF (.NOT.CYLINDRICAL) THEN
         !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
         DO K=1,M%KBAR
            DO J=1,M%JBAR
               DO I=1,M%IBAR
                  TRM1 = (M%FVX(I-1,J,K)-M%FVX(I,J,K))*M%RDX(I)
                  TRM2 = (M%FVY(I,J-1,K)-M%FVY(I,J,K))*M%RDY(J)
                  TRM3 = (M%FVZ(I,J,K-1)-M%FVZ(I,J,K))*M%RDZ(K)
                  TRM4 = -M%DDDT(I,J,K)
                  M%PRHS(I,J,K) = TRM1 + TRM2 + TRM3 + TRM4
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO

      ENDIF

   CASE(2)  ! Switch x and y
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               TRM1 = (M%FVX(I-1,J,K)-M%FVX(I,J,K))*M%RDX(I)
               TRM2 = (M%FVY(I,J-1,K)-M%FVY(I,J,K))*M%RDY(J)
               TRM3 = (M%FVZ(I,J,K-1)-M%FVZ(I,J,K))*M%RDZ(K)
               TRM4 = -M%DDDT(I,J,K)
               M%PRHS(J,I,K) = TRM1 + TRM2 + TRM3 + TRM4
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   CASE(3,6)  ! Switch x and z
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               TRM1 = (M%FVX(I-1,J,K)-M%FVX(I,J,K))*M%RDX(I)
               TRM2 = (M%FVY(I,J-1,K)-M%FVY(I,J,K))*M%RDY(J)
               TRM3 = (M%FVZ(I,J,K-1)-M%FVZ(I,J,K))*M%RDZ(K)
               TRM4 = -M%DDDT(I,J,K)
               M%PRHS(K,J,I) = TRM1 + TRM2 + TRM3 + TRM4
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

   CASE(5)  ! Switch y and z
      !$OMP DO PRIVATE(TRM1,TRM2,TRM3,TRM4)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               TRM1 = (M%FVX(I-1,J,K)-M%FVX(I,J,K))*M%RDX(I)
               TRM2 = (M%FVY(I,J-1,K)-M%FVY(I,J,K))*M%RDY(J)
               TRM3 = (M%FVZ(I,J,K-1)-M%FVZ(I,J,K))*M%RDZ(K)
               TRM4 = -M%DDDT(I,J,K)
               M%PRHS(I,K,J) = TRM1 + TRM2 + TRM3 + TRM4
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO

END SELECT

!$OMP END PARALLEL

END SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS


SUBROUTINE PRESSURE_SOLVER_FFT(M,NM)

USE POIS, ONLY: H3CZSS,H2CZSS,H2CYSS,H3CSSS

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP
INTEGER :: I,J,K

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

IF (PREDICTOR) THEN
   HP => M%H
ELSE
   HP => M%HS
ENDIF

! Call the Poisson solver

SELECT CASE(M%IPS)
   CASE(:1)
      IF (.NOT.TWO_D) THEN
         CALL H3CZSS(M%BXS,M%BXF,M%BYS,M%BYF,M%BZS,M%BZF,&
                      M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HX)
      ELSE
         IF (.NOT.CYLINDRICAL) &
            CALL H2CZSS(M%BXS,M%BXF,M%BZS,M%BZF,&
                         M%ITRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HX)
         IF (     CYLINDRICAL) &
            CALL H2CYSS(M%BXS,M%BXF,M%BZS,M%BZF,&
                         M%ITRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK)
      ENDIF
   CASE(2)
      M%BZST = TRANSPOSE(M%BZS)
      M%BZFT = TRANSPOSE(M%BZF)
      CALL H3CZSS(M%BYS,M%BYF,M%BXS,M%BXF,M%BZST,M%BZFT,&
                   M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HY)
   CASE(3)
      IF (.NOT.TWO_D) THEN
         M%BXST = TRANSPOSE(M%BXS)
         M%BXFT = TRANSPOSE(M%BXF)
         M%BYST = TRANSPOSE(M%BYS)
         M%BYFT = TRANSPOSE(M%BYF)
         M%BZST = TRANSPOSE(M%BZS)
         M%BZFT = TRANSPOSE(M%BZF)
         CALL H3CZSS(M%BZST,M%BZFT,M%BYST,M%BYFT,M%BXST,M%BXFT,&
                      M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HZ)
      ELSE
         CALL H2CZSS(M%BZS,M%BZF,M%BXS,M%BXF,&
                      M%ITRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HZ)
      ENDIF
   CASE(4)
      CALL H3CSSS(M%BXS,M%BXF,M%BYS,M%BYF,M%BZS,M%BZF,&
                   M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HX,M%HY)
   CASE(5)
      IF (.NOT.TWO_D) THEN
         M%BXST = TRANSPOSE(M%BXS)
         M%BXFT = TRANSPOSE(M%BXF)
         CALL H3CSSS(M%BXST,M%BXFT,M%BZS,M%BZF,M%BYS,M%BYF,&
                      M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HX,M%HZ)
      ELSE
         CALL H2CZSS(M%BZS,M%BZF,M%BXS,M%BXF,&
                      M%ITRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HZ)
      ENDIF
   CASE(6)
      M%BXST = TRANSPOSE(M%BXS)
      M%BXFT = TRANSPOSE(M%BXF)
      M%BYST = TRANSPOSE(M%BYS)
      M%BYFT = TRANSPOSE(M%BYF)
      M%BZST = TRANSPOSE(M%BZS)
      M%BZFT = TRANSPOSE(M%BZF)
      CALL H3CSSS(M%BZST,M%BZFT,M%BYST,M%BYFT,M%BXST,M%BXFT,&
                   M%ITRN,M%JTRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HZ,M%HY)
   CASE(7)
      CALL H2CZSS(M%BXS,M%BXF,M%BYS,M%BYF,&
                   M%ITRN,M%PRHS,M%POIS_PTB,M%SAVE1,M%WORK,M%HX)
END SELECT

!$OMP PARALLEL

SELECT CASE(M%IPS)
   CASE(:1,4,7)
      !$OMP DO
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(2)
      !$OMP DO
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(J,I,K)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(3,6)
      !$OMP DO
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(K,J,I)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
   CASE(5)
      !$OMP DO
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
END SELECT

! For the special case of tunnels, add back 1-D global pressure solution

IF (TUNNEL_PRECONDITIONER) THEN
   !$OMP MASTER
   DO I=1,M%IBAR
      HP(I,1:M%JBAR,1:M%KBAR) = HP(I,1:M%JBAR,1:M%KBAR) + H_BAR(I_OFFSET(NM)+I)
   ENDDO
   M%BXS = M%BXS + M%BXS_BAR
   M%BXF = M%BXF + M%BXF_BAR
   !$OMP END MASTER
   !$OMP BARRIER
ENDIF

! Apply boundary conditions to H

!$OMP DO
DO K=1,M%KBAR
   DO J=1,M%JBAR
      IF (M%LBC==3 .OR. M%LBC==4)              HP(0,J,K)       = HP(1,J,K)       - M%DXI*M%BXS(J,K)
      IF (M%LBC==3 .OR. M%LBC==2 .OR. M%LBC==6) HP(M%IBP1,J,K) = HP(M%IBAR,J,K) + M%DXI*M%BXF(J,K)
      IF (M%LBC==1 .OR. M%LBC==2)              HP(0,J,K)       =-HP(1,J,K)       + 2._EB*M%BXS(J,K)
      IF (M%LBC==1 .OR. M%LBC==4 .OR. M%LBC==5) HP(M%IBP1,J,K) =-HP(M%IBAR,J,K) + 2._EB*M%BXF(J,K)
      IF (M%LBC==5 .OR. M%LBC==6)              HP(0,J,K)       = HP(1,J,K)
      IF (M%LBC==0) THEN
         HP(0,J,K) = HP(M%IBAR,J,K)
         HP(M%IBP1,J,K) = HP(1,J,K)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO

!$OMP DO
DO K=1,M%KBAR
   DO I=1,M%IBAR
      IF (M%MBC==3 .OR. M%MBC==4) HP(I,0,K)       = HP(I,1,K)       - M%DETA*M%BYS(I,K)
      IF (M%MBC==3 .OR. M%MBC==2) HP(I,M%JBP1,K)  = HP(I,M%JBAR,K) + M%DETA*M%BYF(I,K)
      IF (M%MBC==1 .OR. M%MBC==2) HP(I,0,K)       =-HP(I,1,K)       + 2._EB*M%BYS(I,K)
      IF (M%MBC==1 .OR. M%MBC==4) HP(I,M%JBP1,K)  =-HP(I,M%JBAR,K) + 2._EB*M%BYF(I,K)
      IF (M%MBC==0) THEN
         HP(I,0,K) = HP(I,M%JBAR,K)
         HP(I,M%JBP1,K) = HP(I,1,K)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO

!$OMP DO
DO J=1,M%JBAR
   DO I=1,M%IBAR
      IF (M%NBC==3 .OR. M%NBC==4)  HP(I,J,0)       = HP(I,J,1)       - M%DZETA*M%BZS(I,J)
      IF (M%NBC==3 .OR. M%NBC==2)  HP(I,J,M%KBP1)  = HP(I,J,M%KBAR) + M%DZETA*M%BZF(I,J)
      IF (M%NBC==1 .OR. M%NBC==2)  HP(I,J,0)       =-HP(I,J,1)       + 2._EB*M%BZS(I,J)
      IF (M%NBC==1 .OR. M%NBC==4)  HP(I,J,M%KBP1)  =-HP(I,J,M%KBAR) + 2._EB*M%BZF(I,J)
      IF (M%NBC==0) THEN
         HP(I,J,0) = HP(I,J,M%KBAR)
         HP(I,J,M%KBP1) = HP(I,J,1)
      ENDIF
   ENDDO
ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE PRESSURE_SOLVER_FFT


SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS(M,NM)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP,P,RESIDUAL
INTEGER :: I,J,K
REAL(EB) :: LHSS,RHSS

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

IF (PREDICTOR) THEN
   HP => M%H
   RHOP => M%RHO
ELSE
   HP => M%HS
   RHOP => M%RHOS
ENDIF

! Optional check of the accuracy of the separable pressure solution

IF (CHECK_POISSON) THEN
   RESIDUAL => M%WORK8(1:M%IBAR,1:M%JBAR,1:M%KBAR)
   !$OMP PARALLEL DO PRIVATE(I,J,K,RHSS,LHSS) SCHEDULE(STATIC)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RHSS = ( M%R(I-1)*M%FVX(I-1,J,K) - M%R(I)*M%FVX(I,J,K) )*M%RDX(I)*M%RRN(I) &
                 + (          M%FVY(I,J-1,K) -        M%FVY(I,J,K) )*M%RDY(J)        &
                 + (          M%FVZ(I,J,K-1) -        M%FVZ(I,J,K) )*M%RDZ(K)        &
                 - M%DDDT(I,J,K)
            LHSS = ((HP(I+1,J,K)-HP(I,J,K))*M%RDXN(I)*M%R(I) - &
                    (HP(I,J,K)-HP(I-1,J,K))*M%RDXN(I-1)*M%R(I-1) )*M%RDX(I)*M%RRN(I) &
                 + ((HP(I,J+1,K)-HP(I,J,K))*M%RDYN(J)      - &
                    (HP(I,J,K)-HP(I,J-1,K))*M%RDYN(J-1)        )*M%RDY(J)        &
                 + ((HP(I,J,K+1)-HP(I,J,K))*M%RDZN(K)      - &
                    (HP(I,J,K)-HP(I,J,K-1))*M%RDZN(K-1)        )*M%RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END PARALLEL DO
   M%POIS_ERR = MAXVAL(RESIDUAL)
ENDIF

! Mandatory check of inseparable Poisson equation

IF (ITERATE_BAROCLINIC_TERM) THEN

   P => M%WORK7
   RESIDUAL => M%WORK8(1:M%IBAR,1:M%JBAR,1:M%KBAR)

   !$OMP PARALLEL

   !$OMP DO SCHEDULE(STATIC)
   DO K=0,M%KBP1
      DO J=0,M%JBP1
         DO I=0,M%IBP1
            P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-M%KRES(I,J,K))
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) PRIVATE(I,J,K,RHSS,LHSS)
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RHSS = ( M%R(I-1)*(M%FVX(I-1,J,K)-M%FVX_B(I-1,J,K)) - &
                     M%R(I)  *(M%FVX(I,J,K)  -M%FVX_B(I,J,K)) )*M%RDX(I)*M%RRN(I) &
                 + ( (M%FVY(I,J-1,K)-M%FVY_B(I,J-1,K)) - &
                     (M%FVY(I,J,K)  -M%FVY_B(I,J,K)) )*M%RDY(J) &
                 + ( (M%FVZ(I,J,K-1)-M%FVZ_B(I,J,K-1)) - &
                     (M%FVZ(I,J,K)  -M%FVZ_B(I,J,K)) )*M%RDZ(K) &
                 - M%DDDT(I,J,K)
            LHSS = ((P(I+1,J,K)-P(I,J,K))*M%RDXN(I)*M%R(I) &
                    *2._EB/(RHOP(I+1,J,K)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I-1,J,K))*M%RDXN(I-1)*M%R(I-1) &
                    *2._EB/(RHOP(I-1,J,K)+RHOP(I,J,K)))*M%RDX(I)*M%RRN(I) &
                 + ((P(I,J+1,K)-P(I,J,K))*M%RDYN(J) &
                    *2._EB/(RHOP(I,J+1,K)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I,J-1,K))*M%RDYN(J-1) &
                    *2._EB/(RHOP(I,J-1,K)+RHOP(I,J,K)))*M%RDY(J) &
                 + ((P(I,J,K+1)-P(I,J,K))*M%RDZN(K) &
                    *2._EB/(RHOP(I,J,K+1)+RHOP(I,J,K)) - &
                    (P(I,J,K)-P(I,J,K-1))*M%RDZN(K-1) &
                    *2._EB/(RHOP(I,J,K-1)+RHOP(I,J,K)))*M%RDZ(K) &
                 + ((M%KRES(I+1,J,K)-M%KRES(I,J,K))*M%RDXN(I)*M%R(I) - &
                    (M%KRES(I,J,K)-M%KRES(I-1,J,K))*M%RDXN(I-1)*M%R(I-1) &
                    )*M%RDX(I)*M%RRN(I) &
                 + ((M%KRES(I,J+1,K)-M%KRES(I,J,K))*M%RDYN(J) - &
                    (M%KRES(I,J,K)-M%KRES(I,J-1,K))*M%RDYN(J-1) &
                    )*M%RDY(J) &
                 + ((M%KRES(I,J,K+1)-M%KRES(I,J,K))*M%RDZN(K) - &
                    (M%KRES(I,J,K)-M%KRES(I,J,K-1))*M%RDZN(K-1) &
                    )*M%RDZ(K)
            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
         ENDDO
      ENDDO
   ENDDO
   !$OMP END DO

   !$OMP END PARALLEL

   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
   IF (STORE_PRESSURE_POISSON_RESIDUAL) &
      M%PP_RESIDUAL(1:M%IBAR,1:M%JBAR,1:M%KBAR) = &
         RESIDUAL(1:M%IBAR,1:M%JBAR,1:M%KBAR)

ENDIF

END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS


END MODULE PRES_KERNELS
