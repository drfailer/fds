!> \brief Pure computation kernels extracted from PRES module
!> These routines take TYPE(MESH_TYPE) as an explicit argument instead of relying on MESH_POINTERS.

MODULE PRES_KERNELS

USE PRECISION_PARAMETERS
USE TYPES
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC PRESSURE_SOLVER_COMPUTE_RHS, PRESSURE_SOLVER_FFT, PRESSURE_SOLVER_CHECK_RESIDUALS, COMPUTE_VELOCITY_ERROR_KERNEL

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


! Apply pressure boundary conditions at external cells.

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

! Compute the RHS of the Poisson equation

SELECT CASE(M%IPS)

   CASE(:1,4,7)
      IF (CYLINDRICAL) THEN
         DO K=1,M%KBAR
            DO I=1,M%IBAR
               TRM1 = (M%R(I-1)*M%FVX(I-1,1,K)-M%R(I)*M%FVX(I,1,K))*M%RDX(I)*M%RRN(I)
               TRM3 = (M%FVZ(I,1,K-1)-M%FVZ(I,1,K))*M%RDZ(K)
               TRM4 = -M%DDDT(I,1,K)
               M%PRHS(I,1,K) = TRM1 + TRM3 + TRM4
            ENDDO
         ENDDO
      ENDIF
      IF (.NOT.CYLINDRICAL) THEN
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

      ENDIF

   CASE(2)  ! Switch x and y
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

   CASE(3,6)  ! Switch x and z
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

   CASE(5)  ! Switch y and z
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

END SELECT


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


SELECT CASE(M%IPS)
   CASE(:1,4,7)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(I,J,K)
            ENDDO
         ENDDO
      ENDDO
   CASE(2)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(J,I,K)
            ENDDO
         ENDDO
      ENDDO
   CASE(3,6)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(K,J,I)
            ENDDO
         ENDDO
      ENDDO
   CASE(5)
      DO K=1,M%KBAR
         DO J=1,M%JBAR
            DO I=1,M%IBAR
               HP(I,J,K) = M%PRHS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
END SELECT

! For the special case of tunnels, add back 1-D global pressure solution

IF (TUNNEL_PRECONDITIONER) THEN
   DO I=1,M%IBAR
      HP(I,1:M%JBAR,1:M%KBAR) = HP(I,1:M%JBAR,1:M%KBAR) + H_BAR(I_OFFSET(NM)+I)
   ENDDO
   M%BXS = M%BXS + M%BXS_BAR
   M%BXF = M%BXF + M%BXF_BAR
ENDIF

! Apply boundary conditions to H

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
   M%POIS_ERR = MAXVAL(RESIDUAL)
ENDIF

! Mandatory check of inseparable Poisson equation

IF (ITERATE_BAROCLINIC_TERM) THEN

   P => M%WORK7
   RESIDUAL => M%WORK8(1:M%IBAR,1:M%JBAR,1:M%KBAR)


   DO K=0,M%KBP1
      DO J=0,M%JBP1
         DO I=0,M%IBP1
            P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-M%KRES(I,J,K))
         ENDDO
      ENDDO
   ENDDO

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


   PRESSURE_ERROR_MAX(NM) = MAXVAL(RESIDUAL)
   PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC(RESIDUAL)
   IF (STORE_PRESSURE_POISSON_RESIDUAL) &
      M%PP_RESIDUAL(1:M%IBAR,1:M%JBAR,1:M%KBAR) = &
         RESIDUAL(1:M%IBAR,1:M%JBAR,1:M%KBAR)

ENDIF

END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS


!> \brief Compute velocity error at solid and interpolated boundaries.
!> \details Thread-safe kernel version of COMPUTE_VELOCITY_ERROR (pres.f90).
!> Compares predicted velocity with boundary conditions and neighbor mesh velocities.
!> Writes VELOCITY_ERROR_MAX(NM), VELOCITY_ERROR_MAX_LOC(:,NM), and WALL_WORK1.
!> \param M Mesh data structure
!> \param DT Time step (s)
!> \param NM Mesh index

SUBROUTINE COMPUTE_VELOCITY_ERROR_KERNEL(M,DT,NM)

USE COMPLEX_GEOMETRY, ONLY: CC_CGSC,CC_GASPHASE

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: NM
INTEGER :: IW,IOR,II,JJ,KK,IIO,JJO,KKO,N_INT_CELLS,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2
REAL(EB) :: UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DUDT,DVDT,DWDT,ITERATIVE_FACTOR,DHFCT
TYPE(OMESH_TYPE), POINTER :: OM
TYPE(MESH_TYPE), POINTER :: M2
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(EXTERNAL_WALL_TYPE), POINTER :: EWC

IF (SOLID_PHASE_ONLY) RETURN
IF (FREEZE_VELOCITY)  RETURN

IF (PREDICTOR) THEN
   ITERATIVE_FACTOR = 0.25_EB
ELSE
   ITERATIVE_FACTOR = 0.50_EB
ENDIF

VELOCITY_ERROR_MAX(NM) = 0._EB
M%WALL_WORK1 = 0._EB

! Loop over wall cells and check velocity error.

CHECK_WALL_LOOP: DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS

   WC=>M%WALL(IW)

   IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY        .AND. &
       WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE CHECK_WALL_LOOP

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
      EWC=>M%EXTERNAL_WALL(IW)
      IF (EWC%AREA_RATIO<0.9_EB) CYCLE CHECK_WALL_LOOP
      OM => M%OMESH(EWC%NOM)
      M2 => MESHES(EWC%NOM)
   ENDIF

   B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
   BC => M%BOUNDARY_COORD(WC%BC_INDEX)

   II  = BC%II
   JJ  = BC%JJ
   KK  = BC%KK
   IOR = BC%IOR

   IF (CC_IBM) THEN
      IF (ANY((/M%CCVAR(BC%IIG,BC%JJG,BC%KKG,CC_CGSC),M%CCVAR(II,JJ,KK,CC_CGSC)/)/=CC_GASPHASE)) CYCLE CHECK_WALL_LOOP
   ENDIF

   DHFCT = 1._EB
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      SELECT CASE(PRES_FLAG)
      CASE(UGLMAT_FLAG,ULMAT_FLAG); DHFCT=0._EB
      CASE(GLMAT_FLAG); IF (IW<=M%N_EXTERNAL_WALL_CELLS) DHFCT=0._EB
      END SELECT
   ENDIF

   ! Update normal component of velocity at the mesh boundary

   IF (PREDICTOR) THEN
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = M%U(II,JJ,KK)   - DT*(M%FVX(II,JJ,KK)   + M%RDXN(II)  *(M%H(II+1,JJ,KK)-M%H(II,JJ,KK))*DHFCT)
         CASE(-1)
            UN_NEW = M%U(II-1,JJ,KK) - DT*(M%FVX(II-1,JJ,KK) + M%RDXN(II-1)*(M%H(II,JJ,KK)-M%H(II-1,JJ,KK))*DHFCT)
         CASE( 2)
            UN_NEW = M%V(II,JJ,KK)   - DT*(M%FVY(II,JJ,KK)   + M%RDYN(JJ)  *(M%H(II,JJ+1,KK)-M%H(II,JJ,KK))*DHFCT)
         CASE(-2)
            UN_NEW = M%V(II,JJ-1,KK) - DT*(M%FVY(II,JJ-1,KK) + M%RDYN(JJ-1)*(M%H(II,JJ,KK)-M%H(II,JJ-1,KK))*DHFCT)
         CASE( 3)
            UN_NEW = M%W(II,JJ,KK)   - DT*(M%FVZ(II,JJ,KK)   + M%RDZN(KK)  *(M%H(II,JJ,KK+1)-M%H(II,JJ,KK))*DHFCT)
         CASE(-3)
            UN_NEW = M%W(II,JJ,KK-1) - DT*(M%FVZ(II,JJ,KK-1) + M%RDZN(KK-1)*(M%H(II,JJ,KK)-M%H(II,JJ,KK-1))*DHFCT)
      END SELECT
   ELSE
      SELECT CASE(IOR)
         CASE( 1)
            UN_NEW = 0.5_EB*(M%U(II,JJ,KK)+M%US(II,JJ,KK) &
                     - DT*(M%FVX(II,JJ,KK)+M%RDXN(II)*(M%HS(II+1,JJ,KK)-M%HS(II,JJ,KK))*DHFCT))
         CASE(-1)
            UN_NEW = 0.5_EB*(M%U(II-1,JJ,KK)+M%US(II-1,JJ,KK) &
                     - DT*(M%FVX(II-1,JJ,KK)+M%RDXN(II-1)*(M%HS(II,JJ,KK)-M%HS(II-1,JJ,KK))*DHFCT))
         CASE( 2)
            UN_NEW = 0.5_EB*(M%V(II,JJ,KK)+M%VS(II,JJ,KK) &
                     - DT*(M%FVY(II,JJ,KK)+M%RDYN(JJ)*(M%HS(II,JJ+1,KK)-M%HS(II,JJ,KK))*DHFCT))
         CASE(-2)
            UN_NEW = 0.5_EB*(M%V(II,JJ-1,KK)+M%VS(II,JJ-1,KK) &
                     - DT*(M%FVY(II,JJ-1,KK)+M%RDYN(JJ-1)*(M%HS(II,JJ,KK)-M%HS(II,JJ-1,KK))*DHFCT))
         CASE( 3)
            UN_NEW = 0.5_EB*(M%W(II,JJ,KK)+M%WS(II,JJ,KK) &
                     - DT*(M%FVZ(II,JJ,KK)+M%RDZN(KK)*(M%HS(II,JJ,KK+1)-M%HS(II,JJ,KK))*DHFCT))
         CASE(-3)
            UN_NEW = 0.5_EB*(M%W(II,JJ,KK-1)+M%WS(II,JJ,KK-1) &
                     - DT*(M%FVZ(II,JJ,KK-1)+M%RDZN(KK-1)*(M%HS(II,JJ,KK)-M%HS(II,JJ,KK-1))*DHFCT))
      END SELECT
   ENDIF

   ! At interpolated boundaries, compare updated normal component of velocity with that of the other mesh

   IF (WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN

      UN_NEW_OTHER = 0._EB

      EWC=>M%EXTERNAL_WALL(IW)
      IIO1 = EWC%IIO_MIN
      JJO1 = EWC%JJO_MIN
      KKO1 = EWC%KKO_MIN
      IIO2 = EWC%IIO_MAX
      JJO2 = EWC%JJO_MAX
      KKO2 = EWC%KKO_MAX

      PREDICTOR_IF: IF (PREDICTOR) THEN
         IOR_SELECT_1: SELECT CASE(IOR)
            CASE( 1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%H(IIO+1,JJO,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO,JJO,KKO)   + DT*DUDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%U(IIO-1,JJO,KKO) + DT*DUDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%H(IIO,JJO+1,KKO)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO,KKO)   + DT*DVDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%V(IIO,JJO-1,KKO) + DT*DVDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%H(IIO,JJO,KKO+1)-OM%H(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO)   + DT*DWDT
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%H(IIO,JJO,KKO)-OM%H(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + OM%W(IIO,JJO,KKO-1) + DT*DWDT
                     ENDDO
                  ENDDO
               ENDDO
         END SELECT IOR_SELECT_1
      ELSE PREDICTOR_IF
         IOR_SELECT_2: SELECT CASE(IOR)
            CASE( 1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO,JJO,KKO)   - M2%RDXN(IIO)  *(OM%HS(IIO+1,JJO,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO,JJO,KKO)+OM%US(IIO,JJO,KKO)     + DT*DUDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-1)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DUDT = -OM%FVX(IIO-1,JJO,KKO) - M2%RDXN(IIO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO-1,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%U(IIO-1,JJO,KKO)+OM%US(IIO-1,JJO,KKO) + DT*DUDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO,KKO)   - M2%RDYN(JJO)  *(OM%HS(IIO,JJO+1,KKO)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO,KKO)+OM%VS(IIO,JJO,KKO)     + DT*DVDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-2)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DVDT = -OM%FVY(IIO,JJO-1,KKO) - M2%RDYN(JJO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO-1,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%V(IIO,JJO-1,KKO)+OM%VS(IIO,JJO-1,KKO) + DT*DVDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE( 3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO)   - M2%RDZN(KKO)  *(OM%HS(IIO,JJO,KKO+1)-OM%HS(IIO,JJO,KKO))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO)+OM%WS(IIO,JJO,KKO)     + DT*DWDT)
                     ENDDO
                  ENDDO
               ENDDO
            CASE(-3)
               DO KKO=KKO1,KKO2
                  DO JJO=JJO1,JJO2
                     DO IIO=IIO1,IIO2
                        DWDT = -OM%FVZ(IIO,JJO,KKO-1) - M2%RDZN(KKO-1)*(OM%HS(IIO,JJO,KKO)-OM%HS(IIO,JJO,KKO-1))
                        UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM%W(IIO,JJO,KKO-1)+OM%WS(IIO,JJO,KKO-1) + DT*DWDT)
                     ENDDO
                  ENDDO
               ENDDO
         END SELECT IOR_SELECT_2
      ENDIF PREDICTOR_IF

      N_INT_CELLS  = (EWC%IIO_MAX-EWC%IIO_MIN+1) * (EWC%JJO_MAX-EWC%JJO_MIN+1) * (EWC%KKO_MAX-EWC%KKO_MIN+1)
      UN_NEW_OTHER = UN_NEW_OTHER/REAL(N_INT_CELLS,EB)

   ENDIF

   ! At solid boundaries, compare updated normal velocity with specified normal velocity

   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      IF (PREDICTOR) THEN
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
      ELSE
         UN_NEW_OTHER = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      ENDIF
   ENDIF

   ! Compute velocity difference

   VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
   B1%VEL_ERR_NEW = VELOCITY_ERROR
   M%WALL_WORK1(IW) = -SIGN(1._EB,REAL(IOR,EB))*ITERATIVE_FACTOR*VELOCITY_ERROR/(B1%RDN*DT)

   ! Save maximum velocity error

   IF (ABS(VELOCITY_ERROR)>VELOCITY_ERROR_MAX(NM)) THEN
      VELOCITY_ERROR_MAX_LOC(1,NM) = II
      VELOCITY_ERROR_MAX_LOC(2,NM) = JJ
      VELOCITY_ERROR_MAX_LOC(3,NM) = KK
      SELECT CASE(IOR)
         CASE(-1) ; VELOCITY_ERROR_MAX_LOC(1,NM) = II-1
         CASE(-2) ; VELOCITY_ERROR_MAX_LOC(2,NM) = JJ-1
         CASE(-3) ; VELOCITY_ERROR_MAX_LOC(3,NM) = KK-1
      END SELECT
      VELOCITY_ERROR_MAX(NM)       = ABS(VELOCITY_ERROR)
   ENDIF

ENDDO CHECK_WALL_LOOP

END SUBROUTINE COMPUTE_VELOCITY_ERROR_KERNEL


END MODULE PRES_KERNELS
