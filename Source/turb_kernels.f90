!> \brief Pure computation kernels extracted from TURBULENCE module.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE TURB_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC WALE_VISCOSITY,WALL_MODEL,TAU_WALL_IJ,TEST_FILTER_LOCAL, &
       EX2G3D_KERNEL,FILL_EDGES_KERNEL,TEST_FILTER_KERNEL, &
       TENSOR_DIFFUSIVITY_MODEL_KERNEL, &
       VARDEN_DYNSMAG_KERNEL

CONTAINS


!> \brief Wall Adapting Local Eddy-viscosity (WALE) model.
!> \param NU_T Eddy viscosity (output)
!> \param G_IJ Velocity gradient tensor (3x3)
!> \param DELTA Filter width

SUBROUTINE WALE_VISCOSITY(NU_T,G_IJ,DELTA)

REAL(EB), INTENT(OUT) :: NU_T
REAL(EB), INTENT(IN) :: G_IJ(3,3),DELTA
REAL(EB) :: S_IJ(3,3),O_IJ(3,3),S2,O2,IV_SO,SD2,DENOM
INTEGER :: I,J,K,L

! compute strain and rotation tensors

DO J=1,3
   DO I=1,3
      S_IJ(I,J) = 0.5_EB * ( G_IJ(I,J) + G_IJ(J,I) )
      O_IJ(I,J) = 0.5_EB * ( G_IJ(I,J) - G_IJ(J,I) )
   ENDDO
ENDDO

! contraction of strain and rotation tensors

S2 = 0._EB
O2 = 0._EB
DO J=1,3
   DO I=1,3
      S2 = S2 + S_IJ(I,J)*S_IJ(I,J)
      O2 = O2 + O_IJ(I,J)*O_IJ(I,J)
   ENDDO
ENDDO

! fourth order contraction

IV_SO = 0._EB
DO L=1,3
   DO K=1,3
      DO J=1,3
         DO I=1,3
            IV_SO = IV_SO + S_IJ(I,K)*S_IJ(K,J)*O_IJ(J,L)*O_IJ(L,I)
         ENDDO
      ENDDO
   ENDDO
ENDDO

! using Caley-Hamilton theorem

SD2 = ONSI*(S2*S2 + O2*O2) + TWTH*S2*O2 + 2._EB*IV_SO
IF (SD2 < 0.0) SD2 = 0._EB

DENOM = S2**2.5_EB + SD2**1.25_EB
IF (DENOM>TWENTY_EPSILON_EB) THEN
   NU_T = (C_WALE*DELTA)**2 * SD2**1.5_EB / DENOM
ELSE
   NU_T = 0._EB
ENDIF

END SUBROUTINE WALE_VISCOSITY


!> \brief Compute wall model slip factor and friction velocity.

SUBROUTINE WALL_MODEL(SLIP_FACTOR,U_TAU,Y_PLUS,NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT,&
                      Y_FORCING_POINT,U_FORCING_POINT,DUDY_FORCING_POINT)

REAL(EB), INTENT(OUT) :: SLIP_FACTOR,U_TAU,Y_PLUS
! S is the "sandgrain" roughness length scale (Pope's notation)
! Y_EXTERNAL_POINT is the distance from the wall to the "external point"; DN/2 for Cartesian grids
! Y_FORCING_POINT  is the distance from the wall to the "forcing point" where the forced velocity lives
REAL(EB), INTENT(IN) :: NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT
REAL(EB), OPTIONAL, INTENT(IN) :: Y_FORCING_POINT
REAL(EB), OPTIONAL, INTENT(OUT) :: U_FORCING_POINT,DUDY_FORCING_POINT

REAL(EB), PARAMETER :: B=5.2_EB,BTILDE_MAX=9.5_EB ! BTILDE_ROUGH=8.5 set in GLOBAL_CONSTANTS; see Pope (2000) pp. 294,297,298
REAL(EB), PARAMETER :: S0=1._EB,S1=5.83_EB,S2=30._EB ! approx piece-wise function for Fig. 7.24, Pope (2000) p. 297
REAL(EB), PARAMETER :: EPS=1.E-10_EB
INTEGER, PARAMETER :: LAMINAR_SMOOTH=1,TURBULENT_SMOOTH=2,TURBULENT_ROUGH=3

REAL(EB) :: U,Y_CELL_CENTER,TAU_W,BTILDE,DELTA_NU,S_PLUS,DUDY,DY,RKAPPA
INTEGER :: ITER,BOUNDARY_LAYER_CODE

! References:
!
! S. B. Pope (2000) Turbulent Flows, Cambridge.

! Step 1: compute laminar (DNS) stress, and initial guess for LES stress

RKAPPA = 1._EB/VON_KARMAN_CONSTANT
DY = 2._EB*Y_EXTERNAL_POINT
Y_CELL_CENTER = Y_EXTERNAL_POINT
U = U_EXTERNAL_POINT
DUDY = ABS(U)/Y_CELL_CENTER
TAU_W = NU*DUDY                         ! actually tau_w/rho
U_TAU = SQRT(ABS(TAU_W))                ! friction velocity
DELTA_NU = NU/(U_TAU+EPS)               ! viscous length scale
Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
SLIP_FACTOR = -1._EB
BOUNDARY_LAYER_CODE = LAMINAR_SMOOTH

! Step 2: compute turbulent (LES) stress

LES_IF: IF (SIM_MODE/=DNS_MODE) THEN

   ! NOTE: 2 iterations converges TAU_W to roughly 5 % residual error
   !       3 iterations converges TAU_W to roughly 1 % residual error

   DO ITER=1,3

      S_PLUS = S/(DELTA_NU+EPS) ! roughness in viscous units

      IF (S_PLUS < S0) THEN
         ! smooth wall
         Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
         IF (Y_PLUS < Y_WERNER_WENGLE) THEN
            ! viscous sublayer
            TAU_W = ( U/Y_PLUS )**2
            U_TAU = SQRT(TAU_W)
            DUDY = ABS(U)/Y_CELL_CENTER
            BOUNDARY_LAYER_CODE = LAMINAR_SMOOTH
         ELSE
            ! log layer
            TAU_W = ( U/(RKAPPA*LOG(Y_PLUS)+B) )**2
            U_TAU = SQRT(TAU_W)
            DUDY = U_TAU*RKAPPA/Y_CELL_CENTER
            BOUNDARY_LAYER_CODE = TURBULENT_SMOOTH
         ENDIF
      ELSE
         ! rough wall
         IF (S_PLUS < S1) THEN
            BTILDE = B + RKAPPA*LOG(S_PLUS) ! Pope (2000) p. 297, Eq. (7.122)
         ELSE IF (S_PLUS < S2) THEN
            BTILDE = BTILDE_MAX ! approximation from Fig. 7.24, Pope (2000) p. 297
         ELSE
            BTILDE = BTILDE_ROUGH ! fully rough
         ENDIF
         Y_PLUS = Y_CELL_CENTER/S
         TAU_W = ( U/(RKAPPA*LOG(Y_PLUS)+BTILDE) )**2  ! Pope (2000) p. 297, Eq. (7.121)
         U_TAU = SQRT(TAU_W)
         DUDY = U_TAU*RKAPPA/Y_CELL_CENTER
         BOUNDARY_LAYER_CODE = TURBULENT_ROUGH
      ENDIF

      DELTA_NU = NU/(U_TAU+EPS)

   ENDDO

   ! NOTE: SLIP_FACTOR is no longer used to compute the wall stress, see VELOCITY_BC.
   ! The stress is taken directly from U_TAU. SLIP_FACTOR is, however, still used to
   ! compute the velocity gradient at the wall that feeds into the wall vorticity.
   ! Since the gradients implied by the wall function can be large and lead to instabilities,
   ! we bound the wall slip between no slip (-1) and free slip (1).

   ! The slip factor (SF) is based on the following approximation to the wall gradient
   ! (note that u0 is the ghost cell value of the streamwise velocity component and
   ! y is the wall-normal direction):
   ! DUDY = (u-u0)/dy = (u-SF*u)/dy = u/dy*(1-SF) => SF = 1 - DUDY*dy/u
   ! In this routine, DUDY is sampled from the wall model at the location y_cell_center.

   SLIP_FACTOR = MAX(-1._EB,MIN(1._EB,1._EB-DUDY*DY/(ABS(U)+EPS))) ! -1.0 <= SLIP_FACTOR <= 1.0

ENDIF LES_IF

! complex geometry

IF (PRESENT(Y_FORCING_POINT)) THEN
   SELECT CASE(BOUNDARY_LAYER_CODE)
      CASE(LAMINAR_SMOOTH)
         DUDY_FORCING_POINT = U_TAU/DELTA_NU
         U_FORCING_POINT = DUDY_FORCING_POINT * Y_FORCING_POINT
      CASE(TURBULENT_SMOOTH)
         DUDY_FORCING_POINT = U_TAU * RKAPPA / Y_FORCING_POINT
         U_FORCING_POINT = U_TAU * (RKAPPA * LOG(Y_FORCING_POINT/DELTA_NU) + B)
      CASE(TURBULENT_ROUGH)
         DUDY_FORCING_POINT = U_TAU * RKAPPA / Y_FORCING_POINT
         U_FORCING_POINT = U_TAU * (RKAPPA * LOG(Y_FORCING_POINT/S) + BTILDE)
   END SELECT
ENDIF

END SUBROUTINE WALL_MODEL


!> \brief Compute wall stress tensor in Cartesian coordinates.

SUBROUTINE TAU_WALL_IJ(TAU_IJ,SS,U_VELO,U_SURF,NN,DN,DIVU,MU,RHO,ROUGHNESS)
USE MATH_FUNCTIONS, ONLY: CROSS_PRODUCT

REAL(EB), INTENT(OUT) :: TAU_IJ(3,3),SS(3)
REAL(EB), INTENT(IN) :: U_VELO(3),U_SURF(3),NN(3),DN,DIVU,MU,RHO,ROUGHNESS
REAL(EB) :: C(3,3),TT(3),U_RELA(3),SLIP_COEF,Y_PLUS,U_STRM,U_ORTH,U_NORM,TAUBAR_IJ(3,3),U_TAU
INTEGER :: K,L,M,N

! Cartesian grid coordinate system orthonormal basis vectors
REAL(EB), DIMENSION(3), PARAMETER :: E1=(/1._EB,0._EB,0._EB/),E2=(/0._EB,1._EB,0._EB/),E3=(/0._EB,0._EB,1._EB/)

! find a vector TT in the tangent plane of the surface and orthogonal to U_VELO-U_SURF
U_RELA = U_VELO-U_SURF
CALL CROSS_PRODUCT(TT,NN,U_RELA) ! TT = NN x U_RELA
IF (ABS(NORM2(TT))<=TWENTY_EPSILON_EB) THEN
   ! tangent vector is completely arbitrary, just perpendicular to NN
   IF (ABS(NN(1))>=TWENTY_EPSILON_EB .OR.  ABS(NN(2))>=TWENTY_EPSILON_EB) TT = (/NN(2),-NN(1),0._EB/)
   IF (ABS(NN(1))<=TWENTY_EPSILON_EB .AND. ABS(NN(2))<=TWENTY_EPSILON_EB) TT = (/NN(3),0._EB,-NN(1)/)
ENDIF
TT = TT/NORM2(TT) ! normalize to unit vector
CALL CROSS_PRODUCT(SS,TT,NN) ! define the streamwise unit vector SS

! directional cosines (see Pope, Eq. A.11)
C(1,1) = DOT_PRODUCT(E1,SS)
C(1,2) = DOT_PRODUCT(E1,TT)
C(1,3) = DOT_PRODUCT(E1,NN)
C(2,1) = DOT_PRODUCT(E2,SS)
C(2,2) = DOT_PRODUCT(E2,TT)
C(2,3) = DOT_PRODUCT(E2,NN)
C(3,1) = DOT_PRODUCT(E3,SS)
C(3,2) = DOT_PRODUCT(E3,TT)
C(3,3) = DOT_PRODUCT(E3,NN)

! transform velocity (see Pope, Eq. A.17)
U_STRM = C(1,1)*U_RELA(1) + C(2,1)*U_RELA(2) + C(3,1)*U_RELA(3)
U_ORTH = C(1,2)*U_RELA(1) + C(2,2)*U_RELA(2) + C(3,2)*U_RELA(3)
U_NORM = C(1,3)*U_RELA(1) + C(2,3)*U_RELA(2) + C(3,3)*U_RELA(3)

! in the streamwise coordinate system, the stress tensor simplifies to the symmetric tensor
! T = [0      0 T(1,3)]
!     [0      0      0]
!     [T(3,1) 0 T(3,3)]

TAUBAR_IJ      = 0._EB
TAUBAR_IJ(3,3) = -2._EB*MU*(U_NORM*2._EB/DN - ONTH*DIVU)

IF (SIM_MODE==DNS_MODE) THEN
   TAUBAR_IJ(1,3) = MU*U_STRM*2._EB/DN
ELSE
   CALL WALL_MODEL(SLIP_COEF,U_TAU,Y_PLUS,MU/RHO,ROUGHNESS,0.5_EB*DN,U_STRM)
   TAUBAR_IJ(1,3) = RHO*U_TAU**2
ENDIF
TAUBAR_IJ(3,1) = TAUBAR_IJ(1,3)

! transform tensors (Pope A.23)
TAU_IJ = 0._EB
DO M=1,3
   DO N=1,3
      ! inner summation for component m,n ---------------------
      DO L=1,3
         DO K=1,3
            TAU_IJ(M,N) = TAU_IJ(M,N) + C(M,K)*C(N,L)*TAUBAR_IJ(K,L)
         ENDDO
      ENDDO
      !--------------------------------------------------------
   ENDDO
ENDDO

RETURN
END SUBROUTINE TAU_WALL_IJ


!> \brief Local box filter of width 2*dx with trapezoid quadrature.

SUBROUTINE TEST_FILTER_LOCAL(HAT,ORIG)

REAL(EB), INTENT(IN) :: ORIG(-1:1,-1:1,-1:1)
REAL(EB), INTENT(OUT) :: HAT
INTEGER :: I, J, K, L, M, N
REAL(EB), PARAMETER :: K1DT(3)=(/1.0_EB,2.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DT(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DT(I)*K1DT(J)*K1DT(K)/64._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))

! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
HAT = 0._EB
DO N = -1,1
   DO M = -1,1
      DO L = -1,1
         HAT = HAT + ORIG(L,M,N) * K3DT(L,M,N)
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE TEST_FILTER_LOCAL


!> \brief Second order extrapolation of 3D array to ghost cells.
!> \param M Mesh data structure
!> \param A Array to extrapolate
!> \param A_MIN Minimum clamp value
!> \param A_MAX Maximum clamp value

SUBROUTINE EX2G3D_KERNEL(M,A,A_MIN,A_MAX)

TYPE(MESH_TYPE), INTENT(IN), TARGET :: M
REAL(EB), INTENT(IN) :: A_MIN,A_MAX
REAL(EB), INTENT(INOUT) :: A(0:M%IBP1,0:M%JBP1,0:M%KBP1)
INTEGER :: IBP1,JBP1,KBP1,IBAR,JBAR,KBAR,IBM1,JBM1,KBM1

IBP1 = M%IBP1; JBP1 = M%JBP1; KBP1 = M%KBP1
IBAR = M%IBAR; JBAR = M%JBAR; KBAR = M%KBAR
IBM1 = M%IBAR-1; JBM1 = M%JBAR-1; KBM1 = M%KBAR-1

A(0,:,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(1,:,:)-A(2,:,:)))
A(:,0,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,1,:)-A(:,2,:)))
A(:,:,0) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,:,1)-A(:,:,2)))

A(IBP1,:,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(IBAR,:,:)-A(IBM1,:,:)))
A(:,JBP1,:) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,JBAR,:)-A(:,JBM1,:)))
A(:,:,KBP1) = MIN(A_MAX,MAX(A_MIN,2._EB*A(:,:,KBAR)-A(:,:,KBM1)))

END SUBROUTINE EX2G3D_KERNEL


!> \brief Extrapolate array A to edges and corners of 3D arrays.
!> \param M Mesh data structure
!> \param A Array to fill

SUBROUTINE FILL_EDGES_KERNEL(M,A)

TYPE(MESH_TYPE), INTENT(IN), TARGET :: M
REAL(EB), INTENT(INOUT) :: A(0:M%IBP1,0:M%JBP1,0:M%KBP1)
INTEGER :: I,J,K,IBP1,JBP1,IBAR,JBAR,KBAR,KBP1

IBP1 = M%IBP1; JBP1 = M%JBP1; KBP1 = M%KBP1
IBAR = M%IBAR; JBAR = M%JBAR; KBAR = M%KBAR

! x edges

J=0; K=0
DO I=1,IBAR
   A(I,J,K) = ( A(I,J+1,K) + A(I,J,K+1) ) - A(I,J+1,K+1)
ENDDO

J=0; K=KBP1
DO I=1,IBAR
   A(I,J,K) = ( A(I,J+1,K) + A(I,J,K-1) ) - A(I,J+1,K-1)
ENDDO

J=JBP1; K=0
DO I=1,IBAR
   A(I,J,K) = ( A(I,J-1,K) + A(I,J,K+1) ) - A(I,J-1,K+1)
ENDDO

J=JBP1; K=KBP1
DO I=1,IBAR
   A(I,J,K) = ( A(I,J-1,K) + A(I,J,K-1) ) - A(I,J-1,K-1)
ENDDO

! y edges

I=0; K=0
DO J=1,JBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J,K+1) ) - A(I+1,J,K+1)
ENDDO

I=0; K=KBP1
DO J=1,JBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J,K-1) ) - A(I+1,J,K-1)
ENDDO

I=IBP1; K=0
DO J=1,JBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J,K+1) ) - A(I-1,J,K+1)
ENDDO

I=IBP1; K=KBP1
DO J=1,JBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J,K-1) ) - A(I-1,J,K-1)
ENDDO

! z edges

I=0; J=0
DO K=1,KBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J+1,K) ) - A(I+1,J+1,K)
ENDDO

I=0; J=JBP1
DO K=1,KBAR
   A(I,J,K) = ( A(I+1,J,K) + A(I,J-1,K) ) - A(I+1,J-1,K)
ENDDO

I=IBP1; J=0
DO K=1,KBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J+1,K) ) - A(I-1,J+1,K)
ENDDO

I=IBP1; J=JBP1
DO K=1,KBAR
   A(I,J,K) = ( A(I-1,J,K) + A(I,J-1,K) ) - A(I-1,J-1,K)
ENDDO

! Corners

A(0,0,0) = 2._EB*A(1,1,1) - A(2,2,2)
A(IBP1,0,0) = 2._EB*A(IBP1-1,1,1) - A(IBP1-2,2,2)
A(0,JBP1,0) = 2._EB*A(1,JBP1-1,1) - A(2,JBP1-2,2)
A(0,0,KBP1) = 2._EB*A(1,1,KBP1-1) - A(2,2,KBP1-2)
A(IBP1,JBP1,0) = 2._EB*A(IBP1-1,JBP1-1,1) - A(IBP1-2,JBP1-2,2)
A(IBP1,0,KBP1) = 2._EB*A(IBP1-1,1,KBP1-1) - A(IBP1-2,2,KBP1-2)
A(0,JBP1,KBP1) = 2._EB*A(1,JBP1-1,KBP1-1) - A(2,JBP1-2,KBP1-2)
A(IBP1,JBP1,KBP1) = 2._EB*A(IBP1-1,JBP1-1,KBP1-1) - A(IBP1-2,JBP1-2,KBP1-2)

END SUBROUTINE FILL_EDGES_KERNEL


!> \brief Apply test filter to 3D array.
!> \param M Mesh data structure
!> \param HAT Filtered output array
!> \param ORIG Original input array

SUBROUTINE TEST_FILTER_KERNEL(M,HAT,ORIG)

TYPE(MESH_TYPE), INTENT(IN), TARGET :: M
REAL(EB), INTENT(IN) :: ORIG(0:M%IBP1,0:M%JBP1,0:M%KBP1)
REAL(EB), INTENT(OUT) :: HAT(0:M%IBP1,0:M%JBP1,0:M%KBP1)
INTEGER :: I, J, K, L, MM, N, IBP1, JBP1, KBP1
REAL(EB), PARAMETER :: K1DM(3)=(/1.0_EB,1.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DM(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DM(I)*K1DM(J)*K1DM(K)/27._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))
REAL(EB), PARAMETER :: K1DT(3)=(/1.0_EB,2.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DT(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DT(I)*K1DT(J)*K1DT(K)/64._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))
REAL(EB), PARAMETER :: K1DS(3)=(/1.0_EB,4.0_EB,1.0_EB/)
REAL(EB), PARAMETER :: K3DS(-1:1,-1:1,-1:1)=RESHAPE((/(((K1DS(I)*K1DS(J)*K1DS(K)/216._EB,I=1,3),J=1,3),K=1,3)/),(/3,3,3/))

IBP1 = M%IBP1; JBP1 = M%JBP1; KBP1 = M%KBP1

! Traverse bulk of mesh

QUADRATURE_SELECT: SELECT CASE(TEST_FILTER_QUADRATURE)

   CASE(TRAPEZOID_QUADRATURE) ! default

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO MM = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+MM,K+N) * K3DT(L,MM,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   CASE(SIMPSON_QUADRATURE)

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO MM = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+MM,K+N) * K3DS(L,MM,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

   CASE(MIDPOINT_QUADRATURE)

      !$OMP PARALLEL
      !$OMP DO SCHEDULE(static)
      DO K = 1,KBP1-1
         DO J = 1,JBP1-1
            DO I = 1,IBP1-1

               ! Apply 3x3x3 Kernel; this is faster than elementwise array multiplication.
               HAT(I,J,K) = 0._EB
               DO N = -1,1
                  DO MM = -1,1
                     DO L = -1,1
                        HAT(I,J,K) = HAT(I,J,K) + ORIG(I+L,J+MM,K+N) * K3DM(L,MM,N)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

END SELECT QUADRATURE_SELECT

! Traverse shell of mesh rather crudely.
! Edges and corners are calculated several times.

!$OMP PARALLEL
!$OMP DO SCHEDULE(static)
DO K = 0,KBP1
   DO J = 0,JBP1
      HAT(0,J,K) = 2._EB * HAT(0+1,J,K) - HAT(0+2,J,K)
      HAT(IBP1,J,K) = 2._EB * HAT(IBP1-1,J,K) - HAT(IBP1-2,J,K)
   END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(static)
DO K = 0,KBP1
   DO I = 0,IBP1
      HAT(I,0,K) = 2._EB * HAT(I,0+1,K) - HAT(I,0+2,K)
      HAT(I,JBP1,K) = 2._EB * HAT(I,JBP1-1,K) - HAT(I,JBP1-2,K)
   END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(static)
DO J = 0,JBP1
   DO I = 0,IBP1
      HAT(I,J,0) = 2._EB * HAT(I,J,0+1) - HAT(I,J,0+2)
      HAT(I,J,KBP1) = 2._EB * HAT(I,J,KBP1-1) - HAT(I,J,KBP1-2)
   END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE TEST_FILTER_KERNEL


!> \brief SGS tensor diffusivity model (nonlinear backscatter and anisotropy).
!> \param M Mesh data structure
!> \param OPT_N Optional species index (if present, compute scalar flux; otherwise thermal flux)

SUBROUTINE TENSOR_DIFFUSIVITY_MODEL_KERNEL(M,OPT_N)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN), OPTIONAL :: OPT_N
INTEGER :: I,J,K,N
REAL(EB) :: DZDX,DZDY,DZDZ,DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,DTDX,DTDY,DTDZ,RHOBAR
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,UU,VV,WW,KDTDX,KDTDY,KDTDZ
REAL(EB), PARAMETER :: C_NL=0.083_EB ! C_NL=1/12, See Pope Exercise 13.28

SCALAR_FLUX_IF: IF (PRESENT(OPT_N)) THEN

   N = OPT_N

   ! SGS scalar flux
   ! CAUTION: The flux arrays must point to the same work arrays used in DIVERGENCE_PART_1
   ! Note: Do not reinitialize!  RHO_D_DZDX, etc., already store molecular diffusive flux
   RHO_D_DZDX=>M%SWORK1
   RHO_D_DZDY=>M%SWORK2
   RHO_D_DZDZ=>M%SWORK3

   IF (PREDICTOR) THEN
      UU=>M%U
      VV=>M%V
      WW=>M%W
      RHOP=>M%RHOS
      ZZP=>M%ZZS
   ELSE
      UU=>M%US
      VV=>M%VS
      WW=>M%WS
      RHOP=>M%RHO
      ZZP=>M%ZZ
   ENDIF

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=0,M%IBAR
            RHOBAR = 0.5_EB*(RHOP(I,J,K)+RHOP(I+1,J,K))

            DUDX = (UU(I+1,J,K)-UU(I-1,J,K))/(M%DX(I)+M%DX(I+1))
            DUDY = (UU(I,J+1,K)-UU(I,J-1,K))/(M%DYN(J-1)+M%DYN(J))
            DUDZ = (UU(I,J,K+1)-UU(I,J,K-1))/(M%DZN(K-1)+M%DZN(K))

            DZDX = M%RDXN(I)*(ZZP(I+1,J,K,N)-ZZP(I,J,K,N))
            DZDY = 0.25_EB*M%RDY(J)*( ZZP(I,J+1,K,N) + ZZP(I+1,J+1,K,N) - ZZP(I,J-1,K,N) - ZZP(I+1,J-1,K,N) )
            DZDZ = 0.25_EB*M%RDZ(K)*( ZZP(I,J,K+1,N) + ZZP(I+1,J,K+1,N) - ZZP(I,J,K-1,N) - ZZP(I+1,J,K-1,N) )

            RHO_D_DZDX(I,J,K,N) = RHO_D_DZDX(I,J,K,N) &
                                + RHOBAR*C_NL*( M%DXN(I)**2*DZDX*DUDX + M%DY(J)**2*DZDY*DUDY + M%DZ(K)**2*DZDZ*DUDZ )
         ENDDO
      ENDDO
   ENDDO

   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            RHOBAR = 0.5_EB*(RHOP(I,J,K)+RHOP(I,J+1,K))

            DVDX = (VV(I+1,J,K)-VV(I-1,J,K))/(M%DXN(I-1)+M%DXN(I))
            DVDY = (VV(I,J+1,K)-VV(I,J-1,K))/(M%DY(J)+M%DY(J+1))
            DVDZ = (VV(I,J,K+1)-VV(I,J,K-1))/(M%DZN(K-1)+M%DZN(K))

            DZDX = 0.25_EB*M%RDX(I)*( ZZP(I+1,J,K,N) + ZZP(I+1,J+1,K,N) - ZZP(I-1,J,K,N) - ZZP(I-1,J+1,K,N) )
            DZDY = M%RDYN(J)*(ZZP(I,J+1,K,N)-ZZP(I,J,K,N))
            DZDZ = 0.25_EB*M%RDZ(K)*( ZZP(I,J,K+1,N) + ZZP(I,J+1,K+1,N) - ZZP(I,J,K-1,N) - ZZP(I,J+1,K-1,N) )

            RHO_D_DZDY(I,J,K,N) = RHO_D_DZDY(I,J,K,N) &
                                + RHOBAR*C_NL*( M%DX(I)**2*DZDX*DVDX + M%DY(J)**2*DZDY*DVDY + M%DZ(K)**2*DZDZ*DVDZ )
         ENDDO
      ENDDO
   ENDDO

   DO K=0,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            RHOBAR = 0.5_EB*(RHOP(I,J,K)+RHOP(I,J,K+1))

            DWDX = (WW(I+1,J,K)-WW(I-1,J,K))/(M%DXN(I-1)+M%DXN(I))
            DWDY = (WW(I,J+1,K)-WW(I,J-1,K))/(M%DYN(J-1)+M%DYN(J))
            DWDZ = (WW(I,J,K+1)-WW(I,J,K-1))/(M%DZ(K)+M%DZ(K+1))

            DZDX = 0.25_EB*M%RDX(I)*( ZZP(I+1,J,K,N) + ZZP(I+1,J,K+1,N) - ZZP(I-1,J,K,N) - ZZP(I-1,J,K+1,N) )
            DZDY = 0.25_EB*M%RDY(J)*( ZZP(I,J+1,K,N) + ZZP(I,J+1,K+1,N) - ZZP(I,J-1,K,N) - ZZP(I,J-1,K+1,N) )
            DZDZ = M%RDZN(K)*(ZZP(I,J,K+1,N)-ZZP(I,J,K,N))

            RHO_D_DZDZ(I,J,K,N) = RHO_D_DZDZ(I,J,K,N) &
                                + RHOBAR*C_NL*( M%DX(I)**2*DZDX*DWDX + M%DY(J)**2*DZDY*DWDY + M%DZ(K)**2*DZDZ*DWDZ )
         ENDDO
      ENDDO
   ENDDO

ELSE SCALAR_FLUX_IF

   ! SGS thermal energy flux
   ! CAUTION: The flux arrays must point to the same work arrays used in DIVERGENCE_PART_1
   ! Note: Do not reinitialize!  KDTDX, etc., already store molecular diffusive flux
   KDTDX=>M%WORK1
   KDTDY=>M%WORK2
   KDTDZ=>M%WORK3

   IF (PREDICTOR) THEN
      UU=>M%U
      VV=>M%V
      WW=>M%W
   ELSE
      UU=>M%US
      VV=>M%VS
      WW=>M%WS
   ENDIF

   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=0,M%IBAR

            DUDX = (UU(I+1,J,K)-UU(I-1,J,K))/(M%DX(I)+M%DX(I+1))
            DUDY = (UU(I,J+1,K)-UU(I,J-1,K))/(M%DYN(J-1)+M%DYN(J))
            DUDZ = (UU(I,J,K+1)-UU(I,J,K-1))/(M%DZN(K-1)+M%DZN(K))

            DTDX = M%RDXN(I)*(M%TMP(I+1,J,K)-M%TMP(I,J,K))
            DTDY = 0.25_EB*M%RDY(J)*( M%TMP(I,J+1,K) + M%TMP(I+1,J+1,K) - M%TMP(I,J-1,K) - M%TMP(I+1,J-1,K) )
            DTDZ = 0.25_EB*M%RDZ(K)*( M%TMP(I,J,K+1) + M%TMP(I+1,J,K+1) - M%TMP(I,J,K-1) - M%TMP(I+1,J,K-1) )

            KDTDX(I,J,K) = KDTDX(I,J,K) + C_NL*(M%DX(I)**2*DUDX*DTDX + M%DY(J)**2*DUDY*DTDY + M%DZ(K)**2*DUDZ*DTDZ)

         ENDDO
      ENDDO
   ENDDO

   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR

            DVDX = (VV(I+1,J,K)-VV(I-1,J,K))/(M%DXN(I-1)+M%DXN(I))
            DVDY = (VV(I,J+1,K)-VV(I,J-1,K))/(M%DY(J)+M%DY(J+1))
            DVDZ = (VV(I,J,K+1)-VV(I,J,K-1))/(M%DZN(K-1)+M%DZN(K))

            DTDX = 0.25_EB*M%RDX(I)*( M%TMP(I+1,J,K) + M%TMP(I+1,J+1,K) - M%TMP(I-1,J,K) - M%TMP(I-1,J+1,K) )
            DTDY = M%RDYN(J)*(M%TMP(I,J+1,K)-M%TMP(I,J,K))
            DTDZ = 0.25_EB*M%RDZ(K)*( M%TMP(I,J,K+1) + M%TMP(I,J+1,K+1) - M%TMP(I,J,K-1) - M%TMP(I,J+1,K-1) )

            KDTDY(I,J,K) = KDTDY(I,J,K) + C_NL*(M%DX(I)**2*DVDX*DTDX + M%DY(J)**2*DVDY*DTDY + M%DZ(K)**2*DVDZ*DTDZ)

         ENDDO
      ENDDO
   ENDDO

   DO K=0,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR

            DWDX = (WW(I+1,J,K)-WW(I-1,J,K))/(M%DXN(I-1)+M%DXN(I))
            DWDY = (WW(I,J+1,K)-WW(I,J-1,K))/(M%DYN(J-1)+M%DYN(J))
            DWDZ = (WW(I,J,K+1)-WW(I,J,K-1))/(M%DZ(K)+M%DZ(K+1))

            DTDX = 0.25_EB*M%RDX(I)*( M%TMP(I+1,J,K) + M%TMP(I+1,J,K+1) - M%TMP(I-1,J,K) - M%TMP(I-1,J,K+1) )
            DTDY = 0.25_EB*M%RDY(J)*( M%TMP(I,J+1,K) + M%TMP(I,J+1,K+1) - M%TMP(I,J-1,K) - M%TMP(I,J-1,K+1) )
            DTDZ = M%RDZN(K)*(M%TMP(I,J,K+1)-M%TMP(I,J,K))

            KDTDZ(I,J,K) = KDTDZ(I,J,K) + C_NL*(M%DX(I)**2*DWDX*DTDX + M%DY(J)**2*DWDY*DTDY + M%DZ(K)**2*DWDZ*DTDZ)

         ENDDO
      ENDDO
   ENDDO

ENDIF SCALAR_FLUX_IF

END SUBROUTINE TENSOR_DIFFUSIVITY_MODEL_KERNEL


!> \brief Dynamic Smagorinsky model (Germano procedure).
!> \param M Mesh data structure

SUBROUTINE VARDEN_DYNSMAG_KERNEL(M)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M

REAL(EB) :: TEMP_TERM,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,ONTHDIV
INTEGER :: I,J,K

REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,UP,VP,WP,RHOP,RHOPHAT,RHOPTMP
REAL(EB), POINTER, DIMENSION(:,:,:) :: S11,S22,S33,S12,S13,S23,SS
REAL(EB), POINTER, DIMENSION(:,:,:) :: SHAT11,SHAT22,SHAT33,SHAT12,SHAT13,SHAT23,SSHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: BETA11,BETA22,BETA33,BETA12,BETA13,BETA23
REAL(EB), POINTER, DIMENSION(:,:,:) :: BETAHAT11,BETAHAT22,BETAHAT33,BETAHAT12,BETAHAT13,BETAHAT23
REAL(EB), POINTER, DIMENSION(:,:,:) :: M11,M22,M33,M12,M13,M23,MM,MMHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: L11,L22,L33,L12,L13,L23,ML,MLHAT

REAL(EB), PARAMETER :: ALPHA = 6.0_EB ! See Lund, 1997 CTR briefs.

! References:
!
! M. Germano, U. Piomelli, P. Moin, and W. Cabot.  A dynamic subgrid-scale eddy viscosity model.
! Phys. Fluids A, 3(7):1760-1765, 1991.
!
! M. Pino Martin, U. Piomelli, and G. Candler. Subgrid-scale models for compressible large-eddy
! simulation. Theoret. Comput. Fluid Dynamics, 13:361-376, 2000.
!
! P. Moin, K. Squires, W. Cabot, and S. Lee.  A dynamic subgrid-scale model for compressible
! turbulence and scalar transport. Phys. Fluids A, 3(11):2746-2757, 1991.
!
! T. S. Lund. On the use of discrete filters for large eddy simulation.  Center for Turbulence
! Research Annual Research Briefs, 1997.
!
! R. McDermott. Variable density formulation of the dynamic Smagorinsky model.
! http://randy.mcdermott.googlepages.com/dynsmag_comp.pdf

! *****************************************************************************
! CAUTION WHEN MODIFYING: The order in which the tensor components are computed
! is important because we overwrite pointers several times to conserve memory.
! *****************************************************************************

IF (PREDICTOR) THEN
   UU=>M%U
   VV=>M%V
   WW=>M%W
   RHOP=>M%RHO
ELSE
   UU=>M%US
   VV=>M%VS
   WW=>M%WS
   RHOP=>M%RHOS
ENDIF

UP => M%TURB_WORK1
VP => M%TURB_WORK2
WP => M%TURB_WORK3

S11 => M%WORK1
S22 => M%WORK2
S33 => M%WORK3
S12 => M%WORK4
S13 => M%WORK5
S23 => M%WORK6
SS  => M%WORK7

DO K = 1,M%KBAR
   DO J = 1,M%JBAR
      DO I = 1,M%IBAR

         UP(I,J,K) = 0.5_EB*(UU(I,J,K) + UU(I-1,J,K))
         VP(I,J,K) = 0.5_EB*(VV(I,J,K) + VV(I,J-1,K))
         WP(I,J,K) = 0.5_EB*(WW(I,J,K) + WW(I,J,K-1))

         S11(I,J,K) = M%RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
         S22(I,J,K) = M%RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
         S33(I,J,K) = M%RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))

         ONTHDIV = ONTH*(S11(I,J,K)+S22(I,J,K)+S33(I,J,K))
         S11(I,J,K) = S11(I,J,K)-ONTHDIV
         S22(I,J,K) = S22(I,J,K)-ONTHDIV
         S33(I,J,K) = S33(I,J,K)-ONTHDIV

         DUDY = 0.25_EB*M%RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
         DUDZ = 0.25_EB*M%RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1))
         DVDX = 0.25_EB*M%RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
         DVDZ = 0.25_EB*M%RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
         DWDX = 0.25_EB*M%RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
         DWDY = 0.25_EB*M%RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
         S12(I,J,K) = 0.5_EB*(DUDY+DVDX)
         S13(I,J,K) = 0.5_EB*(DUDZ+DWDX)
         S23(I,J,K) = 0.5_EB*(DVDZ+DWDY)

         ! calculate magnitude of the grid strain rate

         SS(I,J,K) = RHOP(I,J,K)*M%STRAIN_RATE(I,J,K)

      ENDDO
   ENDDO
ENDDO

! second-order extrapolation to ghost cells

CALL EX2G3D_KERNEL(M,UP,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,VP,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,WP,-1.E10_EB,1.E10_EB)

CALL EX2G3D_KERNEL(M,S11,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,S22,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,S33,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,S12,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,S13,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,S23,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,SS,0._EB,1.E10_EB)

! test filter the strain rate

SHAT11 => M%TURB_WORK4
SHAT22 => M%TURB_WORK5
SHAT33 => M%TURB_WORK6
SHAT12 => M%TURB_WORK7
SHAT13 => M%TURB_WORK8
SHAT23 => M%TURB_WORK9
SSHAT  => M%TURB_WORK10

CALL TEST_FILTER_KERNEL(M,SHAT11,S11)
CALL TEST_FILTER_KERNEL(M,SHAT22,S22)
CALL TEST_FILTER_KERNEL(M,SHAT33,S33)
CALL TEST_FILTER_KERNEL(M,SHAT12,S12)
CALL TEST_FILTER_KERNEL(M,SHAT13,S13)
CALL TEST_FILTER_KERNEL(M,SHAT23,S23)


! calculate magnitude of test filtered strain rate

DO K = 1,M%KBAR
   DO J = 1,M%JBAR
      DO I = 1,M%IBAR
         SSHAT(I,J,K) = SQRT(2._EB*(SHAT11(I,J,K)*SHAT11(I,J,K) + &
                                    SHAT22(I,J,K)*SHAT22(I,J,K) + &
                                    SHAT33(I,J,K)*SHAT33(I,J,K) + &
                             2._EB*(SHAT12(I,J,K)*SHAT12(I,J,K) + &
                                    SHAT13(I,J,K)*SHAT13(I,J,K) + &
                                    SHAT23(I,J,K)*SHAT23(I,J,K)) ) )
      ENDDO
   ENDDO
ENDDO

! calculate the grid filtered stress tensor, beta

BETA11 => M%WORK1
BETA22 => M%WORK2
BETA33 => M%WORK3
BETA12 => M%WORK4
BETA13 => M%WORK5
BETA23 => M%WORK6

BETA11 = SS*S11
BETA22 = SS*S22
BETA33 = SS*S33
BETA12 = SS*S12
BETA13 = SS*S13
BETA23 = SS*S23

! ghost values for beta_ij should be filled already

! test filter the grid filtered stress tensor

BETAHAT11 => M%WORK1
BETAHAT22 => M%WORK2
BETAHAT33 => M%WORK3
BETAHAT12 => M%WORK4
BETAHAT13 => M%WORK5
BETAHAT23 => M%WORK6

M%WORK9=BETA11; CALL TEST_FILTER_KERNEL(M,BETAHAT11,M%WORK9)
M%WORK9=BETA22; CALL TEST_FILTER_KERNEL(M,BETAHAT22,M%WORK9)
M%WORK9=BETA33; CALL TEST_FILTER_KERNEL(M,BETAHAT33,M%WORK9)
M%WORK9=BETA12; CALL TEST_FILTER_KERNEL(M,BETAHAT12,M%WORK9)
M%WORK9=BETA13; CALL TEST_FILTER_KERNEL(M,BETAHAT13,M%WORK9)
M%WORK9=BETA23; CALL TEST_FILTER_KERNEL(M,BETAHAT23,M%WORK9)

! test filter the density

RHOPHAT => M%WORK7
RHOPTMP => M%WORK1
RHOPTMP(0:M%IBP1,0:M%JBP1,0:M%KBP1) = RHOP(0:M%IBP1,0:M%JBP1,0:M%KBP1)
CALL TEST_FILTER_KERNEL(M,RHOPHAT,RHOPTMP)

! calculate the Mij tensor

M11 => M%WORK1
M22 => M%WORK2
M33 => M%WORK3
M12 => M%WORK4
M13 => M%WORK5
M23 => M%WORK6

DO K = 1,M%KBAR
   DO J = 1,M%JBAR
      DO I = 1,M%IBAR
         TEMP_TERM = ALPHA*RHOPHAT(I,J,K)*SSHAT(I,J,K)
         M11(I,J,K) = 2._EB*(BETAHAT11(I,J,K) - TEMP_TERM*SHAT11(I,J,K))
         M22(I,J,K) = 2._EB*(BETAHAT22(I,J,K) - TEMP_TERM*SHAT22(I,J,K))
         M33(I,J,K) = 2._EB*(BETAHAT33(I,J,K) - TEMP_TERM*SHAT33(I,J,K))
         M12(I,J,K) = 2._EB*(BETAHAT12(I,J,K) - TEMP_TERM*SHAT12(I,J,K))
         M13(I,J,K) = 2._EB*(BETAHAT13(I,J,K) - TEMP_TERM*SHAT13(I,J,K))
         M23(I,J,K) = 2._EB*(BETAHAT23(I,J,K) - TEMP_TERM*SHAT23(I,J,K))
      ENDDO
   ENDDO
ENDDO

! calculate the Leonard term, Lij

L11 => M%TURB_WORK4
L22 => M%TURB_WORK5
L33 => M%TURB_WORK6
L12 => M%TURB_WORK7
L13 => M%TURB_WORK8
L23 => M%TURB_WORK9

CALL CALC_VARDEN_LEONARD_TERM_KERNEL

! calculate Mij*Lij & Mij*Mij

MM    => M%TURB_WORK1
MMHAT => M%TURB_WORK1

ML    => M%TURB_WORK2
MLHAT => M%TURB_WORK2

DO K = 1,M%KBAR
   DO J = 1,M%JBAR
      DO I = 1,M%IBAR

         ML(I,J,K) = M11(I,J,K)*L11(I,J,K) + M22(I,J,K)*L22(I,J,K) + M33(I,J,K)*L33(I,J,K) + &
              2._EB*(M12(I,J,K)*L12(I,J,K) + M13(I,J,K)*L13(I,J,K) + M23(I,J,K)*L23(I,J,K))

         MM(I,J,K) = M11(I,J,K)*M11(I,J,K) + M22(I,J,K)*M22(I,J,K) + M33(I,J,K)*M33(I,J,K) + &
              2._EB*(M12(I,J,K)*M12(I,J,K) + M13(I,J,K)*M13(I,J,K) + M23(I,J,K)*M23(I,J,K))

      ENDDO
   ENDDO
ENDDO

! extrapolate to ghost

CALL EX2G3D_KERNEL(M,ML,0._EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,MM,0._EB,1.E10_EB)

! do some smoothing

M%WORK9=ML; CALL TEST_FILTER_KERNEL(M,MLHAT,M%WORK9)
M%WORK9=MM; CALL TEST_FILTER_KERNEL(M,MMHAT,M%WORK9)

DO K = 1,M%KBAR
   DO J = 1,M%JBAR
      DO I = 1,M%IBAR

         ! calculate the local Smagorinsky coefficient

         ! perform "clipping" in case MLij is negative
         IF (MLHAT(I,J,K) < 0._EB) MLHAT(I,J,K) = 0._EB

         ! calculate the effective viscosity

         ! handle the case where we divide by zero, note MMHAT is positive semi-definite
         IF (MMHAT(I,J,K)<=TWENTY_EPSILON_EB) THEN
            M%CSD2(I,J,K) = 0._EB
         ELSE
            M%CSD2(I,J,K) = MLHAT(I,J,K)/MMHAT(I,J,K) ! (Cs*Delta)**2
         ENDIF

      END DO
   END DO
END DO

CONTAINS

SUBROUTINE CALC_VARDEN_LEONARD_TERM_KERNEL

REAL(EB), POINTER, DIMENSION(:,:,:) :: LR11,LR22,LR33,LR12,LR13,LR23
REAL(EB), POINTER, DIMENSION(:,:,:) :: LRHOP,LRHOPHAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: LUP,LVP,LWP
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU,RVV,RWW,RUV,RUW,RVW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU,RV,RW
REAL(EB), POINTER, DIMENSION(:,:,:) :: RUU_HAT,RVV_HAT,RWW_HAT,RUV_HAT,RUW_HAT,RVW_HAT
REAL(EB), POINTER, DIMENSION(:,:,:) :: RU_HAT,RV_HAT,RW_HAT
REAL(EB) :: INV_RHOPHAT
INTEGER :: II,JJ,KK,IBP1L,JBP1L,KBP1L

IBP1L = M%IBP1; JBP1L = M%JBP1; KBP1L = M%KBP1

! *****************************************************************************
! CAUTION WHEN MODIFYING: The order in which the tensor components are computed
! is important because we overwrite pointers several times to conserve memory.
! *****************************************************************************

IF (PREDICTOR) THEN
   LRHOP=>M%RHO
ELSE
   LRHOP=>M%RHOS
ENDIF
LRHOPHAT => M%WORK7

! Compute rho*UiUj

LUP => M%TURB_WORK1 ! will be overwritten by RU
LVP => M%TURB_WORK2
LWP => M%TURB_WORK3

RUU => M%TURB_WORK4 ! will be overwritten by RUU_HAT
RVV => M%TURB_WORK5
RWW => M%TURB_WORK6
RUV => M%TURB_WORK7
RUW => M%TURB_WORK8
RVW => M%TURB_WORK9

RUU(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LUP(0:IBP1L,0:JBP1L,0:KBP1L)*LUP(0:IBP1L,0:JBP1L,0:KBP1L)
RVV(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LVP(0:IBP1L,0:JBP1L,0:KBP1L)*LVP(0:IBP1L,0:JBP1L,0:KBP1L)
RWW(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LWP(0:IBP1L,0:JBP1L,0:KBP1L)*LWP(0:IBP1L,0:JBP1L,0:KBP1L)
RUV(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LUP(0:IBP1L,0:JBP1L,0:KBP1L)*LVP(0:IBP1L,0:JBP1L,0:KBP1L)
RUW(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LUP(0:IBP1L,0:JBP1L,0:KBP1L)*LWP(0:IBP1L,0:JBP1L,0:KBP1L)
RVW(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LVP(0:IBP1L,0:JBP1L,0:KBP1L)*LWP(0:IBP1L,0:JBP1L,0:KBP1L)

! extrapolate to ghost cells

CALL EX2G3D_KERNEL(M,RUU,0.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RVV,0.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RWW,0.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RUV,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RUW,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RVW,-1.E10_EB,1.E10_EB)

! Test filter rho*UiUj

RUU_HAT => M%TURB_WORK4 ! will be overwritten by Lij
RVV_HAT => M%TURB_WORK5
RWW_HAT => M%TURB_WORK6
RUV_HAT => M%TURB_WORK7
RUW_HAT => M%TURB_WORK8
RVW_HAT => M%TURB_WORK9

M%WORK9=RUU; CALL TEST_FILTER_KERNEL(M,RUU_HAT,M%WORK9)
M%WORK9=RVV; CALL TEST_FILTER_KERNEL(M,RVV_HAT,M%WORK9)
M%WORK9=RWW; CALL TEST_FILTER_KERNEL(M,RWW_HAT,M%WORK9)
M%WORK9=RUV; CALL TEST_FILTER_KERNEL(M,RUV_HAT,M%WORK9)
M%WORK9=RUW; CALL TEST_FILTER_KERNEL(M,RUW_HAT,M%WORK9)
M%WORK9=RVW; CALL TEST_FILTER_KERNEL(M,RVW_HAT,M%WORK9)

! Compute rho*Ui

RU => M%TURB_WORK1 ! will be overwritten by RU_HAT
RV => M%TURB_WORK2
RW => M%TURB_WORK3

RU(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LUP(0:IBP1L,0:JBP1L,0:KBP1L)
RV(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LVP(0:IBP1L,0:JBP1L,0:KBP1L)
RW(0:IBP1L,0:JBP1L,0:KBP1L) = LRHOP(0:IBP1L,0:JBP1L,0:KBP1L)*LWP(0:IBP1L,0:JBP1L,0:KBP1L)

! extrapolate to ghost cells

CALL EX2G3D_KERNEL(M,RU,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RV,-1.E10_EB,1.E10_EB)
CALL EX2G3D_KERNEL(M,RW,-1.E10_EB,1.E10_EB)

! Test filter rho*Ui

RU_HAT => M%TURB_WORK1
RV_HAT => M%TURB_WORK2
RW_HAT => M%TURB_WORK3

M%WORK9=RU; CALL TEST_FILTER_KERNEL(M,RU_HAT,M%WORK9)
M%WORK9=RV; CALL TEST_FILTER_KERNEL(M,RV_HAT,M%WORK9)
M%WORK9=RW; CALL TEST_FILTER_KERNEL(M,RW_HAT,M%WORK9)

! Compute variable density Leonard stress

LR11 => M%TURB_WORK4
LR22 => M%TURB_WORK5
LR33 => M%TURB_WORK6
LR12 => M%TURB_WORK7
LR13 => M%TURB_WORK8
LR23 => M%TURB_WORK9

DO KK = 1,M%KBAR
   DO JJ = 1,M%JBAR
      DO II = 1,M%IBAR
         INV_RHOPHAT = 1._EB/LRHOPHAT(II,JJ,KK)
         LR11(II,JJ,KK) = RUU_HAT(II,JJ,KK) - RU_HAT(II,JJ,KK)*RU_HAT(II,JJ,KK)*INV_RHOPHAT
         LR22(II,JJ,KK) = RVV_HAT(II,JJ,KK) - RV_HAT(II,JJ,KK)*RV_HAT(II,JJ,KK)*INV_RHOPHAT
         LR33(II,JJ,KK) = RWW_HAT(II,JJ,KK) - RW_HAT(II,JJ,KK)*RW_HAT(II,JJ,KK)*INV_RHOPHAT
         LR12(II,JJ,KK) = RUV_HAT(II,JJ,KK) - RU_HAT(II,JJ,KK)*RV_HAT(II,JJ,KK)*INV_RHOPHAT
         LR13(II,JJ,KK) = RUW_HAT(II,JJ,KK) - RU_HAT(II,JJ,KK)*RW_HAT(II,JJ,KK)*INV_RHOPHAT
         LR23(II,JJ,KK) = RVW_HAT(II,JJ,KK) - RV_HAT(II,JJ,KK)*RW_HAT(II,JJ,KK)*INV_RHOPHAT
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE CALC_VARDEN_LEONARD_TERM_KERNEL

END SUBROUTINE VARDEN_DYNSMAG_KERNEL


END MODULE TURB_KERNELS
