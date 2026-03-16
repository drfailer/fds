!> \brief Thread-safe computation kernels extracted from PART module.
!> These routines take TYPE(MESH_TYPE) as an argument instead of relying on
!> MESH_POINTERS / POINT_TO_MESH, decoupling computation from global state.

MODULE PART_KERNELS

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES, ONLY: MESH_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC PARTICLE_MOMENTUM_TRANSFER_KERNEL,PARTICLE_MOMENTUM_BLOCK_KERNEL

CONTAINS


!> \brief Add PARTICLE momentum as a force term in momentum equation (thread-safe kernel).
!> \param M Mesh data structure
!> \param DT Current time step (s)
!> \param NM Current mesh number

SUBROUTINE PARTICLE_MOMENTUM_TRANSFER_KERNEL(M,DT,NM)

USE CC_VELOCITY_KERNELS, ONLY: CUTFACE_VELOCITIES
TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: DT
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB) :: RDT,UODT,VODT,WODT
INTEGER :: I,J,K

IF (M%NLP==0) RETURN

RDT = 1._EB/DT

IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
ENDIF

IF (CC_IBM) CALL CUTFACE_VELOCITIES(M,UU,VV,WW, &
   CUTFACES=.TRUE.)

! Add summed particle accelerations to the momentum equation. Limit the value to plus/minus abs(u)/dt to prevent a sudden
! change in gas direction.

DO K=0,M%KBAR
   DO J=0,M%JBAR
      DO I=0,M%IBAR
         UODT = ABS(UU(I,J,K)*RDT)
         VODT = ABS(VV(I,J,K)*RDT)
         WODT = ABS(WW(I,J,K)*RDT)
         M%FVX(I,J,K) = M%FVX(I,J,K) + MIN(UODT,MAX(-UODT,M%FVX_D(I,J,K)))
         M%FVY(I,J,K) = M%FVY(I,J,K) + MIN(VODT,MAX(-VODT,M%FVY_D(I,J,K)))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + MIN(WODT,MAX(-WODT,M%FVZ_D(I,J,K)))
      ENDDO
   ENDDO
ENDDO

IF (CC_IBM) CALL CUTFACE_VELOCITIES(M,UU,VV,WW, &
   CUTFACES=.FALSE.)

END SUBROUTINE PARTICLE_MOMENTUM_TRANSFER_KERNEL


!> \brief Block-decomposed particle momentum: updates FVX/FVY/FVZ for K-range [K1, K2].
!> \details First block (K1==1) extends down to K=0 to cover the full 0:KBAR range.
!> Does NOT include CC_IBM CUTFACE_VELOCITIES calls (mesh-level path handles CC_IBM).
!> \param M Mesh data structure
!> \param DT Time step (s)
!> \param K1 Start of K cell range (1-based inclusive, from decompose state)
!> \param K2 End of K cell range (1-based inclusive)

RECURSIVE SUBROUTINE PARTICLE_MOMENTUM_BLOCK_KERNEL(M,DT,K1,K2)

TYPE(MESH_TYPE), INTENT(INOUT), TARGET :: M
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN) :: K1,K2
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
REAL(EB) :: RDT,UODT,VODT,WODT
INTEGER :: I,J,K,K_START

IF (M%NLP==0) RETURN

RDT = 1._EB/DT

IF (PREDICTOR) THEN
   UU => M%U
   VV => M%V
   WW => M%W
ELSE
   UU => M%US
   VV => M%VS
   WW => M%WS
ENDIF

! First block extends down to K=0 to cover ghost face layer
K_START = K1
IF (K1==1) K_START = 0

DO K=K_START,K2
   DO J=0,M%JBAR
      DO I=0,M%IBAR
         UODT = ABS(UU(I,J,K)*RDT)
         VODT = ABS(VV(I,J,K)*RDT)
         WODT = ABS(WW(I,J,K)*RDT)
         M%FVX(I,J,K) = M%FVX(I,J,K) + MIN(UODT,MAX(-UODT,M%FVX_D(I,J,K)))
         M%FVY(I,J,K) = M%FVY(I,J,K) + MIN(VODT,MAX(-VODT,M%FVY_D(I,J,K)))
         M%FVZ(I,J,K) = M%FVZ(I,J,K) + MIN(WODT,MAX(-WODT,M%FVZ_D(I,J,K)))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE PARTICLE_MOMENTUM_BLOCK_KERNEL

END MODULE PART_KERNELS
