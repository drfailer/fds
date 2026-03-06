!  +++++++++++++++++++++++ CC_VERIFICATION ++++++++++++++++++++++++++

! Rotated cube verification/test routines for the
! cut-cell / immersed-boundary method.

MODULE CC_VERIFICATION

USE CC_SCALARS_DATA
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CC_ROTATED_CUBE_RHS_ZZ, &
          ROTATED_CUBE_ANN_SOLN, ROTATED_CUBE_VELOCITY_FLUX, &
          ROTATED_CUBE_RHS_ZZ

CONTAINS

! ---------------------------- ROTATED_CUBE_VELOCITY_FLUX --------------------------

SUBROUTINE ROTATED_CUBE_VELOCITY_FLUX(NM,TLEVEL)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN):: TLEVEL

! Local Variables:
REAL(EB) :: XGLOB(2,1), XLOC(2,1), FGLOB(2,1), FLOC(2,1), X_I, Y_J, NU
INTEGER  :: I,J,K
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP
REAL(EB) :: SIN_T, COS_T, COS_KX, SIN_KX, COS_KY, SIN_KY

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   RHOP => RHO
ELSE
   RHOP => RHOS
ENDIF

ROTANG = 0._EB
IF(PERIODIC_TEST==21) THEN
   ROTANG = 0._EB ! No rotation.
ELSEIF(PERIODIC_TEST==22) THEN
   ROTANG = ATAN(1._EB/2._EB) ! ~27 Degrees.
ELSEIF(PERIODIC_TEST==23) THEN
   ROTANG = ATAN(1._EB)       ! 45 degrees.
ELSE
   RETURN
ENDIF
ROTMAT(1,1) = COS(ROTANG); ROTMAT(1,2) = -SIN(ROTANG);
ROTMAT(2,1) = SIN(ROTANG); ROTMAT(2,2) =  COS(ROTANG);
TROTMAT = TRANSPOSE(ROTMAT)

COS_T = COS(TLEVEL)
SIN_T = SIN(TLEVEL)

! X Force:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         IF (CC_IBM) THEN
            IF (FCVAR(I,J,K,CC_FGSC,IAXIS) == CC_SOLID) CYCLE    ! Stress method.
         ENDIF
         ! Kinematic Viscosity:
         NU = 0.5_EB*(MU(I,J,K)/RHOP(I,J,K) + MU(I+1,J,K)/RHOP(I+1,J,K))

         ! Global position:
         XGLOB(1:2,1) = (/ X(I), ZC(K) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         X_I = XLOC(1,1); Y_J = XLOC(2,1)

         COS_KX = COS(NWAVE*X_I)
         SIN_KX = SIN(NWAVE*X_I)

         COS_KY = COS(NWAVE*Y_J)
         SIN_KY = SIN(NWAVE*Y_J)

         FLOC(1,1)=2._EB*COS_KY*SIN_KX**2._EB*SIN_KY*COS_T - &
         4._EB*NWAVE*SIN_KY*SIN_T*(COS_KX*SIN_KX**3._EB*SIN_KY*SIN_T + &
         NWAVE*NU*COS_KX**2._EB*COS_KY - 3._EB*NWAVE*NU*COS_KY*SIN_KX**2._EB) + &
         NWAVE*COS_KX*SIN_KX*SIN_T*(2._EB*COS_KY**2._EB + 1._EB);


         FLOC(2,1)= NWAVE*COS_KY*SIN_KY*SIN_T*(2._EB*COS_KX**2._EB + 1._EB) - &
         2._EB*COS_KX*SIN_KX*SIN_KY**2._EB*COS_T - &
         4._EB*NWAVE*SIN_KX*SIN_T*(COS_KY*SIN_KX*SIN_KY**3._EB*SIN_T - &
         NWAVE*NU*COS_KX*COS_KY**2._EB + 3._EB*NWAVE*NU*COS_KX*SIN_KY**2._EB);

         FGLOB        = MATMUL(ROTMAT, FLOC )

         FVX(I,J,K)   = FVX(I,J,K) + 0.5_EB*(RHOP(I,J,K)+RHOP(I+1,J,K))*FGLOB(IAXIS,1)

      ENDDO
   ENDDO
ENDDO


! Z Force:
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CC_IBM) THEN
            IF (FCVAR(I,J,K,CC_FGSC,KAXIS) == CC_SOLID) CYCLE    ! Stress method.
         ENDIF
         ! Kinematic Viscosity:
         NU = 0.5_EB*(MU(I,J,K)/RHOP(I,J,K) + MU(I,J,K+1)/RHOP(I,J,K+1))

         ! Global position:
         XGLOB(1:2,1) = (/ XC(I), Z(K) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         X_I = XLOC(1,1); Y_J = XLOC(2,1)

         COS_KX = COS(NWAVE*X_I)
         SIN_KX = SIN(NWAVE*X_I)

         COS_KY = COS(NWAVE*Y_J)
         SIN_KY = SIN(NWAVE*Y_J)

         FLOC(1,1)=2._EB*COS_KY*SIN_KX**2._EB*SIN_KY*COS_T - &
         4._EB*NWAVE*SIN_KY*SIN_T*(COS_KX*SIN_KX**3._EB*SIN_KY*SIN_T + &
         NWAVE*NU*COS_KX**2._EB*COS_KY - 3._EB*NWAVE*NU*COS_KY*SIN_KX**2._EB) + &
         NWAVE*COS_KX*SIN_KX*SIN_T*(2._EB*COS_KY**2._EB + 1._EB);


         FLOC(2,1)= NWAVE*COS_KY*SIN_KY*SIN_T*(2._EB*COS_KX**2._EB + 1._EB) - &
         2._EB*COS_KX*SIN_KX*SIN_KY**2._EB*COS_T - &
         4._EB*NWAVE*SIN_KX*SIN_T*(COS_KY*SIN_KX*SIN_KY**3._EB*SIN_T - &
         NWAVE*NU*COS_KX*COS_KY**2._EB + 3._EB*NWAVE*NU*COS_KX*SIN_KY**2._EB);

         FGLOB        = MATMUL(ROTMAT, FLOC )

         FVZ(I,J,K)   = FVZ(I,J,K) + 0.5_EB*(RHOP(I,J,K)+RHOP(I,J,K+1))*FGLOB(JAXIS,1)

      ENDDO
   ENDDO
ENDDO


RETURN
END SUBROUTINE ROTATED_CUBE_VELOCITY_FLUX

! ------------------------------ ROTATED_CUBE_ANN_SOLN ----------------------------

SUBROUTINE ROTATED_CUBE_ANN_SOLN(NM,TLEVEL)

INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN):: TLEVEL

! Local Variables:
REAL(EB) :: XGLOB(2,1), XLOC(2,1), UGLOB(2,1), ULOC(2,1), VEL_CF
INTEGER :: I,J,K,NFACE,X1AXIS,ICF,ICF2,ICC,JCC,NCELL

ROTANG = 0._EB
IF(PERIODIC_TEST==21) THEN
   ROTANG = 0._EB ! No rotation.
ELSEIF(PERIODIC_TEST==22) THEN
   ROTANG = ATAN(1._EB/2._EB) ! ~27 Degrees.
ELSEIF(PERIODIC_TEST==23) THEN
   ROTANG = ATAN(1._EB)       ! 45 degrees.
ELSE
   RETURN
ENDIF
ROTMAT(1,1) = COS(ROTANG); ROTMAT(1,2) = -SIN(ROTANG);
ROTMAT(2,1) = SIN(ROTANG); ROTMAT(2,2) =  COS(ROTANG);
TROTMAT = TRANSPOSE(ROTMAT)

CALL POINT_TO_MESH(NM)


! X Velocities:
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBAR
         ! Global position:
         XGLOB(1:2,1) = (/ X(I), ZC(K) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         ! Velocity field in local axes:
         ULOC(IAXIS,1)= -SIN(NWAVE*XLOC(IAXIS,1))**2._EB * SIN(2._EB*NWAVE*XLOC(JAXIS,1))
         ULOC(JAXIS,1)=  SIN(2._EB*NWAVE*XLOC(IAXIS,1))  * SIN(NWAVE*XLOC(JAXIS,1))**2._EB

         ! Velocity field in global axes:
         UGLOB        = MATMUL(ROTMAT,ULOC)
         U(I,J,K)     = SIN(TLEVEL) * UGLOB(IAXIS,1)

      ENDDO
   ENDDO
ENDDO

! Y Velocities:
V(:,:,:)=0._EB

! Z Velocities:
DO K=0,KBAR
   DO J=0,JBP1
      DO I=0,IBP1
         ! Global position:
         XGLOB(1:2,1) = (/ XC(I), Z(K) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         ! Velocity field in local axes:
         ULOC(IAXIS,1)= -SIN(NWAVE*XLOC(IAXIS,1))**2._EB * SIN(2._EB*NWAVE*XLOC(JAXIS,1))
         ULOC(JAXIS,1)=  SIN(2._EB*NWAVE*XLOC(IAXIS,1))  * SIN(NWAVE*XLOC(JAXIS,1))**2._EB

         ! Velocity field in global axes:
         UGLOB        = MATMUL(ROTMAT,ULOC)
         W(I,J,K)     = SIN(TLEVEL) * UGLOB(JAXIS,1)

      ENDDO
   ENDDO
ENDDO

! Now GASPHASE cut-faces:
IF (CC_IBM) THEN
   CUTFACE_LOOP : DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
      NFACE  = CUT_FACE(ICF)%NFACE
      IF (CUT_FACE(ICF)%STATUS == CC_GASPHASE) THEN
         I      = CUT_FACE(ICF)%IJK(IAXIS)
         J      = CUT_FACE(ICF)%IJK(JAXIS)
         K      = CUT_FACE(ICF)%IJK(KAXIS)
         X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
         DO ICF2=1,NFACE
            ! Global position:
            XGLOB(1:2,1) = (/ CUT_FACE(ICF)%XYZCEN(IAXIS,ICF2), CUT_FACE(ICF)%XYZCEN(KAXIS,ICF2) /)

            ! Local position:
            XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
            XLOC         = MATMUL(TROTMAT, XGLOB )
            XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

            ! Velocity field in local axes:
            ULOC(IAXIS,1)= -SIN(NWAVE*XLOC(IAXIS,1))**2._EB * SIN(2._EB*NWAVE*XLOC(JAXIS,1))
            ULOC(JAXIS,1)=  SIN(2._EB*NWAVE*XLOC(IAXIS,1))  * SIN(NWAVE*XLOC(JAXIS,1))**2._EB

            ! Velocity field in global axes:
            UGLOB        = MATMUL(ROTMAT,ULOC)
            SELECT CASE(X1AXIS)
            CASE(IAXIS)
               VEL_CF = SIN(TLEVEL) * UGLOB(IAXIS,1)
            CASE(JAXIS)
               VEL_CF = 0._EB
            CASE(KAXIS)
               VEL_CF = SIN(TLEVEL) * UGLOB(JAXIS,1)
            END SELECT

            CUT_FACE(ICF)%VEL(ICF2)  = VEL_CF
            CUT_FACE(ICF)%VELS(ICF2) = VEL_CF

         ENDDO

      ELSE ! CC_INBOUNDARY

         VEL_CF = 0._EB
         CUT_FACE(ICF)%VEL(1:NFACE)  = VEL_CF
         CUT_FACE(ICF)%VELS(1:NFACE) = VEL_CF

      ENDIF
   ENDDO CUTFACE_LOOP
ENDIF

! Fields for scalars:
! Regular cells:
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
         ! Global position:
         XGLOB(1:2,1) = (/ XC(I), ZC(K) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         ZZ(I,J,K,N_SPEC_NEUMN) = AMP_Z/3._EB*SIN(TLEVEL)*(1._EB-COS(2._EB*NWAVE*(XLOC(IAXIS,1)-GAM))) * &
                                                          (1._EB-COS(2._EB*NWAVE*(XLOC(JAXIS,1)-GAM)))   &
                                                          -AMP_Z/3._EB*SIN(TLEVEL) + MEAN_Z
         ZZ(I,J,K,N_SPEC_BACKG) = 1._EB - ZZ(I,J,K,N_SPEC_NEUMN)

      ENDDO
   ENDDO
ENDDO

! Cut cells:
IF (CC_IBM) THEN
   CUTCELL_LOOP : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      NCELL  = CUT_CELL(ICC)%NCELL
      DO JCC=1,NCELL
         ! Global position:
         XGLOB(1:2,1) = (/ CUT_CELL(ICC)%XYZCEN(IAXIS,JCC), CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) /)

         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)

         CUT_CELL(ICC)%ZZ(N_SPEC_NEUMN,JCC) = AMP_Z/3._EB*SIN(TLEVEL)*(1._EB-COS(2._EB*NWAVE*(XLOC(IAXIS,1)-GAM))) * &
                                                                      (1._EB-COS(2._EB*NWAVE*(XLOC(JAXIS,1)-GAM)))   &
                                                                      -AMP_Z/3._EB*SIN(TLEVEL) + MEAN_Z
         CUT_CELL(ICC)%ZZ(N_SPEC_BACKG,JCC) = 1._EB - CUT_CELL(ICC)%ZZ(N_SPEC_NEUMN,JCC)
      ENDDO
   ENDDO CUTCELL_LOOP
ENDIF


RETURN

END SUBROUTINE ROTATED_CUBE_ANN_SOLN

! ------------------------------ ROTATED_CUBE_RHS_ZZ -----------------------------------

SUBROUTINE ROTATED_CUBE_RHS_ZZ(TLEVEL,DT,NM)

USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM

REAL(EB), INTENT(IN) :: TLEVEL,DT
INTEGER, INTENT(IN)  :: NM

! Local Variables:
INTEGER :: I, J ,K
REAL(EB), POINTER, DIMENSION(:,:,:)   :: RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: D_Z_N(0:MAX_I_MAX_TEMP),XLOC(2,1),XGLOB(2,1),Q_ZN,D_Z_TEMP,RHO_IJK,DTFC
REAL(EB) :: SIN_T, COS_T

ROTANG = 0._EB
IF(PERIODIC_TEST==21) THEN
   ROTANG = 0._EB ! No rotation.
ELSEIF(PERIODIC_TEST==22) THEN
   ROTANG = ATAN(1._EB/2._EB) ! ~27 Degrees.
ELSEIF(PERIODIC_TEST==23) THEN
   ROTANG = ATAN(1._EB)       ! 45 degrees.
ELSE
   RETURN
ENDIF
ROTMAT(1,1) = COS(ROTANG); ROTMAT(1,2) = -SIN(ROTANG);
ROTMAT(2,1) = SIN(ROTANG); ROTMAT(2,2) =  COS(ROTANG);
TROTMAT = TRANSPOSE(ROTMAT)

SIN_T    = SIN(TLEVEL)
COS_T    = COS(TLEVEL)

CALL POINT_TO_MESH(NM)

IF (PREDICTOR) THEN
   RHOP => RHO
   ZZP  => ZZS
   DTFC = 1._EB
ELSE
   RHOP => RHOS
   ZZP  => ZZ
   DTFC = 0.5_EB
ENDIF

D_Z_N=D_Z(:,N_SPEC_NEUMN)

! Add Q_Z on regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
         IF (CC_IBM) THEN
            IF(CCVAR(I,J,K,CC_CGSC)/=CC_GASPHASE) CYCLE
            IF(CCVAR(I,J,K,CC_UNKZ)>0) CYCLE
         ENDIF
         ! Global position:
         XGLOB(1:2,1) = (/ XC(I), ZC(K) /)
         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)
         RHO_IJK = RHOP(I,J,K)
         CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z_N,1),D_Z_N,TMP(I,J,K),D_Z_TEMP)
         CALL ROTATED_CUBE_NEUMN_FZ(SIN_T,COS_T,RHO_IJK,D_Z_TEMP,XLOC,Q_ZN)

         ! Update species:
         ZZP(I,J,K,N_SPEC_NEUMN) = ZZP(I,J,K,N_SPEC_NEUMN) + DTFC*DT*Q_ZN
         ZZP(I,J,K,N_SPEC_BACKG) = ZZP(I,J,K,N_SPEC_BACKG) - DTFC*DT*Q_ZN
      ENDDO
   ENDDO
ENDDO


RETURN
END SUBROUTINE ROTATED_CUBE_RHS_ZZ

! ------------------------- CC_ROTATED_CUBE_RHS_ZZ -------------------------------

SUBROUTINE CC_ROTATED_CUBE_RHS_ZZ(TLEVEL,N)

USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM

REAL(EB), INTENT(IN) :: TLEVEL
INTEGER, INTENT(IN)  :: N

! Local Variables:
INTEGER :: NM, I, J ,K, IROW, ICC, JCC
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP
REAL(EB) :: PREDFCT,D_Z_N(0:MAX_I_MAX_TEMP),XLOC(2,1),XGLOB(2,1),Q_Z,D_Z_TEMP,RHO_IJK,FCT
REAL(EB) :: SIN_T, COS_T

ROTANG = 0._EB
IF(PERIODIC_TEST==21) THEN
   ROTANG = 0._EB ! No rotation.
ELSEIF(PERIODIC_TEST==22) THEN
   ROTANG = ATAN(1._EB/2._EB) ! ~27 Degrees.
ELSEIF(PERIODIC_TEST==23) THEN
   ROTANG = ATAN(1._EB)       ! 45 degrees.
ELSE
   RETURN
ENDIF
ROTMAT(1,1) = COS(ROTANG); ROTMAT(1,2) = -SIN(ROTANG);
ROTMAT(2,1) = SIN(ROTANG); ROTMAT(2,2) =  COS(ROTANG);
TROTMAT = TRANSPOSE(ROTMAT)

SIN_T    = SIN(TLEVEL)
COS_T    = COS(TLEVEL)

MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   IF (PREDICTOR) THEN
      PREDFCT=1._EB
      RHOP => RHO
   ELSE
      PREDFCT=0._EB
      RHOP => RHOS
   ENDIF

   D_Z_N = D_Z(:,N)

   ! First add Q_Z on regular cells to source F_Z:
   IF (N==N_SPEC_BACKG) THEN
      FCT=-1._EB
   ELSEIF (N==N_SPEC_NEUMN) THEN
      FCT=1._EB
   ENDIF

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
            ! Global position:
            XGLOB(1:2,1) = (/ XC(I), ZC(K) /)
            ! Local position:
            XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
            XLOC         = MATMUL(TROTMAT, XGLOB )
            XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)
            CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z_N,1),D_Z_N,TMP(I,J,K),D_Z_TEMP)
            RHO_IJK = RHOP(I,J,K)
            CALL ROTATED_CUBE_NEUMN_FZ(SIN_T,COS_T,RHO_IJK,D_Z_TEMP,XLOC,Q_Z)
            F_Z(IROW) = F_Z(IROW) - FCT*Q_Z*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   ! Then add Cut-cell contributions to F_Z:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      DO JCC=1,CUT_CELL(ICC)%NCELL
         IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         ! Global position:
         XGLOB(1:2,1) = (/ CUT_CELL(ICC)%XYZCEN(IAXIS,JCC), CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) /)
         ! Local position:
         XGLOB(1:2,1) = XGLOB(1:2,1) - DISPL
         XLOC         = MATMUL(TROTMAT, XGLOB )
         XLOC( 1:2,1) = XLOC( 1:2,1) - DISPXY(1:2,1)
         CALL INTERPOLATE1D_UNIFORM(LBOUND(D_Z_N,1),D_Z_N,CUT_CELL(ICC)%TMP(JCC),D_Z_TEMP)
         RHO_IJK = (1._EB - PREDFCT)*CUT_CELL(ICC)%RHOS(JCC) + PREDFCT*CUT_CELL(ICC)%RHO(JCC)
         CALL ROTATED_CUBE_NEUMN_FZ(SIN_T,COS_T,RHO_IJK,D_Z_TEMP,XLOC,Q_Z)
         F_Z(IROW) = F_Z(IROW) - FCT*Q_Z*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO

ENDDO MESH_LOOP


RETURN
END SUBROUTINE CC_ROTATED_CUBE_RHS_ZZ

! ---------------------------------- ROTATED_CUBE_NEUMN_FZ ---------------------------------

SUBROUTINE ROTATED_CUBE_NEUMN_FZ(SIN_T,COS_T,RHO_IJK,D_Z_TEMP,XLOC,Q_Z)

REAL(EB), INTENT(IN) :: SIN_T,COS_T, RHO_IJK, D_Z_TEMP, XLOC(2,1)
REAL(EB), INTENT(OUT):: Q_Z

! Local Variables:
REAL(EB) :: COS_2KGX, COS_2KGZ, SIN_2KGX, SIN_2KGZ
REAL(EB) :: COS_KX, COS_KZ, SIN_KX, SIN_KZ
REAL(EB) :: SIN_2KX, SIN_2KZ

Q_Z=0._EB

COS_2KGX = COS(2._EB*NWAVE*(GAM - XLOC(IAXIS,1)))
COS_2KGZ = COS(2._EB*NWAVE*(GAM - XLOC(JAXIS,1)))

SIN_2KGX = SIN(2._EB*NWAVE*(GAM - XLOC(IAXIS,1)))
SIN_2KGZ = SIN(2._EB*NWAVE*(GAM - XLOC(JAXIS,1)))

COS_KX   = COS(NWAVE*XLOC(IAXIS,1))
COS_KZ   = COS(NWAVE*XLOC(JAXIS,1))
SIN_KX   = SIN(NWAVE*XLOC(IAXIS,1))
SIN_KZ   = SIN(NWAVE*XLOC(JAXIS,1))

SIN_2KX  = SIN(2._EB*NWAVE*XLOC(IAXIS,1))
SIN_2KZ  = SIN(2._EB*NWAVE*XLOC(JAXIS,1))

Q_Z = (4._EB*AMP_Z*D_Z_TEMP*NWAVE**2*RHO_IJK*COS_2KGX*SIN_T*(COS_2KGZ - 1._EB))/3._EB - &
       RHO_IJK*((AMP_Z*COS_T)/3._EB - (AMP_Z*COS_T*(COS_2KGX - 1._EB)*(COS_2KGZ - 1._EB))/3._EB) + &
      (4._EB*AMP_Z*D_Z_TEMP*NWAVE**2*RHO_IJK*COS_2KGZ*SIN_T*(COS_2KGX - 1._EB))/3._EB - &
       2._EB*NWAVE*RHO_IJK*COS_KX*SIN_KX*SIN_2KZ*SIN_T*(MEAN_Z - &
      (AMP_Z*SIN_T)/3._EB + (AMP_Z*SIN_T*(COS_2KGX - 1._EB)*(COS_2KGZ - 1._EB))/3._EB) + &
       2._EB*NWAVE*RHO_IJK*COS_KZ*SIN_2KX*SIN_KZ*SIN_T*(MEAN_Z - (AMP_Z*SIN_T)/3._EB + &
      (AMP_Z*SIN_T*(COS_2KGX - 1._EB)*(COS_2KGZ - 1._EB))/3._EB) + &
      (2._EB*AMP_Z*NWAVE*RHO_IJK*SIN_2KX*SIN_KZ**2*SIN_2KGZ*SIN_T**2*(COS_2KGX - 1._EB))/3._EB - &
      (2._EB*AMP_Z*NWAVE*RHO_IJK*SIN_KX**2*SIN_2KZ*SIN_2KGX*SIN_T**2*(COS_2KGZ - 1._EB))/3._EB


RETURN
END SUBROUTINE ROTATED_CUBE_NEUMN_FZ

END MODULE CC_VERIFICATION
