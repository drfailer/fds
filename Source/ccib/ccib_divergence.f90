!  +++++++++++++++++++++++ CC_DIVERGENCE ++++++++++++++++++++++++++

! Divergence computation routines for the cut-cell / immersed-boundary method.

MODULE CC_DIVERGENCE

USE CC_SCALARS_DATA
USE CC_DIVERGENCE_KERNELS
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CC_DIVERGENCE_PART_1, CC_CHECK_DIVERGENCE, &
          CC_DIFFUSIVE_MASS_FLUXES, &
          GET_CC_CELL_DIFFUSIVITY, GET_CC_CELL_CONDUCTIVITY, &
          GET_VELOC_DIVERGENCE_CUTCELL, GET_FN_DIVERGENCE_CUTCELL, &
          SET_EXIMADVFLX_3D, SET_EXIMDIFFLX_3D, &
          SET_EXIMRHOHSLIM_3D, SET_EXIMRHOZZLIM_3D

CONTAINS


! -------------------------- CC_DIVERGENCE_PART_1 --------------------------

SUBROUTINE CC_DIVERGENCE_PART_1(T,DT,NM)

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_Z, &
                              GET_SENSIBLE_ENTHALPY,GET_VISCOSITY,GET_MOLECULAR_WEIGHT
USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_Z_SRC

REAL(EB), INTENT(IN) :: T,DT
INTEGER,  INTENT(IN) :: NM
! Recompute divergence terms in cut-cell region and surrounding cells.
! Use velocity divergence equivalence to define divergence on cut-cell underlying Cartesian cells.

! Local Variables:
INTEGER :: N,I,J,K,X1AXIS,ISIDE,IFACE,ICC,JCC,ICF
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP,DPVOL,RHOP,RTRM,CP,R_H_G,U_DOT_DEL_RHO_Z_VOL
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: RDT,CCM1,CCP1,IDX,AF,TMP_G,H_S,TNOW,RHOPV(-2:1),TMPV(-1:0),X1F,PRFCT,PRFCTV, &
            CPV(-1:0),FCT,MUV(-1:0),MU_DNSV(-1:0)
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM

REAL(EB) :: VCELL, VCCELL, DIVVOL, DUMMY, RTRMVOL, CCVOL

LOGICAL, PARAMETER :: DO_CONDUCTION_HEAT_FLUX=.TRUE.
INTEGER :: DIFFHFLX_IND, JFLX_IND

LOGICAL, PARAMETER :: SET_DIV_TO_ZERO  = .FALSE.
LOGICAL, PARAMETER :: SET_CCDIV_TO_ZERO= .FALSE.
LOGICAL, PARAMETER :: FIX_DIFF_FLUXES  = .TRUE.

REAL(EB), ALLOCATABLE, DIMENSION(:) :: DIVRG_VEC , RTRM_VEC, VOLDVRG
INTEGER :: INDZ

! Pressure sums re-integration vars:
INTEGER :: IW,IND1,IND2

! Shunn MMS test case vars:
REAL(EB) :: XHAT, ZHAT, Q_Z, TT

! Dummy on T:
DUMMY = T

! Check whether to skip this routine

IF (SOLID_PHASE_ONLY) RETURN

TNOW=CURRENT_TIME()

DIFFHFLX_IND = LOW_IND  ! -rho Da Grad(Za)
JFLX_IND     = LOW_IND

CALL POINT_TO_MESH(NM)

RDT = 1._EB/DT

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      DP     => DS
      PBAR_P => PBAR_S
      RHOP   => RHOS
      PRFCT  = 0._EB ! Use star cut-cell quantities.
   CASE(.FALSE.)
      DP     => D
      PBAR_P => PBAR
      RHOP   => RHO
      PRFCT  = 1._EB ! Use end of step cut-cell quantities.
END SELECT


R_PBAR = 1._EB/PBAR_P
DPVOL  => DP
RTRM   => WORK1

! Set DP to zero in Cartesian cells of type: CC_SOLID, CC_CUTCFE, and CC_GASPHASE where CC_UNKZ > 0:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF ((CCVAR(I,J,K,CC_CGSC) == CC_GASPHASE) .AND. (CCVAR(I,J,K,CC_UNKZ) <= 0)) CYCLE
         DPVOL(I,J,K) = 0._EB
         DEL_RHO_D_DEL_Z(I,J,K,1:N_TOTAL_SCALARS) = 0._EB
      ENDDO
   ENDDO
ENDDO
DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   CUT_CELL(ICC)%DVOL(1:CUT_CELL(ICC)%NCELL)= 0._EB
   CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,1:CUT_CELL(ICC)%NCELL)=0._EB
ENDDO
IF (CORRECTOR) THEN
   IF (ALLOCATED(MESHES(NM)%D_SOURCE)) THEN
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         I     = CUT_CELL(ICC)%IJK(IAXIS)
         J     = CUT_CELL(ICC)%IJK(JAXIS)
         K     = CUT_CELL(ICC)%IJK(KAXIS)
         VCELL = DX(I)*DY(J)*DZ(K)

         ! Up to here in D_SOURCE(I,J,K), M_DOT_PPP(I,J,K,1:N_TOTAL_SCALARS) we have contributions by particles.
         ! Add these contributions in corresponding cut-cells:
         ! NOTE : Assumes the source from particles is distributed evenly over CCs of the Cartesian cell.
         VCCELL = SUM(CUT_CELL(ICC)%VOLUME(1:CUT_CELL(ICC)%NCELL))
         DO JCC=1,CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%D_SOURCE(JCC) = CUT_CELL(ICC)%D_SOURCE(JCC) + D_SOURCE(I,J,K)
            CUT_CELL(ICC)%M_DOT_PPP(1:N_TOTAL_SCALARS,JCC) = &
            CUT_CELL(ICC)%M_DOT_PPP(1:N_TOTAL_SCALARS,JCC) + M_DOT_PPP(I,J,K,1:N_TOTAL_SCALARS)
         ENDDO
      ENDDO
   ENDIF
ENDIF

IF (SET_DIV_TO_ZERO) THEN
   DP = 0._EB ! Set to zero divg on all cells.
   RETURN
ENDIF
IF (SET_CCDIV_TO_ZERO) RETURN

! Point to corresponding ZZ array:
SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      ZZP => ZZS
   CASE(.FALSE.)
      ZZP => ZZ
END SELECT

ALLOCATE(ZZ_GET(N_TRACKED_SPECIES))

! Add species diffusion terms to divergence expression and compute diffusion term for species equations
SPECIES_GT_1_IF: IF (N_TOTAL_SCALARS>1) THEN

   ! 1. Diffusive Heat flux = - Grad dot (h_s rho D Grad Z_n):
   ! In FV form: use faces to add corresponding face integral terms, for face k
   ! (sum_a{h_{s,a} rho D_a Grad z_a) dot \hat{n}_k A_k, where \hat{n}_k is the versor outside of cell
   ! at face k.
   CALL CC_DIFFUSIVE_MASS_FLUXES(NM)

   ! Ensure RHO_D terms sum to zero over all species.  Gather error into largest mass fraction present.
   IF (FIX_DIFF_FLUXES) CALL FIX_CC_DIFF_MASS_FLUXES

   ! Zero out DEL_RHO_D_DEL_Z for impregion regular cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            DEL_RHO_D_DEL_Z(I,J,K,1:N_TOTAL_SCALARS) = 0._EB
         ENDDO
      ENDDO
   ENDDO

   ! 1. Diffusive heat flux  = - hs,a (Da Grad(rho*Ya) - Da/rho Grad(rho) (rho Ya)):
   CALL CC_DIFFUSIVE_HEAT_FLUXES

ENDIF SPECIES_GT_1_IF


CONDUCTION_HEAT_IF : IF( DO_CONDUCTION_HEAT_FLUX ) THEN
   ! 2. Conduction heat flux = - k Grad(T):
   CALL CC_CONDUCTION_HEAT_FLUX
ENDIF CONDUCTION_HEAT_IF


! Add \dot{q}''' and QR to DP:
! Regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         ! Add \dot{q}''' and QR to DP*Vii:
         DPVOL(I,J,K) = DPVOL(I,J,K) + (Q(I,J,K) + QR(I,J,K)) * DX(I)*DY(J)*DZ(K)
      ENDDO
   ENDDO
ENDDO

! HERE Cut-cells \dot{q}'''*VOL and QR*VOL:
DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   DO JCC=1,CUT_CELL(ICC)%NCELL
      CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC)+(CUT_CELL(ICC)%Q(JCC)+CUT_CELL(ICC)%QR(JCC))*CUT_CELL(ICC)%VOLUME(JCC)
   ENDDO
ENDDO

! 3. Enthalpy advection term = - \bar{ u dot Grad (rho h_s) }:
! R_H_G = 1/(Cp * T)
! RTRM  = 1/(rho * Cp * T)
! Point to the appropriate velocity components

IF (PREDICTOR) THEN
   UU=>U
   VV=>V
   WW=>W
   PRFCTV = 1._EB
ELSE
   UU=>US
   VV=>VS
   WW=>WS
   PRFCTV = 0._EB
ENDIF

CONST_GAMMA_IF_1: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
   CALL CCENTHALPY_ADVECTION ! Compute u dot grad rho h_s in FV form and add to DP in regular + cut-cells.
ENDIF CONST_GAMMA_IF_1


! Loop through regular cells in the implicit region, as well as cut-cells and compute R_H_G, and RTRM:
CP    => WORK5
R_H_G => WORK9
RTRM  => WORK1
! Regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K,1:N_TRACKED_SPECIES)
         CALL GET_SPECIFIC_HEAT(ZZ_GET,CP(I,J,K),TMP(I,J,K))
         R_H_G(I,J,K) = 1._EB/(CP(I,J,K)*TMP(I,J,K))
         RTRM(I,J,K)  = R_H_G(I,J,K)/RHOP(I,J,K)
         DPVOL(I,J,K) = RTRM(I,J,K)*DPVOL(I,J,K)
      ENDDO
   ENDDO
ENDDO

! Cut-cells:
DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
   DO JCC=1,CC%NCELL
      TMPV(0) = CC%TMP(JCC)
      ZZ_GET(1:N_TRACKED_SPECIES) = PRFCT*CC%ZZ(1:N_TRACKED_SPECIES,JCC) + (1._EB-PRFCT)*CC%ZZS(1:N_TRACKED_SPECIES,JCC)
      CALL GET_SPECIFIC_HEAT(ZZ_GET,CPV(0),TMPV(0))
      CC%R_H_G(JCC) = 1._EB/(CPV(0)*TMPV(0))
      RHOPV(0) = PRFCT *CC%RHO(JCC) + (1._EB-PRFCT)*CC%RHOS(JCC)
      CC%RTRM(JCC) = CC%R_H_G(JCC)/RHOPV(0)
      CC%DVOL(JCC) = CC%RTRM(JCC)*CC%DVOL(JCC)
   ENDDO
ENDDO


! 4. Enthalpy flux due to mass diffusion and advection:
! sum_n [\bar{W}/W_n - h_{s,n}*R_H_G] ( Grad dot (rho D_\alpha Grad Z_n) - \bar{u dot Grad (rho Z_n)})

CONST_GAMMA_IF_2: IF (.NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN

   SPECIES_LOOP: DO N=1,N_TRACKED_SPECIES

      CALL CCSPECIES_ADVECTION ! Compute u dot grad rho Z_n

      SM  => SPECIES_MIXTURE(N)

      ! Regular cells:
      ICC = 0
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP(I,J,K),H_S)
               DPVOL(I,J,K) = DPVOL(I,J,K) + (SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K))* &
                                             (DEL_RHO_D_DEL_Z(I,J,K,N) - U_DOT_DEL_RHO_Z_VOL(I,J,K))/RHOP(I,J,K)
               ! Values of DEL_RHO_D_DEL_Z(I,J,K,N) have been filled previously.
               ! RSUM was computed in the implicit region advance routine for scalars CCDENSITY.
               ICC = ICC + 1
            ENDDO
         ENDDO
      ENDDO

      ! Cut-cells:
      IF (PREDICTOR) THEN
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL
               TMPV(0) = CC%TMP(JCC)
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMPV(0),H_S)
               CC%DVOL(JCC) = CC%DVOL(JCC) + (SM%RCON/CC%RSUM(JCC) - H_S*CC%R_H_G(JCC))/CC%RHOS(JCC) * &
                                             (CC%DEL_RHO_D_DEL_Z_VOL(N,JCC)- CC%U_DOT_DEL_RHO_Z_VOL(N,JCC))
            ENDDO
         ENDDO
      ELSE
         DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
            CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
            DO JCC=1,CC%NCELL
               TMPV(0) = CC%TMP(JCC)
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMPV(0),H_S)
               CC%DVOL(JCC) = CC%DVOL(JCC) + (SM%RCON/CC%RSUM(JCC) - H_S*CC%R_H_G(JCC))/CC%RHO(JCC) * &
                                             (CC%DEL_RHO_D_DEL_Z_VOL(N,JCC)- CC%U_DOT_DEL_RHO_Z_VOL(N,JCC))
            ENDDO
         ENDDO
      ENDIF

   ENDDO SPECIES_LOOP

ENDIF CONST_GAMMA_IF_2

! Add contribution of reactions

IF (ALLOCATED(MESHES(NM)%D_SOURCE)) THEN

   ! Regular Cells on the implicit region:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            DPVOL(I,J,K) = DPVOL(I,J,K) + D_SOURCE(I,J,K)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   ! Cut cells:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         CC%DVOL(JCC) = CC%DVOL(JCC) + CC%D_SOURCE(JCC)*CC%VOLUME(JCC)
      ENDDO
   ENDDO

ENDIF

! Atmospheric stratification term

IF (STRATIFICATION) THEN
   ! Regular Cells on the implicit region:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            DPVOL(I,J,K) = DPVOL(I,J,K) + &
            RTRM(I,J,K)*0.5_EB*(WW(I,J,K)+WW(I,J,K-1))* &
            RHO_0_CV( CCVAR(I,J,K,CC_UNKZ)-UNKZ_IND(NM_START) )*GVEC(KAXIS)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         ! D = D + w*rho_0*g/(rho*Cp*T)*Vii
         CC%DVOL(JCC) = CC%DVOL(JCC) + CC%RTRM(JCC)*CC%WVEL(JCC)*CC%RHO_0(JCC)*GVEC(KAXIS)*CC%VOLUME(JCC)
      ENDDO
   ENDDO
ENDIF

! Manufactured solution

MMS_IF: IF (PERIODIC_TEST==7) THEN
   IF (PREDICTOR) TT=T+DT
   IF (CORRECTOR) TT=T
   ! Regular cells on cut-cell region:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            ! this term is similar to D_REACTION from fire
            XHAT = XC(I) - UF_MMS*TT
            ZHAT = ZC(K) - WF_MMS*TT
            DO N=1,N_TRACKED_SPECIES
               SM => SPECIES_MIXTURE(N)
               SELECT CASE(N)
                  CASE(1); Q_Z = -VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
                  CASE(2); Q_Z =  VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
               END SELECT
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP(I,J,K),H_S)
               DPVOL(I,J,K) = DPVOL(I,J,K) + ( SM%RCON/RSUM(I,J,K) - H_S*R_H_G(I,J,K) )*Q_Z/RHOP(I,J,K)*DX(I)*DY(J)*DZ(K)
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   ! Cut-cells:
   IF (PREDICTOR) THEN
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         DO JCC=1,CUT_CELL(ICC)%NCELL
            ! this term is similar to D_REACTION from fire
            XHAT = CUT_CELL(ICC)%XYZCEN(IAXIS,JCC) - UF_MMS*TT
            ZHAT = CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) - WF_MMS*TT
            TMPV(0) = CUT_CELL(ICC)%TMP(JCC)
            DO N=1,N_TRACKED_SPECIES
               SM => SPECIES_MIXTURE(N)
               SELECT CASE(N)
                  CASE(1); Q_Z = -VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
                  CASE(2); Q_Z =  VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
               END SELECT
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMPV(0),H_S)
               CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) +  &
               (SM%RCON/CUT_CELL(ICC)%RSUM(JCC) - H_S*CUT_CELL(ICC)%R_H_G(JCC)) * &
               Q_Z/CUT_CELL(ICC)%RHOS(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
            ENDDO
         ENDDO
      ENDDO
   ELSE ! CORRECTOR
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         DO JCC=1,CUT_CELL(ICC)%NCELL
            ! this term is similar to D_REACTION from fire
            XHAT = CUT_CELL(ICC)%XYZCEN(IAXIS,JCC) - UF_MMS*TT
            ZHAT = CUT_CELL(ICC)%XYZCEN(KAXIS,JCC) - WF_MMS*TT
            TMPV(0) = CUT_CELL(ICC)%TMP(JCC)
            DO N=1,N_TRACKED_SPECIES
               SM => SPECIES_MIXTURE(N)
               SELECT CASE(N)
                  CASE(1); Q_Z = -VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
                  CASE(2); Q_Z =  VD2D_MMS_Z_SRC(XHAT,ZHAT,TT)
               END SELECT
               CALL GET_SENSIBLE_ENTHALPY_Z(N,TMPV(0),H_S)
               CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) +  &
               (SM%RCON/CUT_CELL(ICC)%RSUM(JCC) - H_S*CUT_CELL(ICC)%R_H_G(JCC)) * &
               Q_Z/CUT_CELL(ICC)%RHO(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF MMS_IF

! Assign divergence and 1/(rho*Cp*T) on Cartesian Cells:
! Average divergence on linked cells:
ALLOCATE ( DIVRG_VEC(1:NUNKZ_LOCAL) , VOLDVRG(1:NUNKZ_LOCAL), RTRM_VEC(1:NUNKZ_LOCAL) )
DIVRG_VEC(:) = 0._EB
VOLDVRG(:)   = 0._EB
RTRM_VEC(:)  = 0._EB

! Add div*vol for all cells and cut-cells on implicit region:
! Regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         ! Unknown number:
         INDZ  = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
         DIVRG_VEC(INDZ) =  DIVRG_VEC(INDZ) + DPVOL(I,J,K)
         RTRM_VEC(INDZ)  =  RTRM_VEC(INDZ)  + RTRM(I,J,K)*(DX(I)*DY(J)*DZ(K))
         VOLDVRG(INDZ)   =  VOLDVRG(INDZ)   + (DX(I)*DY(J)*DZ(K))
      ENDDO
   ENDDO
ENDDO

If (PREDICTOR) THEN
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I      = CUT_CELL(ICC)%IJK(IAXIS)
      J      = CUT_CELL(ICC)%IJK(JAXIS)
      K      = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO JCC=1,CUT_CELL(ICC)%NCELL
         INDZ = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         DIVRG_VEC(INDZ) =  DIVRG_VEC(INDZ) + CUT_CELL(ICC)%DVOL(JCC)
         RTRM_VEC(INDZ)  =  RTRM_VEC(INDZ)  + CUT_CELL(ICC)%RTRM(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
         VOLDVRG(INDZ)   =  VOLDVRG(INDZ)   + CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO
ELSE ! CORRECTOR
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I      = CUT_CELL(ICC)%IJK(IAXIS)
      J      = CUT_CELL(ICC)%IJK(JAXIS)
      K      = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO JCC=1,CUT_CELL(ICC)%NCELL
         INDZ = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         DIVRG_VEC(INDZ) =  DIVRG_VEC(INDZ) + CUT_CELL(ICC)%DVOL(JCC)
         RTRM_VEC(INDZ)  =  RTRM_VEC(INDZ)  + CUT_CELL(ICC)%RTRM(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
         VOLDVRG(INDZ)   =  VOLDVRG(INDZ)   + CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO
   ENDDO
ENDIF

! Here there should be a mesh exchange (add) of div*vol for cases where cut-cells are linked to cells
! that belong to other meshes.

! Compute final divergence:
DO INDZ=UNKZ_ILC(NM)+1,UNKZ_ILC(NM)+NUNKZ_LOC(NM)
   DIVRG_VEC(INDZ)=DIVRG_VEC(INDZ)/VOLDVRG(INDZ)
   RTRM_VEC(INDZ) = RTRM_VEC(INDZ)/VOLDVRG(INDZ)
ENDDO

! Finally load final thermodynamic divergence to corresponding cells:
! Regular cells:
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
         INDZ  = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
         DP(I,J,K)   = DIVRG_VEC(INDZ) ! Previously divided by VOL.
         RTRM(I,J,K) = RTRM_VEC(INDZ)  ! Previously divided by VOL.
         DEL_RHO_D_DEL_Z(I,J,K,1:N_TRACKED_SPECIES) = DEL_RHO_D_DEL_Z(I,J,K,1:N_TRACKED_SPECIES)/(DX(I)*DY(J)*DZ(K))
      ENDDO
   ENDDO
ENDDO

IF (PREDICTOR) THEN
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I      = CUT_CELL(ICC)%IJK(IAXIS)
      J      = CUT_CELL(ICC)%IJK(JAXIS)
      K      = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DIVVOL = 0._EB
      RTRMVOL= 0._EB
      CCVOL  = 0._EB
      DO JCC=1,CUT_CELL(ICC)%NCELL
         INDZ = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         CUT_CELL(ICC)%DVOL(JCC)= DIVRG_VEC(INDZ)*CUT_CELL(ICC)%VOLUME(JCC)
         CUT_CELL(ICC)%DS(JCC)  = DIVRG_VEC(INDZ)
         CUT_CELL(ICC)%RTRM(JCC)= RTRM_VEC(INDZ)
         DIVVOL = DIVVOL + CUT_CELL(ICC)%DVOL(JCC)
         RTRMVOL= RTRMVOL+ CUT_CELL(ICC)%RTRM(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
         CCVOL  = CCVOL  + CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO

      ! Now get sum(un*ACFace) and add to divergence:
      DP(I,J,K)  = DIVVOL/(DX(I)*DY(J)*DZ(K)) ! Now push Divergence to underlying Cartesian cell.
      RTRM(I,J,K)= RTRMVOL/(DX(I)*DY(J)*DZ(K))
   ENDDO
ELSE ! CORRECTOR
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I      = CUT_CELL(ICC)%IJK(IAXIS)
      J      = CUT_CELL(ICC)%IJK(JAXIS)
      K      = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DIVVOL = 0._EB
      RTRMVOL= 0._EB
      CCVOL  = 0._EB
      DO JCC=1,CUT_CELL(ICC)%NCELL
         INDZ = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         CUT_CELL(ICC)%DVOL(JCC)= DIVRG_VEC(INDZ)*CUT_CELL(ICC)%VOLUME(JCC)
         CUT_CELL(ICC)%D(JCC)   = DIVRG_VEC(INDZ)
         CUT_CELL(ICC)%RTRM(JCC)= RTRM_VEC(INDZ)
         DIVVOL = DIVVOL + CUT_CELL(ICC)%DVOL(JCC)
         RTRMVOL= RTRMVOL+ CUT_CELL(ICC)%RTRM(JCC)*CUT_CELL(ICC)%VOLUME(JCC)
         CCVOL  = CCVOL  + CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO

      ! Now get sum(un*ACFace) and add to divergence:
      DP(I,J,K) = DIVVOL/(DX(I)*DY(J)*DZ(K))
      RTRM(I,J,K)= RTRMVOL/(DX(I)*DY(J)*DZ(K))
   ENDDO
ENDIF
DEALLOCATE ( DIVRG_VEC , VOLDVRG, RTRM_VEC )
DEALLOCATE(ZZ_GET)

T_USED(14) = T_USED(14) + CURRENT_TIME() - TNOW
IF (TIME_CC_IBM) &
   T_CC_USED(CC_DIVERGENCE_PART_1_TIME_INDEX) = T_CC_USED(CC_DIVERGENCE_PART_1_TIME_INDEX) + CURRENT_TIME() - TNOW
RETURN

CONTAINS

! --------------------------- FIX_CC_DIFF_MASS_FLUXES -------------------------

SUBROUTINE FIX_CC_DIFF_MASS_FLUXES

REAL(EB) :: ZZ_FACE(1:N_TRACKED_SPECIES)

! IAXIS faces:
X1AXIS = IAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC
   ! Cycle if boundary condition other then INTERPOLATED, OPEN or PERIODIC, already done in GET_BBREGFACE_RHO_D_DZDN.
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY         .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(KAXIS)

   ZZ_FACE(1:N_TRACKED_SPECIES) = 0.5_EB*(ZZP(I+1,J,K,1:N_TRACKED_SPECIES) + &
                                          ZZP(I  ,J,K,1:N_TRACKED_SPECIES))

   N=MAXLOC(ZZ_FACE(1:N_TRACKED_SPECIES),1)

   CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N) = &
   -(SUM(CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(1:N_TRACKED_SPECIES))-CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N))

ENDDO

! JAXIS faces:
X1AXIS = JAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC
   ! Cycle if boundary condition other then INTERPOLATED, OPEN or PERIODIC, already done in GET_BBREGFACE_RHO_D_DZDN.
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY         .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(KAXIS)

   ZZ_FACE(1:N_TRACKED_SPECIES) = 0.5_EB*(ZZP(I,J+1,K,1:N_TRACKED_SPECIES) + &
                                          ZZP(I,J  ,K,1:N_TRACKED_SPECIES))

   N=MAXLOC(ZZ_FACE(1:N_TRACKED_SPECIES),1)

   CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N) = &
   -(SUM(CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(1:N_TRACKED_SPECIES))-CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N))

ENDDO

! KAXIS faces:
X1AXIS = KAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC
   ! Cycle if boundary condition other then INTERPOLATED, OPEN or PERIODIC, already done in GET_BBREGFACE_RHO_D_DZDN.
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY         .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(KAXIS)

   ZZ_FACE(1:N_TRACKED_SPECIES) = 0.5_EB*(ZZP(I,J,K+1,1:N_TRACKED_SPECIES) + &
                                          ZZP(I,J,K  ,1:N_TRACKED_SPECIES))

   N=MAXLOC(ZZ_FACE(1:N_TRACKED_SPECIES),1)

   CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N) = &
   -(SUM(CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(1:N_TRACKED_SPECIES))-CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N))

ENDDO

! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z
   IW = MESHES(NM)%RC_FACE(IFACE)%IWC
   ! Cycle if boundary condition other then INTERPOLATED, OPEN or PERIODIC, already done in GET_BBRCFACE_RHO_D_DZDN.
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY         .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   ZZ_FACE(1:N_TRACKED_SPECIES) = RC_FACE(IFACE)%ZZ_FACE(1:N_TRACKED_SPECIES)
   N=MAXLOC(ZZ_FACE(1:N_TRACKED_SPECIES),1)

   RC_FACE(IFACE)%RHO_D_DZDN(N) = -(SUM(RC_FACE(IFACE)%RHO_D_DZDN(1:N_TRACKED_SPECIES))- &
                                         RC_FACE(IFACE)%RHO_D_DZDN(N))

ENDDO


! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE
   IW = MESHES(NM)%CUT_FACE(ICF)%IWC
   ! Cycle if boundary condition other then INTERPOLATED, OPEN or PERIODIC, already done in GET_BBCUTFACE_RHO_D_DZDN.
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY         .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   DO IFACE=1,CUT_FACE(ICF)%NFACE
      ZZ_FACE(1:N_TRACKED_SPECIES) = CUT_FACE(ICF)%ZZ_FACE(1:N_TRACKED_SPECIES,IFACE)

      N=MAXLOC(ZZ_FACE(1:N_TRACKED_SPECIES),1)
      CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE) = &
      -(SUM(CUT_FACE(ICF)%RHO_D_DZDN(1:N_TRACKED_SPECIES,IFACE))-CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE))

   ENDDO ! IFACE
ENDDO ! ICF

END SUBROUTINE FIX_CC_DIFF_MASS_FLUXES


! ---------------------------- CCSPECIES_ADVECTION ------------------------------

SUBROUTINE CCSPECIES_ADVECTION


! Computes FV version of flux limited \bar{u dot Grad rho Yalpha} in faces near IB
! region and adds components to thermodynamic divergence.

! Local Variables:
REAL(EB) :: RHO_Z_PV(-2:1), VELC, FN_ZZ, ZZ_GET_N
REAL(EB), PARAMETER :: SGNFCT=1._EB
INTEGER :: IOR, ICFA
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(CC_REGFACEZ_TYPE),  POINTER, DIMENSION(:) :: RGF

U_DOT_DEL_RHO_Z_VOL=>WORK7
U_DOT_DEL_RHO_Z_VOL=0._EB
U_TEMP => U_WORK
F_TEMP => F_WORK
Z_TEMP => Z_WORK

! Zero out  for species N in cut-cells:
DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
   CC => CUT_CELL(ICC); CC%U_DOT_DEL_RHO_Z_VOL(N,1:CC%NCELL) = 0._EB
ENDDO

! IAXIS faces:
X1AXIS = IAXIS; RGF => CC_REGFACE_IAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW=RGF(IFACE)%IWC; I =RGF(IFACE)%IJK(IAXIS); J=RGF(IFACE)%IJK(JAXIS); K =RGF(IFACE)%IJK(KAXIS); AF = DY(J)*DZ(K)
   RG_ON_WC_IF_1 : IF((IW>0).AND. .NOT.ANY(WALL(IW)%BOUNDARY_TYPE==(/INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/))) THEN
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
      !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
      ISIDE = -1 + (SIGN(1,IOR)+1) / 2
      RHOPV(ISIDE)    = RHOP(I+1+ISIDE,J,K)
      RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZP(I+1+ISIDE,J,K,N)
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT
            VELC = UU(I,J,K)
         CASE(SOLID_BOUNDARY)
            IF (PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            IF (CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      END SELECT
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K) = U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K) - &
                                           SIGN(1._EB,REAL(IOR,EB))*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
   ELSE RG_ON_WC_IF_1
      VELC=UU(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
      ! Get rho*zz on cells at both sides of IFACE:
      RHO_Z_PV(-1:0) = RHOP(I:I+1,J,K)*ZZP(I:I+1,J,K,N)
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      IF(RGF(IFACE)%DO_LO_IND) U_DOT_DEL_RHO_Z_VOL(I  ,J,K) = U_DOT_DEL_RHO_Z_VOL(I  ,J,K) + &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV(-1))*VELC*AF !+ve dot
      IF(RGF(IFACE)%DO_HI_IND) U_DOT_DEL_RHO_Z_VOL(I+1,J,K) = U_DOT_DEL_RHO_Z_VOL(I+1,J,K) - &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV( 0))*VELC*AF !-ve dot
   ENDIF RG_ON_WC_IF_1
ENDDO

! JAXIS faces:
X1AXIS = JAXIS; RGF => CC_REGFACE_JAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW=RGF(IFACE)%IWC; I =RGF(IFACE)%IJK(IAXIS); J=RGF(IFACE)%IJK(JAXIS); K =RGF(IFACE)%IJK(KAXIS); AF = DX(I)*DZ(K)
   RG_ON_WC_IF_2 : IF((IW>0).AND. .NOT.ANY(WALL(IW)%BOUNDARY_TYPE==(/INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/))) THEN
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC    => BOUNDARY_COORD(WC%BC_INDEX)
      IOR = BC%IOR
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
      !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
      ISIDE = -1 + (SIGN(1,IOR)+1) / 2
      RHOPV(ISIDE)    = RHOP(I,J+1+ISIDE,K)
      RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZP(I,J+1+ISIDE,K,N)
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT
            VELC = VV(I,J,K)
         CASE(SOLID_BOUNDARY)
            IF (PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            IF (CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      END SELECT
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K) = U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K) - &
                                           SIGN(1._EB,REAL(IOR,EB))*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
   ELSE RG_ON_WC_IF_2
      VELC=VV(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
      ! Get rho*zz on cells at both sides of IFACE:
      RHO_Z_PV(-1:0) = RHOP(I,J:J+1,K)*ZZP(I,J:J+1,K,N)
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      IF(RGF(IFACE)%DO_LO_IND) U_DOT_DEL_RHO_Z_VOL(I,J  ,K) = U_DOT_DEL_RHO_Z_VOL(I,J  ,K) + &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV(-1))*VELC*AF !+ve dot
      IF(RGF(IFACE)%DO_HI_IND) U_DOT_DEL_RHO_Z_VOL(I,J+1,K) = U_DOT_DEL_RHO_Z_VOL(I,J+1,K) - &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV( 0))*VELC*AF !-ve dot
   ENDIF RG_ON_WC_IF_2
ENDDO

! KAXIS faces:
X1AXIS = KAXIS; RGF => CC_REGFACE_KAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW=RGF(IFACE)%IWC; I =RGF(IFACE)%IJK(IAXIS); J=RGF(IFACE)%IJK(JAXIS); K =RGF(IFACE)%IJK(KAXIS); AF = DX(I)*DY(J)
   RG_ON_WC_IF_3 : IF((IW>0).AND. .NOT.ANY(WALL(IW)%BOUNDARY_TYPE==(/INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/))) THEN
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR
      ! This expression is such that when sign of IOR is -1 ->G use Low Side cell  -> ISIDE=-1,
      !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
      ISIDE = -1 + (SIGN(1,IOR)+1) / 2
      RHOPV(ISIDE)    = RHOP(I,J,K+1+ISIDE)
      RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZP(I,J,K+1+ISIDE,N)
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT
            VELC = WW(I,J,K)
         CASE(SOLID_BOUNDARY)
            IF (PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            IF (CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
      END SELECT
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE) = U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE) - &
                                           SIGN(1._EB,REAL(IOR,EB))*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
   ELSE RG_ON_WC_IF_3
      VELC=WW(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
      ! Get rho*zz on cells at both sides of IFACE:
      RHO_Z_PV(-1:0) = RHOP(I,J,K:K+1)*ZZP(I,J,K:K+1,N)
      ! Add: -(bar{rho*zz} u dot n - (rho*zz) u dot n) to corresponding cell DP:
      IF(RGF(IFACE)%DO_LO_IND) U_DOT_DEL_RHO_Z_VOL(I,J,K  ) = U_DOT_DEL_RHO_Z_VOL(I,J,K  ) + &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV(-1))*VELC*AF !+ve dot
      IF(RGF(IFACE)%DO_HI_IND) U_DOT_DEL_RHO_Z_VOL(I,J,K+1) = U_DOT_DEL_RHO_Z_VOL(I,J,K+1) - &
                                                         SGNFCT*(RGF(IFACE)%FN_ZZ(N)-RHO_Z_PV( 0))*VELC*AF !-ve dot
   ENDIF RG_ON_WC_IF_3
ENDDO

! Regular faces connecting gasphase - cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z
   RCF => RC_FACE(IFACE); I =RCF%IJK(IAXIS); J =RCF%IJK(JAXIS); K =RCF%IJK(KAXIS); X1AXIS =RCF%IJK(KAXIS+1); IW =RCF%IWC
   RCF_ON_WALL_CELL_IF : IF((IW > 0)) THEN ! INTERPOLATED or PERIODIC treated through RHO_F, ZZ_F(N).
      IF(WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
      !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
      ISIDE = -1 + (SIGN(1,IOR)+1) / 2
      FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
      ! First (rho hs)_i,j,k:
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF = DY(J)*DZ(K)
         VELC = UU(I,J,K)
         RHOPV(ISIDE)      = RHOP(I+1+ISIDE,J,K)
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
            ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
            ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
         END SELECT
      CASE(JAXIS)
         AF = DX(I)*DZ(K)
         VELC = VV(I,J,K)
         RHOPV(ISIDE)      = RHOP(I,J+1+ISIDE,K)
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
            ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
            ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
         END SELECT
      CASE(KAXIS)
         AF = DX(I)*DY(J)
         VELC = WW(I,J,K)
         RHOPV(ISIDE)      = RHOP(I,J,K+1+ISIDE)
         SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
            ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = RCF%CELL_LIST(2,ISIDE+2)
            JCC = RCF%CELL_LIST(3,ISIDE+2)
            RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
            ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
         END SELECT
      END SELECT
      RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
      FN_ZZ           = B1%RHO_F*B1%ZZ_F(N)
      SELECT CASE(WC%BOUNDARY_TYPE)
         CASE DEFAULT
            ! Already filled in previous X1AXIS select case.
         CASE(SOLID_BOUNDARY)
            IF (PREDICTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL_S
            IF (CORRECTOR) VELC = -SIGN(1._EB,REAL(IOR,EB))*B1%U_NORMAL
         CASE(INTERPOLATED_BOUNDARY)
            VELC = UVW_SAVE(IW)
      END SELECT
      SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
      CASE(CC_FTYPE_RGGAS)
        SELECT CASE(X1AXIS)
        CASE(IAXIS)
        U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K)=U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K)+FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
        CASE(JAXIS)
        U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K)=U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K)+FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
        CASE(KAXIS)
        U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE)=U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE)+FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
        END SELECT
      CASE(CC_FTYPE_CFGAS) ! Cut-cell
        ICC = RCF%CELL_LIST(2,ISIDE+2)
        JCC = RCF%CELL_LIST(3,ISIDE+2)
        CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
                                                   FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
      END SELECT

   ELSE RCF_ON_WALL_CELL_IF
      RHO_Z_PV(-2:1) = 0._EB
      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            AF = DY(J)*DZ(K)
            RHOPV(-2:1)      = RHOP(I-1:I+2,J,K)
            ! First two cells surrounding face:
            DO ISIDE=-1,0
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = UU(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)

            DO ISIDE=-1,0
               FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell
                  U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K) = U_DOT_DEL_RHO_Z_VOL(I+1+ISIDE,J,K) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               CASE(CC_FTYPE_CFGAS) ! Cut-cell
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               END SELECT
            ENDDO
         CASE(JAXIS)
            AF = DX(I)*DZ(K)
            RHOPV(-2:1)      = RHOP(I,J-1:J+2,K)
            DO ISIDE=-1,0
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = VV(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)
            DO ISIDE=-1,0
               FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell
                  U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K) = U_DOT_DEL_RHO_Z_VOL(I,J+1+ISIDE,K) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               CASE(CC_FTYPE_CFGAS) ! Cut-cell
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               END SELECT
            ENDDO
         CASE(KAXIS)
            AF = DX(I)*DY(J)
            RHOPV(-2:1)      = RHOP(I,J,K-1:K+2)
            DO ISIDE=-1,0
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)* CUT_CELL(ICC)%RHOS(JCC)
                  ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDDO
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
            ENDIF
            ! Now Godunov flux limited value of rho*zz on face:
            VELC = WW(I,J,K)
            ! bar{rho*zz}:
            Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
            U_TEMP(1,1,1) = VELC
            CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
            FN_ZZ = F_TEMP(1,1,1)
            DO ISIDE=-1,0
               FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
               SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell
                  U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE) = U_DOT_DEL_RHO_Z_VOL(I,J,K+1+ISIDE) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               CASE(CC_FTYPE_CFGAS) ! Cut-cell
                  ICC = RCF%CELL_LIST(2,ISIDE+2)
                  IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
                  JCC = RCF%CELL_LIST(3,ISIDE+2)
                  CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
                  FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
               END SELECT
            ENDDO
      ENDSELECT
   ENDIF RCF_ON_WALL_CELL_IF
ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   CF =>  CUT_FACE(ICF);  IF ( CF%STATUS /= CC_GASPHASE ) CYCLE
   IW = CF%IWC; I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   CF_ON_WALL_CELL_IF : IF (IW > 0) THEN
      IF(WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
      WC => WALL(IW)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR
      FN_ZZ = B1%RHO_F*B1%ZZ_F(N)
      ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
      !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
      ISIDE = -1 + (SIGN(1,IOR)+1) / 2
      FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
      DO IFACE=1,CF%NFACE
         AF   = CF%AREA(IFACE)
         IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN ! Here used the CFA corresponding U_NORMAL.
            FCT  = 1._EB
            B1 => BOUNDARY_PROP1(CFACE(CF%CFACE_INDEX(IFACE))%B1_INDEX)
            VELC =         PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
         ELSEIF(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
            VELC = CUT_FACE(ICF)%VEL_SAVE(IFACE)
         ELSE
            ! Last known cut-face velocity.
            VELC = (1._EB-PRFCT)*CF%VEL(IFACE) + PRFCT*CF%VELS(IFACE)
         ENDIF
         ! First (rho hs)_i,j,k:
         IF (CF%CELL_LIST(1,ISIDE+2,IFACE) == CC_FTYPE_CFGAS) THEN
           ICC = CF%CELL_LIST(2,ISIDE+2,IFACE)
           JCC = CF%CELL_LIST(3,ISIDE+2,IFACE)
           RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
           ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
           RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
           CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
                                                      FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
         ENDIF
      ENDDO ! IFACE

   ELSE CF_ON_WALL_CELL_IF
      DO IFACE=1,CF%NFACE
         AF = CF%AREA(IFACE)
         RHOPV(-1:0)    = -1._EB
         RHO_Z_PV(-1:0) =  0._EB
         DO ISIDE=-1,0
            SELECT CASE(CF%CELL_LIST(1,ISIDE+2,IFACE))
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = CF%CELL_LIST(2,ISIDE+2,IFACE)
               JCC = CF%CELL_LIST(3,ISIDE+2,IFACE)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
         ENDDO
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I+1+ISIDE,J,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*ZZ_GET_N
            ENDIF
         CASE(JAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J+1+ISIDE,K,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*ZZ_GET_N
            ENDIF
         CASE(KAXIS)
            ! Lower cell:
            ISIDE=-2
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE+1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
            ENDIF
            ! Upper cell:
            ISIDE=1
            IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
               RHO_Z_PV(ISIDE) = RHO_Z_PV(ISIDE-1) ! Use center cell.
            ELSE
               ZZ_GET_N = ZZP(I,J,K+1+ISIDE,N)
               RHO_Z_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*ZZ_GET_N
            ENDIF
         END SELECT
         VELC  = PRFCTV *CF%VEL(IFACE) + (1._EB-PRFCTV)*CF%VELS(IFACE)
         ! bar{rho*zz}:
         Z_TEMP(0:3,1,1) = RHO_Z_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_ZZ = F_TEMP(1,1,1)
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(CF%CELL_LIST(1,ISIDE+2,IFACE))
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = CF%CELL_LIST(2,ISIDE+2,IFACE)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = CF%CELL_LIST(3,ISIDE+2,IFACE)
               CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) = CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + &
               FCT*SGNFCT*(FN_ZZ-RHO_Z_PV(ISIDE))*VELC * AF ! +ve or -ve dot
            END SELECT
         ENDDO
      ENDDO ! IFACE
   ENDIF CF_ON_WALL_CELL_IF
ENDDO ! ICF

! INBOUNDARY cut-faces:
! Species advection due to INBOUNDARY cut-faces (CFACE):
ISIDE=-1
CFACE_LOOP : DO ICFA=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA => CFACE(ICFA)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   ! Find associated cut-cell:
   IND1=CFA%CUT_FACE_IND1
   IND2=CFA%CUT_FACE_IND2
   ICC = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2)
   JCC = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
   AF  = CFA%AREA
   VELC= PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   ! Takes place of flux limited interpolation:
   FN_ZZ        = B1%RHO_F * B1%ZZ_F(N)
   RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT) *CUT_CELL(ICC)%RHOS(JCC)
   ZZ_GET_N     = PRFCT*CUT_CELL(ICC)%ZZ(N,JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
   ! Cut-cell value of rho*Z:
   RHO_Z_PV(ISIDE) = RHOPV(ISIDE)*ZZ_GET_N
   ! Add to U_DOT_DEL_RHO_Z:                                        ! (\bar{rho*Z}_CFACE - (rho*Z)_CC)*VELOUT*AF
   CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC)=CUT_CELL(ICC)%U_DOT_DEL_RHO_Z_VOL(N,JCC) + (FN_ZZ-RHO_Z_PV(ISIDE))*VELC*AF
ENDDO CFACE_LOOP

RETURN
END SUBROUTINE CCSPECIES_ADVECTION


! ---------------------------- CCENTHALPY_ADVECTION -----------------------------

SUBROUTINE CCENTHALPY_ADVECTION


! Computes FV version of flux limited \bar{ u dot Grad rho hs} in faces of near IB
! region and adds components to thermodynamic divergence.


! Local Variables:
REAL(EB) :: RHO_H_S_PV(-2:1), VELC, VELC2, FN_H_S, TMP_F_GAS
INTEGER  :: IOR, ICFA
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(CC_REGFACEZ_TYPE),  POINTER, DIMENSION(:) :: REGFACE_Z
LOGICAL :: DO_LO, DO_HI

U_TEMP => U_WORK
F_TEMP => F_WORK
Z_TEMP => Z_WORK

! IAXIS faces:
X1AXIS = IAXIS
REGFACE_Z => CC_REGFACE_IAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I    = REGFACE_Z(IFACE)%IJK(IAXIS)
   J    = REGFACE_Z(IFACE)%IJK(JAXIS)
   K    = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO= REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI= REGFACE_Z(IFACE)%DO_HI_IND
   RHOPV(-1:0)      = RHOP(I:I+1,J,K)
   TMPV(-1:0)       =  TMP(I:I+1,J,K)
   RHO_H_S_PV(-1:0) = 0._EB
   ! Get rho*hs on cells at both sides of IFACE:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
      RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ENDDO
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF  = DY(J)*DZ(K)
   VELC=UU(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
   IF(DO_LO) DPVOL(I  ,J,K) = DPVOL(I  ,J,K) + (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV(-1))*VELC*AF ! +ve dot
   IF(DO_HI) DPVOL(I+1,J,K) = DPVOL(I+1,J,K) - (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV( 0))*VELC*AF ! -ve dot
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
REGFACE_Z => CC_REGFACE_JAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I    = REGFACE_Z(IFACE)%IJK(IAXIS)
   J    = REGFACE_Z(IFACE)%IJK(JAXIS)
   K    = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO= REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI= REGFACE_Z(IFACE)%DO_HI_IND
   RHOPV(-1:0)      = RHOP(I,J:J+1,K)
   TMPV(-1:0)       =  TMP(I,J:J+1,K)
   RHO_H_S_PV(-1:0) = 0._EB
   ! Get rho*hs on cells at both sides of IFACE:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
      RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ENDDO
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF = DX(I)*DZ(K)
   VELC=VV(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
   IF(DO_LO) DPVOL(I,J  ,K) = DPVOL(I,J  ,K) + (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV(-1))*VELC*AF ! +ve dot
   IF(DO_HI) DPVOL(I,J+1,K) = DPVOL(I,J+1,K) - (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV( 0))*VELC*AF ! -ve dot
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
REGFACE_Z => CC_REGFACE_KAXIS_Z
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = REGFACE_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I    = REGFACE_Z(IFACE)%IJK(IAXIS)
   J    = REGFACE_Z(IFACE)%IJK(JAXIS)
   K    = REGFACE_Z(IFACE)%IJK(KAXIS)
   DO_LO= REGFACE_Z(IFACE)%DO_LO_IND
   DO_HI= REGFACE_Z(IFACE)%DO_HI_IND
   RHOPV(-1:0)      = RHOP(I,J,K:K+1)
   TMPV(-1:0)       =  TMP(I,J,K:K+1)
   RHO_H_S_PV(-1:0) = 0._EB
   ! Get rho*hs on cells at both sides of IFACE:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
      CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
      RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ENDDO
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF = DX(I)*DY(J)
   VELC=WW(I,J,K); IF(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) VELC=UVW_SAVE(IW)
   IF(DO_LO) DPVOL(I,J,K  ) = DPVOL(I,J,K  ) + (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV(-1))*VELC*AF ! +ve dot
   IF(DO_HI) DPVOL(I,J,K+1) = DPVOL(I,J,K+1) - (-1._EB)*(REGFACE_Z(IFACE)%FN_H_S-RHO_H_S_PV( 0))*VELC*AF ! -ve dot
ENDDO

! Regular faces connecting gasphase - cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z
   IW = RC_FACE(IFACE)%IWC
   IF( IW > 0 ) CYCLE
   I      = RC_FACE(IFACE)%IJK(IAXIS)
   J      = RC_FACE(IFACE)%IJK(JAXIS)
   K      = RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)
   RHO_H_S_PV(-2:1) = 0._EB
   SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF = DY(J)*DZ(K)
         RHOPV(-1:0)      = RHOP(I:I+1,J,K)
         TMPV(-1:0)       =  TMP(I:I+1,J,K)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                       (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
            RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
         ENDDO
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I+1+ISIDE,J,K))
            RHO_H_S_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I+1+ISIDE,J,K))
            RHO_H_S_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*H_S
         ENDIF
         ! Now Godunov flux limited value of rho*hs on face:
         VELC = UU(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_H_S_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_H_S = F_TEMP(1,1,1)
         ! Add contribution to DP:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I+1+ISIDE,J,K)=DPVOL(I+1+ISIDE,J,K)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF !+ve/-ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
            END SELECT
         ENDDO

      CASE(JAXIS)
         AF = DX(I)*DZ(K)
         RHOPV(-1:0)      = RHOP(I,J:J+1,K)
         TMPV(-1:0)       =  TMP(I,J:J+1,K)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                       (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
            RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
         ENDDO
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J+1+ISIDE,K))
            RHO_H_S_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J+1+ISIDE,K))
            RHO_H_S_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*H_S
         ENDIF
         ! Now Godunov flux limited value of rho*hs on face:
         VELC = VV(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_H_S_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_H_S = F_TEMP(1,1,1)
         ! Add contribution to DP:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I,J+1+ISIDE,K)=DPVOL(I,J+1+ISIDE,K)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF !+ve/-ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
            END SELECT
         ENDDO

      CASE(KAXIS)
         AF = DX(I)*DY(J)
         RHOPV(-1:0)      = RHOP(I,J,K:K+1)
         TMPV(-1:0)       =  TMP(I,J,K:K+1)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                       (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
            RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
         ENDDO
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K+1+ISIDE))
            RHO_H_S_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K+1+ISIDE))
            RHO_H_S_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*H_S
         ENDIF
         ! Now Godunov flux limited value of rho*hs on face:
         VELC = WW(I,J,K)
         Z_TEMP(0:3,1,1) = RHO_H_S_PV(-2:1)
         U_TEMP(1,1,1) = VELC
         CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
         FN_H_S = F_TEMP(1,1,1)
         ! Add contribution to DP:
         ! Low side cell:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I,J,K+1+ISIDE)=DPVOL(I,J,K+1+ISIDE)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF !+ve/-ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
            END SELECT
         ENDDO

   ENDSELECT

ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE
   IW = CUT_FACE(ICF)%IWC
   IF(IW > 0) CYCLE
   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
   DO IFACE=1,CUT_FACE(ICF)%NFACE
      AF = CUT_FACE(ICF)%AREA(IFACE)
      RHOPV(-1:0)      = -1._EB
      TMPV(-1:0)       = -1._EB
      RHO_H_S_PV(-2:1) =  0._EB
      DO ISIDE=-1,0
         ZZ_GET = 0._EB
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
            RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
            ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                    (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
         END SELECT
         CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
         RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
      ENDDO
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I+1+ISIDE,J,K))
            RHO_H_S_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I+1+ISIDE,J,K))%SOLID .OR. CCVAR(I+1+ISIDE,J,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I+1+ISIDE,J,K))
            RHO_H_S_PV(ISIDE) = RHOP(I+1+ISIDE,J,K)*H_S
         ENDIF
      CASE(JAXIS)
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J+1+ISIDE,K))
            RHO_H_S_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J+1+ISIDE,K))%SOLID .OR. CCVAR(I,J+1+ISIDE,K,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J+1+ISIDE,K))
            RHO_H_S_PV(ISIDE) = RHOP(I,J+1+ISIDE,K)*H_S
         ENDIF
      CASE(KAXIS)
         ! Lower cell:
         ISIDE=-2
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE+1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K+1+ISIDE))
            RHO_H_S_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*H_S
         ENDIF
         ! Upper cell:
         ISIDE=1
         IF (CELL(CELL_INDEX(I,J,K+1+ISIDE))%SOLID .OR. CCVAR(I,J,K+1+ISIDE,CC_CGSC)==CC_SOLID) THEN
            RHO_H_S_PV(ISIDE) = RHO_H_S_PV(ISIDE-1) ! Use center cell.
         ELSE
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP(I,J,K+1+ISIDE))
            RHO_H_S_PV(ISIDE) = RHOP(I,J,K+1+ISIDE)*H_S
         ENDIF
      END SELECT
      VELC    = PRFCTV *CUT_FACE(ICF)%VEL(IFACE) + (1._EB-PRFCTV)*CUT_FACE(ICF)%VELS(IFACE)
      Z_TEMP(0:3,1,1) = RHO_H_S_PV(-2:1)
      U_TEMP(1,1,1) = VELC
      CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,1,I_FLUX_LIMITER)
      FN_H_S = F_TEMP(1,1,1)
      ! Add to divergence integral of surrounding cut-cells:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
         END SELECT
      ENDDO
   ENDDO ! IFACE
ENDDO ! ICF

! Now work with boundary faces:
! IAXIS faces:
X1AXIS = IAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(KAXIS)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   RHOPV(ISIDE)      = RHOP(I+1+ISIDE,J,K)
   TMPV(ISIDE)       =  TMP(I+1+ISIDE,J,K)
   ! Get rho*hs on cells at both sides of IFACE:
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
   RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.
   VELC      = UU(I,J,K)
   VELC2     = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_F_GAS = B1%TMP_F
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      VELC = -SIGN(1._EB,REAL(IOR,EB))*VELC2
      IF (VELC2>0._EB) TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
   ENDIF
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF = DY(J)*DZ(K)
   DPVOL(I+1+ISIDE,J,K) = DPVOL(I+1+ISIDE,J,K) + SIGN(1._EB,REAL(IOR,EB))*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(KAXIS)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   RHOPV(ISIDE)      = RHOP(I,J+1+ISIDE,K)
   TMPV(ISIDE)       =  TMP(I,J+1+ISIDE,K)
   ! Get rho*hs on cells at both sides of IFACE:
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
   RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.
   VELC      = VV(I,J,K)
   VELC2     = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_F_GAS = B1%TMP_F
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      VELC = -SIGN(1._EB,REAL(IOR,EB))*VELC2
      IF (VELC2>0._EB) TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
   ENDIF
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF = DX(I)*DZ(K)
   DPVOL(I,J+1+ISIDE,K) = DPVOL(I,J+1+ISIDE,K) + SIGN(1._EB,REAL(IOR,EB))*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   I  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(KAXIS)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   RHOPV(ISIDE)      = RHOP(I,J,K+1+ISIDE)
   TMPV(ISIDE)       =  TMP(I,J,K+1+ISIDE)
   ! Get rho*hs on cells at both sides of IFACE:
   ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
   RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.
   VELC      = WW(I,J,K)
   VELC2     = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_F_GAS = B1%TMP_F
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN
      VELC = -SIGN(1._EB,REAL(IOR,EB))*VELC2
      IF (VELC2>0._EB) TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
   ENDIF
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}
   ! Add: -(bar{rho*hs} u dot n - (rho*hs) u dot n) to corresponding cell DP:
   AF = DX(I)*DY(J)
   DPVOL(I,J,K+1+ISIDE) = DPVOL(I,J,K+1+ISIDE) + SIGN(1._EB,REAL(IOR,EB))*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
ENDDO

! Regular Faces connecting gasphase cells to cut-cells:
DO IFACE=1,MESHES(NM)%CC_NBBRCFACE_Z
   IW = RC_FACE(IFACE)%IWC
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE
   I      = RC_FACE(IFACE)%IJK(IAXIS)
   J      = RC_FACE(IFACE)%IJK(JAXIS)
   K      = RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   FCT   = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
   ! First (rho hs)_i,j,k:
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      AF = DY(J)*DZ(K)
      VELC = UU(I,J,K)
      RHOPV(ISIDE)      = RHOP(I+1+ISIDE,J,K)
      TMPV(ISIDE)       =  TMP(I+1+ISIDE,J,K)
      SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
      CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
      CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
         RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
         ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
      END SELECT
   CASE(JAXIS)
      AF = DX(I)*DZ(K)
      VELC = VV(I,J,K)
      RHOPV(ISIDE)      = RHOP(I,J+1+ISIDE,K)
      TMPV(ISIDE)       =  TMP(I,J+1+ISIDE,K)
      SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
      CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
      CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
         RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
         ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
      END SELECT
   CASE(KAXIS)
      AF = DX(I)*DY(J)
      VELC = WW(I,J,K)
      RHOPV(ISIDE)      = RHOP(I,J,K+1+ISIDE)
      TMPV(ISIDE)       =  TMP(I,J,K+1+ISIDE)
      SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
      CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
         ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
      CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
         RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
         ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
      END SELECT
   END SELECT
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
   RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S

   ! Flux limited face value bar{rho*hs}_F
   ! Calculate the sensible enthalpy at the boundary. If the boundary is solid
   ! and the gas is flowing out, use the gas temperature for the calculation.
   VELC2     = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_F_GAS = B1%TMP_F
   SELECT CASE(WC%BOUNDARY_TYPE)
      CASE DEFAULT
         ! No need to do anything, populated before.
      CASE(SOLID_BOUNDARY)
         VELC = -SIGN(1._EB,REAL(IOR,EB))*VELC2
         IF (VELC2>0._EB) TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
      CASE(INTERPOLATED_BOUNDARY)
         VELC = UVW_SAVE(IW)
   END SELECT
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}
   ! Finally add to Div:
   SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
   CASE(CC_FTYPE_RGGAS) ! Regular cell
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         DPVOL(I+1+ISIDE,J,K)=DPVOL(I+1+ISIDE,J,K)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF !+ve/-ve dot
      CASE(JAXIS)
         DPVOL(I,J+1+ISIDE,K)=DPVOL(I,J+1+ISIDE,K)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
      CASE(KAXIS)
         DPVOL(I,J,K+1+ISIDE)=DPVOL(I,J,K+1+ISIDE)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
      END SELECT
   CASE(CC_FTYPE_CFGAS) ! Cut-cell
      ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
      JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
      CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF
   END SELECT
ENDDO

! Finally Gasphase cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE .OR. MESHES(NM)%CUT_FACE(ICF)%IWC<1) CYCLE
   IW = CUT_FACE(ICF)%IWC
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY ) CYCLE
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IOR = BC%IOR
   ! Flux limited face value bar{rho*hs}_F, the P1 variable values fo TMP, RHOP, ZZ and RSUM have been averaged to
   ! the cartesian cell location in CC_DENSITY:
   VELC2     = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_F_GAS = B1%TMP_F
   IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. VELC2>0._EB) TMP_F_GAS = TMP(BC%IIG,BC%JJG,BC%KKG)
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
   DO IFACE=1,CUT_FACE(ICF)%NFACE
      AF   = CUT_FACE(ICF)%AREA(IFACE)
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN ! Here used the CFA corresponding U_NORMAL.
         FCT  = 1._EB
         B1 => BOUNDARY_PROP1(CFACE(CUT_FACE(ICF)%CFACE_INDEX(IFACE))%B1_INDEX)
         VELC =         PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
      ELSEIF(WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) THEN
         VELC = CUT_FACE(ICF)%VEL_SAVE(IFACE)
      ELSE
         ! Last known cut-face velocity.
         VELC = (1._EB-PRFCT)*CUT_FACE(ICF)%VEL(IFACE) + PRFCT*CUT_FACE(ICF)%VELS(IFACE)
      ENDIF
      ! Here if INTERPOLATED_BOUNDARY we might need UVW_SAVE.
      ! First (rho hs)_i,j,k:
      SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
      CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
         ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
         JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
         TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
         RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
         ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
         CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
         RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
         CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF ! +ve or -ve dot
      END SELECT
   ENDDO ! IFACE
ENDDO ! ICF

! Enthalpy advection due to INBOUNDARY cut-faces (CFACE):
ISIDE=-1
CFACE_LOOP : DO ICFA=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA => CFACE(ICFA)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   ! Find associated cut-cell:
   IND1=CFA%CUT_FACE_IND1
   IND2=CFA%CUT_FACE_IND2
   ICC = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2)
   JCC = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
   AF = CFA%AREA
   VELC      = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S ! Contains AREA_ADJUST for the CFACE.
   TMP_F_GAS = B1%TMP_F
   IF (VELC>0._EB) TMP_F_GAS = B1%TMP_G ! CUT_CELL(ICC)%TMP(JCC)
   ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMP_F_GAS)
   FN_H_S = B1%RHO_F*H_S ! bar{rho*hs}

   TMPV(ISIDE)  = CUT_CELL(ICC)%TMP(JCC)
   RHOPV(ISIDE) = PRFCT *CUT_CELL(ICC)%RHO(JCC) + (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
   ZZ_GET(1:N_TRACKED_SPECIES) = PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                          (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
   CALL GET_SENSIBLE_ENTHALPY(ZZ_GET,H_S,TMPV(ISIDE))
   RHO_H_S_PV(ISIDE) = RHOPV(ISIDE)*H_S
   CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+(-1._EB)*(FN_H_S-RHO_H_S_PV(ISIDE))*VELC*AF ! +ve or -ve dot
ENDDO CFACE_LOOP

RETURN
END SUBROUTINE CCENTHALPY_ADVECTION

! ----------------------- CC_DIFFUSIVE_HEAT_FLUXES ------------------------

SUBROUTINE CC_DIFFUSIVE_HEAT_FLUXES

! NOTE: this routine assumes POINT_TO_MESH(NM) has been previously called.

! Local Variables:
INTEGER :: IIG, JJG, KKG , IOR
REAL(EB) :: UN_P
LOGICAL :: DO_LO, DO_HI

SPECIES_LOOP1: DO N=1,N_TOTAL_SCALARS

   ! IAXIS faces:
   X1AXIS = IAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

      I     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(IAXIS)
      J     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(JAXIS)
      K     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(KAXIS)
      DO_LO = CC_REGFACE_IAXIS_Z(IFACE)%DO_LO_IND
      DO_HI = CC_REGFACE_IAXIS_Z(IFACE)%DO_HI_IND

      ! H_RHO_D_DZDN
      TMP_G = 0.5_EB*(TMP(I+1,J,K)+TMP(I,J,K))
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      AF = DY(J)*DZ(K)
      IF (DO_LO) THEN
         DPVOL(I  ,J,K) = DPVOL(I  ,J,K) + CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
         DEL_RHO_D_DEL_Z(I  ,J,K,N)=DEL_RHO_D_DEL_Z(I  ,J,K,N)+CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !+ dot
      ENDIF
      IF (DO_HI) THEN
         DPVOL(I+1,J,K) = DPVOL(I+1,J,K) - CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
         DEL_RHO_D_DEL_Z(I+1,J,K,N)=DEL_RHO_D_DEL_Z(I+1,J,K,N)-CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !- dot
      ENDIF
   ENDDO

   ! JAXIS faces:
   X1AXIS = JAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

      I     = CC_REGFACE_JAXIS_Z(IFACE)%IJK(IAXIS)
      J     = CC_REGFACE_JAXIS_Z(IFACE)%IJK(JAXIS)
      K     = CC_REGFACE_JAXIS_Z(IFACE)%IJK(KAXIS)
      DO_LO = CC_REGFACE_JAXIS_Z(IFACE)%DO_LO_IND
      DO_HI = CC_REGFACE_JAXIS_Z(IFACE)%DO_HI_IND

      ! H_RHO_D_DZDN
      TMP_G = 0.5_EB*(TMP(I,J+1,K)+TMP(I,J,K))
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      AF = DX(I)*DZ(K)
      IF (DO_LO) THEN
         DPVOL(I,J  ,K) = DPVOL(I,J  ,K) + CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
         DEL_RHO_D_DEL_Z(I,J  ,K,N)=DEL_RHO_D_DEL_Z(I,J  ,K,N)+CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !+ dot
      ENDIF
      IF (DO_HI) THEN
         DPVOL(I,J+1,K) = DPVOL(I,J+1,K) - CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
         DEL_RHO_D_DEL_Z(I,J+1,K,N)=DEL_RHO_D_DEL_Z(I,J+1,K,N)-CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !- dot
      ENDIF
   ENDDO

   ! KAXIS faces:
   X1AXIS = KAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
      I     = CC_REGFACE_KAXIS_Z(IFACE)%IJK(IAXIS)
      J     = CC_REGFACE_KAXIS_Z(IFACE)%IJK(JAXIS)
      K     = CC_REGFACE_KAXIS_Z(IFACE)%IJK(KAXIS)
      DO_LO = CC_REGFACE_KAXIS_Z(IFACE)%DO_LO_IND
      DO_HI = CC_REGFACE_KAXIS_Z(IFACE)%DO_HI_IND

      ! H_RHO_D_DZDN
      TMP_G = 0.5_EB*(TMP(I,J,K+1)+TMP(I,J,K))
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      AF = DX(I)*DY(J)
      IF (DO_LO) THEN
         DPVOL(I,J,K  ) = DPVOL(I,J,K  ) + CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
         DEL_RHO_D_DEL_Z(I,J,K  ,N)=DEL_RHO_D_DEL_Z(I,J,K  ,N)+CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !+ dot
      ENDIF
      IF (DO_HI) THEN
         DPVOL(I,J,K+1) = DPVOL(I,J,K+1) - CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
         DEL_RHO_D_DEL_Z(I,J,K+1,N)=DEL_RHO_D_DEL_Z(I,J,K+1,N)-CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)*AF !- dot
      ENDIF
   ENDDO

ENDDO SPECIES_LOOP1

! Regular faces connecting gasphase- cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z
   IW = RC_FACE(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I      = RC_FACE(IFACE)%IJK(IAXIS)
   J      = RC_FACE(IFACE)%IJK(JAXIS)
   K      = RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)
   TMP_G  = RC_FACE(IFACE)%TMP_FACE
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      AF = DY(J)*DZ(K)
      ! H_RHO_D_DZDN
      DO N=1,N_TOTAL_SCALARS
         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
         RC_FACE(IFACE)%H_RHO_D_DZDN(N) = H_S*RC_FACE(IFACE)%RHO_D_DZDN(N)
      ENDDO
      ! Add contribution to DP:
      ! Low side cell:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell
         DPVOL(I+1+ISIDE,J,K)=DPVOL(I+1+ISIDE,J,K)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF
         ! +ve or -ve dot
         DO N=1,N_TOTAL_SCALARS
         DEL_RHO_D_DEL_Z(I+1+ISIDE,J,K,N)=DEL_RHO_D_DEL_Z(I+1+ISIDE,J,K,N)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(N)*AF
         ENDDO
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF !+/- dot
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)= &
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(1:N_TOTAL_SCALARS)*AF
         END SELECT
      ENDDO
   CASE(JAXIS)
      AF = DX(I)*DZ(K)
      ! H_RHO_D_DZDN
      DO N=1,N_TOTAL_SCALARS
         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
         RC_FACE(IFACE)%H_RHO_D_DZDN(N) = H_S*RC_FACE(IFACE)%RHO_D_DZDN(N)
      ENDDO
      ! Add contribution to DP:
      ! Low side cell:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell
         DPVOL(I,J+1+ISIDE,K)=DPVOL(I,J+1+ISIDE,K)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF
         ! +ve or -ve dot
         DO N=1,N_TOTAL_SCALARS
         DEL_RHO_D_DEL_Z(I,J+1+ISIDE,K,N)=DEL_RHO_D_DEL_Z(I,J+1+ISIDE,K,N)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(N)*AF
         ENDDO
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF !+/- dot
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)= &
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(1:N_TOTAL_SCALARS)*AF
         END SELECT
      ENDDO
   CASE(KAXIS)
      AF = DX(I)*DY(J)
      ! H_RHO_D_DZDN
      DO N=1,N_TOTAL_SCALARS
         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
         RC_FACE(IFACE)%H_RHO_D_DZDN(N) = H_S*RC_FACE(IFACE)%RHO_D_DZDN(N)
      ENDDO
      ! Add contribution to DP:
      ! Low side cell:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
         CASE(CC_FTYPE_RGGAS) ! Regular cell
         DPVOL(I,J,K+1+ISIDE)=DPVOL(I,J,K+1+ISIDE)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF
         ! +ve or -ve dot
         DO N=1,N_TOTAL_SCALARS
         DEL_RHO_D_DEL_Z(I,J,K+1+ISIDE,N)=DEL_RHO_D_DEL_Z(I,J,K+1+ISIDE,N)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(N)*AF
         ENDDO
         CASE(CC_FTYPE_CFGAS) ! Cut-cell
         ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
         IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE
         JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
         CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)+FCT*SUM(RC_FACE(IFACE)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF !+/- dot
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)= &
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)+FCT*RC_FACE(IFACE)%RHO_D_DZDN(1:N_TOTAL_SCALARS)*AF
         END SELECT
      ENDDO
   END SELECT
ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE
   IW = CUT_FACE(ICF)%IWC
   ! Note: for cut-faces open boundaries are dealt with below in external BC loops:
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
   DO IFACE=1,CUT_FACE(ICF)%NFACE
      AF = CUT_FACE(ICF)%AREA(IFACE)
      ! H_RHO_D_DZDN
      TMP_G = CUT_FACE(ICF)%TMP_FACE(IFACE)
      DO N=1,N_TOTAL_SCALARS
         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
         CUT_FACE(ICF)%H_RHO_D_DZDN(N,IFACE) = H_S*CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE)
      ENDDO
      ! Add to divergence integral of surrounding cut-cells:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
         ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
         IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE
         JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
         CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC)+FCT*SUM(CUT_FACE(ICF)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS,IFACE))*AF !+/- dot
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)= &
         CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC)+FCT*CUT_FACE(ICF)%RHO_D_DZDN(1:N_TOTAL_SCALARS,IFACE)*AF
         END SELECT
      ENDDO
   ENDDO ! IFACE
ENDDO ! ICF


! Now define diffussive heat flux components in Boundaries:
! CFACES:
ISIDE=-1
CFACE_LOOP : DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA => CFACE(ICF)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   IND1 = CFA%CUT_FACE_IND1;                         IND2 = CFA%CUT_FACE_IND2
   ICC  = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2);  JCC  = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
   ! H_RHO_D_DZDN
   UN_P =  PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
   TMP_G = B1%TMP_F
   IF (UN_P>0._EB) TMP_G = B1%TMP_G
   DO N=1,N_TOTAL_SCALARS
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CUT_FACE(IND1)%H_RHO_D_DZDN(N,IND2) = H_S*B1%RHO_D_DZDN_F(N)
   ENDDO
   AF = CFA%AREA  ! No need for B1%AREA_ADJUST, RHO_D_DZDN_F is already area djusted. Same for Domain Boundaries below.
   ! Add diffusive mass flux enthalpy contribution to cut-cell thermo divg:
   CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) - SUM(CUT_FACE(IND1)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS,IND2))*AF
   CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) = &
   CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) - B1%RHO_D_DZDN_F(1:N_TOTAL_SCALARS)*AF
ENDDO CFACE_LOOP

! Domain boundaries:
SPECIES_LOOP2: DO N=1,N_TOTAL_SCALARS

   ! IAXIS faces:
   X1AXIS = IAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF( ANY(WC%BOUNDARY_TYPE==(/NULL_BOUNDARY,INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/)) ) CYCLE
      I  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(IAXIS)
      J  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(JAXIS)
      K  = CC_REGFACE_IAXIS_Z(IFACE)%IJK(KAXIS)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      UN_P = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
      TMP_G = B1%TMP_F
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) TMP_G = TMP(BC%IIG,BC%JJG,BC%KKG)
      ! H_RHO_D_DZDN
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP:
      AF = DY(J)*DZ(K)
      SELECT CASE(BOUNDARY_COORD(WC%BC_INDEX)%IOR)
      CASE(-IAXIS) ! Low side cell. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      IF (.NOT.CC_REGFACE_IAXIS_Z(IFACE)%DO_LO_IND) CYCLE
      DPVOL(I  ,J,K)             =             DPVOL(I  ,J,K) + CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
      DEL_RHO_D_DEL_Z(I  ,J,K,N) = DEL_RHO_D_DEL_Z(I  ,J,K,N) + CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! +ve dot
      CASE( IAXIS) ! High side cell.
      IF (.NOT.CC_REGFACE_IAXIS_Z(IFACE)%DO_HI_IND) CYCLE
      DPVOL(I+1,J,K)             =             DPVOL(I+1,J,K) - CC_REGFACE_IAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
      DEL_RHO_D_DEL_Z(I+1,J,K,N) = DEL_RHO_D_DEL_Z(I+1,J,K,N) - CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! -ve dot
      END SELECT
   ENDDO

   ! JAXIS faces:
   X1AXIS = JAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF( ANY(WC%BOUNDARY_TYPE==(/NULL_BOUNDARY,INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/)) ) CYCLE
      I  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(IAXIS)
      J  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(JAXIS)
      K  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(KAXIS)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      UN_P = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
      TMP_G = B1%TMP_F
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) TMP_G = TMP(BC%IIG,BC%JJG,BC%KKG)
      ! H_RHO_D_DZDN
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP:
      AF = DX(I)*DZ(K)
      SELECT CASE(BOUNDARY_COORD(WC%BC_INDEX)%IOR)
      CASE(-JAXIS) ! Low side cell. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      IF (.NOT.CC_REGFACE_JAXIS_Z(IFACE)%DO_LO_IND) CYCLE
      DPVOL(I,J  ,K)             =             DPVOL(I,J  ,K) + CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
      DEL_RHO_D_DEL_Z(I,J  ,K,N) = DEL_RHO_D_DEL_Z(I,J  ,K,N) + CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! +ve dot
      CASE( JAXIS) ! High side cell.
      IF (.NOT.CC_REGFACE_JAXIS_Z(IFACE)%DO_HI_IND) CYCLE
      DPVOL(I,J+1,K)             =             DPVOL(I,J+1,K) - CC_REGFACE_JAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
      DEL_RHO_D_DEL_Z(I,J+1,K,N) = DEL_RHO_D_DEL_Z(I,J+1,K,N) - CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! -ve dot
      END SELECT
   ENDDO

   ! KAXIS faces:
   X1AXIS = KAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE; WC => WALL(IW)
      IF( ANY(WC%BOUNDARY_TYPE==(/NULL_BOUNDARY,INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/)) ) CYCLE
      I  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(IAXIS)
      J  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(JAXIS)
      K  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(KAXIS)
      B1 => BOUNDARY_PROP1(WC%B1_INDEX)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      UN_P = PRFCT*B1%U_NORMAL + (1._EB-PRFCT)*B1%U_NORMAL_S
      TMP_G = B1%TMP_F
      IF (WC%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. UN_P>0._EB) TMP_G = TMP(BC%IIG,BC%JJG,BC%KKG)
      ! H_RHO_D_DZDN
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) = H_S*CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)

      ! Add H_RHO_D_DZDN dot n to corresponding cell DP:
      AF = DX(I)*DY(J)
      SELECT CASE(BOUNDARY_COORD(WC%BC_INDEX)%IOR)
      CASE(-KAXIS) ! Low side cell. Add to int(DEL_RHO_D_DEL_Z)dv in FV form:
      IF (.NOT.CC_REGFACE_KAXIS_Z(IFACE)%DO_LO_IND) CYCLE
      DPVOL(I,J,K  )             =             DPVOL(I,J,K  ) + CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! +ve dot
      DEL_RHO_D_DEL_Z(I,J,K  ,N) = DEL_RHO_D_DEL_Z(I,J,K  ,N) + CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! +ve dot
      CASE( KAXIS) ! High side cell.
      IF (.NOT.CC_REGFACE_KAXIS_Z(IFACE)%DO_HI_IND) CYCLE
      DPVOL(I,J,K+1)             =             DPVOL(I,J,K+1) - CC_REGFACE_KAXIS_Z(IFACE)%H_RHO_D_DZDN(N) * AF ! -ve dot
      DEL_RHO_D_DEL_Z(I,J,K+1,N) = DEL_RHO_D_DEL_Z(I,J,K+1,N) - CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N)   * AF ! -ve dot
      END SELECT
   ENDDO

ENDDO SPECIES_LOOP2

! Regular faces connecting gasphase- cut-cells:
DO IFACE=1,MESHES(NM)%CC_NBBRCFACE_Z
   RCF => RC_FACE(IFACE); IW = RCF%IWC; WC => WALL(IW)
   IF( ANY(WC%BOUNDARY_TYPE==(/NULL_BOUNDARY,INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/)) ) CYCLE
   I   = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1); TMP_G = RCF%TMP_FACE
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG; JJG = BC%JJG; KKG = BC%KKG; IOR = BC%IOR
   ! H_RHO_D_DZDN
   DO N=1,N_TOTAL_SCALARS
      CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
      RCF%H_RHO_D_DZDN(N) = H_S*RCF%RHO_D_DZDN(N)
   ENDDO
   SELECT CASE(X1AXIS)
   CASE(IAXIS); AF = DY(J)*DZ(K)
   CASE(JAXIS); AF = DX(I)*DZ(K)
   CASE(KAXIS); AF = DX(I)*DY(J)
   END SELECT
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
   SELECT CASE(RCF%CELL_LIST(1,ISIDE+2))
   CASE(CC_FTYPE_RGGAS) ! Regular cell
   DPVOL(IIG,JJG,KKG)=DPVOL(IIG,JJG,KKG)+FCT*SUM(RCF%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF ! +ve or -ve dot
   DO N=1,N_TOTAL_SCALARS
      DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N)=DEL_RHO_D_DEL_Z(IIG,JJG,KKG,N)+FCT*RCF%RHO_D_DZDN(N)*AF
   ENDDO
   CASE(CC_FTYPE_CFGAS) ! Cut-cell
   ICC = RCF%CELL_LIST(2,ISIDE+2); JCC = RCF%CELL_LIST(3,ISIDE+2)
   CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*SUM(RCF%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))*AF !+/- dot
   CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) = &
   CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) + FCT*RCF%RHO_D_DZDN(1:N_TOTAL_SCALARS) * AF
   END SELECT
ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   CF => CUT_FACE(ICF); IF ( CF%STATUS /= CC_GASPHASE .OR. CF%IWC<1) CYCLE
   IW = CF%IWC; WC => WALL(IW)
   IF( ANY(WC%BOUNDARY_TYPE==(/NULL_BOUNDARY,INTERPOLATED_BOUNDARY,PERIODIC_BOUNDARY/)) ) CYCLE
   I  = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   IOR= BOUNDARY_COORD(WC%BC_INDEX)%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   DO IFACE=1,CF%NFACE
      AF = CF%AREA(IFACE)
      ! H_RHO_D_DZDN
      TMP_G = CF%TMP_FACE(IFACE)
      DO N=1,N_TOTAL_SCALARS
         CALL GET_SENSIBLE_ENTHALPY_Z(N,TMP_G,H_S)
         CF%H_RHO_D_DZDN(N,IFACE) = H_S*CF%RHO_D_DZDN(N,IFACE)
      ENDDO
      ! Add to divergence integral of surrounding cut-cells:
      FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
      SELECT CASE(CF%CELL_LIST(1,ISIDE+2,IFACE))
      CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
      ICC = CF%CELL_LIST(2,ISIDE+2,IFACE); JCC = CF%CELL_LIST(3,ISIDE+2,IFACE)
      CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*SUM(CF%H_RHO_D_DZDN(1:N_TOTAL_SCALARS,IFACE)) * AF !+/- dot
      CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) = &
      CUT_CELL(ICC)%DEL_RHO_D_DEL_Z_VOL(1:N_TOTAL_SCALARS,JCC) + FCT*CF%RHO_D_DZDN(1:N_TOTAL_SCALARS,IFACE)*AF
      END SELECT
   ENDDO ! IFACE
ENDDO ! ICF

RETURN
END SUBROUTINE CC_DIFFUSIVE_HEAT_FLUXES

! ----------------------- CC_CONDUCTION_HEAT_FLUX --------------------------

SUBROUTINE CC_CONDUCTION_HEAT_FLUX

INTEGER :: IIG, JJG, KKG, IOR
REAL(EB):: KPDTDN=0._EB,KPV(-1:0)=0._EB

! IAXIS faces:
X1AXIS = IAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)

   IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   I     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(IAXIS)
   J     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(JAXIS)
   K     = CC_REGFACE_IAXIS_Z(IFACE)%IJK(KAXIS)

   ! K*DTDN:
   TMPV(-1:0)  = TMP(I:I+1,J,K)
   ! KP on low-high side cells:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
      CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I+1+ISIDE,J,K),&
                                             MU_DNS(I+1+ISIDE,J,K),TMPV(ISIDE),KPV(ISIDE))
   ENDDO
   KPDTDN = 0.5_EB*(KPV(-1)+KPV(0)) * (TMPV(0)-TMPV(-1)) / DX(I)

   ! Add K*DTDN dot n to corresponding cell DP:
   AF = DY(J)*DZ(K)
   IF(CC_REGFACE_IAXIS_Z(IFACE)%DO_LO_IND) DPVOL(I  ,J,K) = DPVOL(I  ,J,K) + KPDTDN * AF ! +ve dot
   IF(CC_REGFACE_IAXIS_Z(IFACE)%DO_HI_IND) DPVOL(I+1,J,K) = DPVOL(I+1,J,K) - KPDTDN * AF ! -ve dot
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)

   IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   I  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_JAXIS_Z(IFACE)%IJK(KAXIS)

   ! K*DTDN:
   TMPV(-1:0)  = TMP(I,J:J+1,K)
   ! KP on low-high side cells:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
      CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I,J+1+ISIDE,K),&
                                             MU_DNS(I,J+1+ISIDE,K),TMPV(ISIDE),KPV(ISIDE))
   ENDDO
   KPDTDN = 0.5_EB*(KPV(-1)+KPV(0)) * (TMPV(0)-TMPV(-1)) / DY(J)

   ! Add K*DTDN dot n to corresponding cell DP:
   AF = DX(I)*DZ(K)
   IF(CC_REGFACE_JAXIS_Z(IFACE)%DO_LO_IND) DPVOL(I,J  ,K) = DPVOL(I,J  ,K) + KPDTDN * AF ! +ve dot
   IF(CC_REGFACE_JAXIS_Z(IFACE)%DO_HI_IND) DPVOL(I,J+1,K) = DPVOL(I,J+1,K) - KPDTDN * AF ! -ve dot
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)

   IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   I  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(IAXIS)
   J  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(JAXIS)
   K  = CC_REGFACE_KAXIS_Z(IFACE)%IJK(KAXIS)

   ! K*DTDN:
   TMPV(-1:0)  = TMP(I,J,K:K+1)
   ! KP on low-high side cells:
   DO ISIDE=-1,0
      ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
      CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I,J,K+1+ISIDE),&
                                             MU_DNS(I,J,K+1+ISIDE),TMPV(ISIDE),KPV(ISIDE))
   ENDDO
   KPDTDN = 0.5_EB*(KPV(-1)+KPV(0)) * (TMPV(0)-TMPV(-1)) / DZ(K)

   ! Add K*DTDN dot n to corresponding cell DP:
   AF = DX(I)*DY(J)
   IF(CC_REGFACE_KAXIS_Z(IFACE)%DO_LO_IND) DPVOL(I,J,K  ) = DPVOL(I,J,K  ) + KPDTDN * AF ! +ve dot
   IF(CC_REGFACE_KAXIS_Z(IFACE)%DO_HI_IND) DPVOL(I,J,K+1) = DPVOL(I,J,K+1) - KPDTDN * AF ! -ve dot
ENDDO


! Regular faces connecting gasphase - cut-cells:
DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z

   IW = RC_FACE(IFACE)%IWC
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

   I      = RC_FACE(IFACE)%IJK(IAXIS)
   J      = RC_FACE(IFACE)%IJK(JAXIS)
   K      = RC_FACE(IFACE)%IJK(KAXIS)
   X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)

   SELECT CASE(X1AXIS)
      CASE(IAXIS)
         AF = DY(J)*DZ(K)
         X1F= MESHES(NM)%X(I)
         IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                         RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
         ! Linear interpolation coefficients:
         CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
         CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

         TMPV(-1:0)  = TMP(I:I+1,J,K)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I+1+ISIDE,J,K,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  &
                      PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
               (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            ! KP on low-high side cells:
            CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I+1+ISIDE,J,K),&
                                                   MU_DNS(I+1+ISIDE,J,K),TMPV(ISIDE),KPV(ISIDE))
         ENDDO

         KPDTDN = (CCM1*KPV(-1)+CCP1*KPV(0)) * (TMPV(0)-TMPV(-1)) * IDX

         ! Add contribution to DP:
         ! Low side cell:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I+1+ISIDE,J,K) = DPVOL(I+1+ISIDE,J,K) + FCT*KPDTDN * AF ! +ve or -ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*KPDTDN * AF ! +ve or -ve dot
            END SELECT
         ENDDO

      CASE(JAXIS)
         AF = DX(I)*DZ(K)
         X1F= MESHES(NM)%Y(J)
         IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                         RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
         ! Linear interpolation coefficients:
         CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
         CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

         TMPV(-1:0)  = TMP(I,J:J+1,K)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J+1+ISIDE,K,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  &
                      PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
               (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            ! KP on low-high side cells:
            CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I,J+1+ISIDE,K),&
                                                   MU_DNS(I,J+1+ISIDE,K),TMPV(ISIDE),KPV(ISIDE))
         ENDDO

         KPDTDN = (CCM1*KPV(-1)+CCP1*KPV(0)) * (TMPV(0)-TMPV(-1)) * IDX

         ! Add contribution to DP:
         ! Low side cell:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I,J+1+ISIDE,K) = DPVOL(I,J+1+ISIDE,K) + FCT*KPDTDN * AF ! +ve or -ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*KPDTDN * AF ! +ve or -ve dot
            END SELECT
         ENDDO

      CASE(KAXIS)
         AF = DX(I)*DY(J)
         X1F= MESHES(NM)%Z(K)
         IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                         RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
         ! Linear interpolation coefficients:
         CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
         CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

         TMPV(-1:0)  = TMP(I,J,K:K+1)
         DO ISIDE=-1,0
            ZZ_GET = 0._EB
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(I,J,K+1+ISIDE,1:N_TRACKED_SPECIES)
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               ZZ_GET(1:N_TRACKED_SPECIES) =  &
                      PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
               (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            END SELECT
            ! KP on low-high side cells:
            CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MU(I,J,K+1+ISIDE),&
                                                   MU_DNS(I,J,K+1+ISIDE),TMPV(ISIDE),KPV(ISIDE))
         ENDDO

         KPDTDN = (CCM1*KPV(-1)+CCP1*KPV(0)) * (TMPV(0)-TMPV(-1)) * IDX

         ! Add contribution to DP:
         ! Low side cell:
         DO ISIDE=-1,0
            FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
            SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
            CASE(CC_FTYPE_RGGAS) ! Regular cell
               DPVOL(I,J,K+1+ISIDE) = DPVOL(I,J,K+1+ISIDE) + FCT*KPDTDN * AF ! +ve or -ve dot
            CASE(CC_FTYPE_CFGAS) ! Cut-cell
               ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
               IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
               JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
               CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*KPDTDN * AF ! +ve or -ve dot
            END SELECT
         ENDDO

   ENDSELECT

ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE
   IW = CUT_FACE(ICF)%IWC
   ! Note: for cut-faces open boundaries are dealt with below in external BC loops:
   IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                           WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE
   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
   SELECT CASE(X1AXIS)
   CASE(IAXIS)
      MUV(-1:0)    = MU(I:I+1,J,K)
      MU_DNSV(-1:0)= MU_DNS(I:I+1,J,K)
   CASE(JAXIS)
      MUV(-1:0)    = MU(I,J:J+1,K)
      MU_DNSV(-1:0)= MU_DNS(I,J:J+1,K)
   CASE(KAXIS)
      MUV(-1:0)    = MU(I,J,K:K+1)
      MU_DNSV(-1:0)= MU_DNS(I,J,K:K+1)
   END SELECT
   DO IFACE=1,CUT_FACE(ICF)%NFACE
      AF = CUT_FACE(ICF)%AREA(IFACE)
      X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
      IDX= 1._EB/ ( CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE) - &
                    CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE) )
      CCM1= IDX*(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F)
      CCP1= IDX*(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))
      ! Interpolate D_Z to the face, linear interpolation:
      TMPV(-1:0)  = -1._EB
      DO ISIDE=-1,0
         ZZ_GET = 0._EB
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
            ZZ_GET(1:N_TRACKED_SPECIES) =  &
                   PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
            (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
         END SELECT
         ! KP on low-high side cells:
         CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MUV(ISIDE),MU_DNSV(ISIDE),TMPV(ISIDE),KPV(ISIDE))
      ENDDO
      KPDTDN = (CCM1*KPV(-1)+CCP1*KPV(0)) * (TMPV(0)-TMPV(-1)) * IDX
      ! Add to divergence integral of surrounding cut-cells:
      DO ISIDE=-1,0
         FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*KPDTDN * AF ! +ve or -ve dot
         END SELECT
      ENDDO
   ENDDO ! IFACE
ENDDO ! ICF

! Now do Boundary conditions for Conductive Heat Flux:
! IAXIS faces:
X1AXIS = IAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE ! Already done on previous loops.
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   AF  = DY(JJG)*DZ(KKG)
   ! Q_LEAK accounts for enthalpy moving through leakage paths
   DPVOL(IIG,JJG,KKG) = DPVOL(IIG,JJG,KKG) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF  + B1%Q_LEAK * (DX(IIG)*DY(JJG)*DZ(KKG))
ENDDO

! JAXIS faces:
X1AXIS = JAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE ! Already done on previous loops.
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   AF  = DX(IIG)*DZ(KKG)
   ! Q_LEAK accounts for enthalpy moving through leakage paths
   DPVOL(IIG,JJG,KKG) = DPVOL(IIG,JJG,KKG) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF  + B1%Q_LEAK * (DX(IIG)*DY(JJG)*DZ(KKG))
ENDDO

! KAXIS faces:
X1AXIS = KAXIS
DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
   IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE ! Already done on previous loops.
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   AF  = DX(IIG)*DY(JJG)
   ! Q_LEAK accounts for enthalpy moving through leakage paths
   DPVOL(IIG,JJG,KKG) = DPVOL(IIG,JJG,KKG) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF  + B1%Q_LEAK * (DX(IIG)*DY(JJG)*DZ(KKG))
ENDDO

! Regular faces connecting gasphase - cut-cells:
DO IFACE=1,MESHES(NM)%CC_NBBRCFACE_Z
   IW = RC_FACE(IFACE)%IWC
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
   X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   BC => BOUNDARY_COORD(WC%BC_INDEX)
   IIG = BC%IIG
   JJG = BC%JJG
   KKG = BC%KKG
   IOR = BC%IOR
   SELECT CASE(X1AXIS)
       CASE(IAXIS)
          AF=DY(JJG)*DZ(KKG)
       CASE(JAXIS)
          AF=DX(IIG)*DZ(KKG)
       CASE(KAXIS)
          AF=DX(IIG)*DY(JJG)
   END SELECT
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
   CASE(CC_FTYPE_RGGAS) ! Regular cell.
      ! Q_LEAK accounts for enthalpy moving through leakage paths
      DPVOL(IIG,JJG,KKG) = DPVOL(IIG,JJG,KKG) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF  + B1%Q_LEAK * (DX(IIG)*DY(JJG)*DZ(KKG))
   CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
      ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
      IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
      JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
      CUT_CELL(ICC)%DVOL(JCC) = &
      CUT_CELL(ICC)%DVOL(JCC) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF + B1%Q_LEAK * CUT_CELL(ICC)%VOLUME(JCC) ! Qconf +ve sign is
                                                                                                       ! outwards of cut-cell.
   END SELECT
ENDDO

! GASPHASE cut-faces:
DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
   IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE .OR. MESHES(NM)%CUT_FACE(ICF)%IWC<1) CYCLE
   IW = MESHES(NM)%CUT_FACE(ICF)%IWC
   IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
       WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
       WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE

   I = CUT_FACE(ICF)%IJK(IAXIS)
   J = CUT_FACE(ICF)%IJK(JAXIS)
   K = CUT_FACE(ICF)%IJK(KAXIS)
   X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
   WC => WALL(IW)
   B1 => BOUNDARY_PROP1(WC%B1_INDEX)
   IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR
   ! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
   !                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
   ISIDE = -1 + (SIGN(1,IOR)+1) / 2
   ! External boundary cut-cells of type OPEN_BOUNDARY:
   GASBOUND_IF : IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
      FCT = -REAL(2*ISIDE+1,EB) ! Factor to set +ve or -ve sign of dot with normal outside.
      DO IFACE=1,CUT_FACE(ICF)%NFACE
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            AF = CUT_FACE(ICF)%AREA(IFACE)
            X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
            IF (IOR > 0) THEN
               IDX= 0.5_EB/(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F) ! Assumes DX twice the distance from WALL_CELL
                                                                      ! to internal cut-cell centroid.
            ELSE
               IDX= 0.5_EB/(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))
            ENDIF
            SELECT CASE(X1AXIS)
            CASE(IAXIS)
               MUV(-1:0)    =     MU(I:I+1,J,K)
               MU_DNSV(-1:0)= MU_DNS(I:I+1,J,K)
               KPV(-1:0)    =     MU(I:I+1,J,K)*CPOPR
               TMPV(-1:0)   =    TMP(I:I+1,J,K)
            CASE(JAXIS)
               MUV(-1:0)    =     MU(I,J:J+1,K)
               MU_DNSV(-1:0)= MU_DNS(I,J:J+1,K)
               KPV(-1:0)    =     MU(I,J:J+1,K)*CPOPR
               TMPV(-1:0)   =    TMP(I,J:J+1,K)
            CASE(KAXIS)
               MUV(-1:0)    =     MU(I,J,K:K+1)
               MU_DNSV(-1:0)= MU_DNS(I,J,K:K+1)
               KPV(-1:0)    =     MU(I,J,K:K+1)*CPOPR
               TMPV(-1:0)   =    TMP(I,J,K:K+1)
            END SELECT
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
            ZZ_GET(1:N_TRACKED_SPECIES) =  PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                                    (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
            CALL GET_CC_CELL_CONDUCTIVITY(ZZ_GET,MUV(ISIDE),MU_DNSV(ISIDE),TMPV(ISIDE),KPV(ISIDE))
            KPDTDN = 0.5_EB*(KPV(-1)+KPV(0)) * (TMPV(0)-TMPV(-1)) * IDX
            CUT_CELL(ICC)%DVOL(JCC) = CUT_CELL(ICC)%DVOL(JCC) + FCT*KPDTDN * AF ! +ve or -ve dot
         END SELECT
      ENDDO

   ELSE
      ! Other boundary conditions:
      DO IFACE=1,CUT_FACE(ICF)%NFACE
         AF = CUT_FACE(ICF)%AREA(IFACE)
         SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
         CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
            ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
            IF (ICC > MESHES(NM)%N_CUTCELL_MESH) CYCLE ! Cut-cell is guard-cell cc.
            JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
            CUT_CELL(ICC)%DVOL(JCC) = &
            CUT_CELL(ICC)%DVOL(JCC) - ( B1%AREA_ADJUST*B1%Q_CON_F ) * AF + B1%Q_LEAK * CUT_CELL(ICC)%VOLUME(JCC) !Qconf +ve sgn
                                                                                                     ! is outwards of cut-cell.
         END SELECT
      ENDDO
   ENDIF GASBOUND_IF
ENDDO

! INBOUNDARY cut-faces, loop on CFACE to add BC defined at SOLID phase:
IF (PREDICTOR) THEN
  DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
     CFA  => CFACE(ICF)
     B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
     IND1 = CFA%CUT_FACE_IND1;                         IND2 = CFA%CUT_FACE_IND2
     ICC  = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2);  JCC  = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
     CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)-( B1%AREA_ADJUST*B1%Q_CON_F ) * CUT_FACE(IND1)%AREA(IND2) !QCONF(+) into solid.
  ENDDO
ELSE
  DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
     CFA  => CFACE(ICF)
     B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
     IND1 = CFA%CUT_FACE_IND1;                         IND2 = CFA%CUT_FACE_IND2
     ICC  = CUT_FACE(IND1)%CELL_LIST(2,LOW_IND,IND2);  JCC  = CUT_FACE(IND1)%CELL_LIST(3,LOW_IND,IND2)
     CUT_CELL(ICC)%DVOL(JCC)=CUT_CELL(ICC)%DVOL(JCC)-( B1%AREA_ADJUST*B1%Q_CON_F ) * CUT_FACE(IND1)%AREA(IND2) !QCONF(+) into solid.
  ENDDO
ENDIF

RETURN
END SUBROUTINE CC_CONDUCTION_HEAT_FLUX


END SUBROUTINE CC_DIVERGENCE_PART_1


! ----------------------- CC_DIFFUSIVE_MASS_FLUXES -------------------------

SUBROUTINE CC_DIFFUSIVE_MASS_FLUXES(NM)

INTEGER, INTENT(IN) :: NM

! NOTE: this routine assumes POINT_TO_MESH(NM) has been previously called.

! Local Variables:
INTEGER :: N,I,J,K,X1AXIS,ISIDE,IFACE,ICC,JCC,ICF
REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP
REAL(EB) :: D_Z_N(0:MAX_I_MAX_TEMP),CCM1,CCP1,IDX,DIFF_FACE,D_Z_TEMP(-1:0),MUV(-1:0),MU_DNSV(-1:0), &
            RHOPV(-1:0),TMPV(-1:0),ZZPV(-1:0),X1F,PRFCT
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET,RHO_D_DZDN_GET
INTEGER,  ALLOCATABLE, DIMENSION(:) :: N_ZZ_MAX_V
INTEGER :: N_LOOKUP, IW
REAL(EB) :: RHO_D_DZDN, ZZ_FACE, TMP_FACE

SELECT CASE(PREDICTOR)
   CASE(.TRUE.)
      ZZP => ZZS
      RHOP => RHOS
      PRFCT = 0._EB ! Use star cut-cell quantities.
   CASE(.FALSE.)
      ZZP => ZZ
      RHOP => RHO
      PRFCT = 1._EB ! Use end of step cut-cell quantities.
END SELECT

ALLOCATE(ZZ_GET(N_TRACKED_SPECIES),RHO_D_DZDN_GET(N_TRACKED_SPECIES))

! Define species index of max CFACE mass fraction.
ALLOCATE(N_ZZ_MAX_V(N_EXTERNAL_CFACE_CELLS+N_INTWALL_CFACE_CELLS+N_INTERNAL_CFACE_CELLS))
DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
   CFA => CFACE(ICF)
   B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
   N_ZZ_MAX_V(ICF)=MAXLOC(B1%ZZ_F(1:N_TRACKED_SPECIES),1)
ENDDO

! 1. Diffusive Heat flux = - Grad dot (h_s rho D Grad Z_n):
! In FV form: use faces to add corresponding face integral terms, for face k
! (sum_a{h_{s,a} rho D_a Grad z_a) dot \hat{n}_k A_k, where \hat{n}_k is the versor outside of cell
! at face k.
DIFFUSIVE_FLUX_LOOP: DO N=1,N_TOTAL_SCALARS

   ! Diffusivity lookup table for species N:
   N_LOOKUP = N
   D_Z_N(:) = D_Z(:,N_LOOKUP)

   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   DO IFACE=1,MESHES(NM)%CC_NRCFACE_Z

      IW = RC_FACE(IFACE)%IWC
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

      I      = RC_FACE(IFACE)%IJK(IAXIS)
      J      = RC_FACE(IFACE)%IJK(JAXIS)
      K      = RC_FACE(IFACE)%IJK(KAXIS)
      X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)

      SELECT CASE(X1AXIS)
         CASE(IAXIS)
            X1F= MESHES(NM)%X(I)
            IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                            RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
            ! Linear interpolation coefficients:
            CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
            CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

            TMPV(-1:0)  = TMP(I:I+1,J,K)
            RHOPV(-1:0) = RHOP(I:I+1,J,K)
            ZZPV(-1:0)  = ZZP(I:I+1,J,K,N)
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ! TMPV(ISIDE) = TMPV(ISIDE)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
                  RHOPV(ISIDE)=        PRFCT *CUT_CELL(ICC)%RHO(JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
                  ZZPV(ISIDE) =        PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               CALL GET_CC_CELL_DIFFUSIVITY(RHOPV(ISIDE),D_Z_N,MU(I+1+ISIDE,J,K),&
                                                MU_DNS(I+1+ISIDE,J,K),TMPV(ISIDE),D_Z_TEMP(ISIDE))
            ENDDO

         CASE(JAXIS)
            X1F= MESHES(NM)%Y(J)
            IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                            RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
            ! Linear interpolation coefficients:
            CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
            CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

            TMPV(-1:0)  = TMP(I,J:J+1,K)
            RHOPV(-1:0) = RHOP(I,J:J+1,K)
            ZZPV(-1:0)  = ZZP(I,J:J+1,K,N)
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ! TMPV(ISIDE) = TMPV(ISIDE)
                  ! RHOPV(ISIDE)= RHOPV(ISIDE)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
                  RHOPV(ISIDE)=        PRFCT *CUT_CELL(ICC)%RHO(JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
                  ZZPV(ISIDE) =        PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               CALL GET_CC_CELL_DIFFUSIVITY(RHOPV(ISIDE),D_Z_N,MU(I,J+1+ISIDE,K),&
                                                MU_DNS(I,J+1+ISIDE,K),TMPV(ISIDE),D_Z_TEMP(ISIDE))
            ENDDO

         CASE(KAXIS)
            X1F= MESHES(NM)%Z(K)
            IDX = 1._EB / ( RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND) - &
                            RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND) )
            ! Linear interpolation coefficients:
            CCM1 = IDX*(RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
            CCP1 = IDX*(X1F -RC_FACE(IFACE)%XCEN(X1AXIS,LOW_IND))

            TMPV(-1:0)  = TMP(I,J,K:K+1)
            RHOPV(-1:0) = RHOP(I,J,K:K+1)
            ZZPV(-1:0)  = ZZP(I,J,K:K+1,N)
            DO ISIDE=-1,0
               SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
               CASE(CC_FTYPE_RGGAS) ! Regular cell -> use stored TMPV from TMP array.
                  ! TMPV(ISIDE) = TMPV(ISIDE)
                  ! RHOPV(ISIDE)= RHOPV(ISIDE)
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
                  JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
                  TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
                  RHOPV(ISIDE)=        PRFCT *CUT_CELL(ICC)%RHO(JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
                  ZZPV(ISIDE) =        PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               CALL GET_CC_CELL_DIFFUSIVITY(RHOPV(ISIDE),D_Z_N,MU(I,J,K+1+ISIDE),&
                                                MU_DNS(I,J,K+1+ISIDE),TMPV(ISIDE),D_Z_TEMP(ISIDE))
            ENDDO

      ENDSELECT

      ! One Term defined flux:
      DIFF_FACE = CCM1*RHOPV(-1)*D_Z_TEMP(-1) + CCP1*RHOPV(0)*D_Z_TEMP(0)
      RC_FACE(IFACE)%RHO_D_DZDN(N) = DIFF_FACE*IDX*(ZZPV(0) - ZZPV(-1) ) ! + rho D_a Grad(Y_a)
      RC_FACE(IFACE)%ZZ_FACE(N) = CCM1*ZZPV(-1) + CCP1*ZZPV(0) ! Linear interpolation of ZZ to the face.
      RC_FACE(IFACE)%TMP_FACE = CCM1*TMPV(-1) + CCP1*TMPV(0)   ! Linear interpolation of Temp to the face.

   ENDDO


   ! GASPHASE cut-faces:
   DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH

      IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE ) CYCLE
      IW = CUT_FACE(ICF)%IWC
      ! Note: for cut-faces open boundaries are dealt with below in external BC loops:
      IF((IW > 0) .AND. .NOT.(WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
                              WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY)) CYCLE

      I = CUT_FACE(ICF)%IJK(IAXIS)
      J = CUT_FACE(ICF)%IJK(JAXIS)
      K = CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)

      DO IFACE=1,CUT_FACE(ICF)%NFACE

         !AF = CUT_FACE(ICF)%AREA(IFACE)
         X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
         IDX= 1._EB/ ( CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE) - &
                       CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE) )
         CCM1= IDX*(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F)
         CCP1= IDX*(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))

         SELECT CASE (X1AXIS)
         CASE(IAXIS)
            MUV(-1:0)     =     MU(I:I+1,J,K)
            MU_DNSV(-1:0) = MU_DNS(I:I+1,J,K)
         CASE(JAXIS)
            MUV(-1:0)     =     MU(I,J:J+1,K)
            MU_DNSV(-1:0) = MU_DNS(I,J:J+1,K)
         CASE(KAXIS)
            MUV(-1:0)     =     MU(I,J,K:K+1)
            MU_DNSV(-1:0) = MU_DNS(I,J,K:K+1)
         END SELECT

         ! Interpolate D_Z to the face, linear interpolation:
         TMPV(-1:0)  = -1._EB; RHOPV(-1:0) = -1._EB; ZZPV(-1:0)  = -1._EB
         DO ISIDE=-1,0
            SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
            CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
               ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
               JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
               TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
               RHOPV(ISIDE)=        PRFCT *CUT_CELL(ICC)%RHO(JCC) + &
                             (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
               ZZPV(ISIDE) =        PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                             (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
            END SELECT
            CALL GET_CC_CELL_DIFFUSIVITY(RHOPV(ISIDE),D_Z_N,MUV(ISIDE),MU_DNSV(ISIDE),TMPV(ISIDE),D_Z_TEMP(ISIDE))
         ENDDO

         ! One Term defined flux:
         DIFF_FACE = CCM1*RHOPV(-1)*D_Z_TEMP(-1) + CCP1*RHOPV(0)*D_Z_TEMP(0)
         CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE) = DIFF_FACE*IDX*(ZZPV(0) - ZZPV(-1) ) ! rho D_a Grad(Y_a)
         CUT_FACE(ICF)%ZZ_FACE(N,IFACE) = CCM1*ZZPV(-1) + CCP1*ZZPV(0) ! Linear interpolation of ZZ to the face.
         CUT_FACE(ICF)%TMP_FACE(IFACE)  = CCM1*TMPV(-1) + CCP1*TMPV(0) ! Linear interpolation of TMP to the face.

      ENDDO ! IFACE

   ENDDO ! ICF

   ! Now Boundary Conditions:
   ! CFACES:
   ISIDE=-1
   CFACE_LOOP : DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
      CFA => CFACE(ICF)
      B1 => BOUNDARY_PROP1(CFA%B1_INDEX)
      ! Use external Gas point data for ZZ_G estimation, consistent with CFA%B1%RDN in the finite difference.
      ! Flux fixing done here for CFACEs:
      RHO_D_DZDN = 2._EB*B1%RHO_D_F(N)*(B1%ZZ_G(N)-B1%ZZ_F(N))*B1%RDN
      IF (N==N_ZZ_MAX_V(ICF)) THEN
         ZZ_GET(1:N_TRACKED_SPECIES) = B1%ZZ_G(1:N_TRACKED_SPECIES)
         RHO_D_DZDN_GET(1:N_TRACKED_SPECIES) = &
         2._EB*B1%RHO_D_F(1:N_TRACKED_SPECIES)*( ZZ_GET(1:N_TRACKED_SPECIES) - B1%ZZ_F(1:N_TRACKED_SPECIES))*B1%RDN
         RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(1:N_TRACKED_SPECIES))-RHO_D_DZDN)
      ENDIF
      B1%RHO_D_DZDN_F(N) = RHO_D_DZDN

      ! Now add variables from CFACES to INBOUNDARY cut-faces containers:
      CUT_FACE(CFA%CUT_FACE_IND1)%RHO_D_DZDN(N,CFA%CUT_FACE_IND2) = RHO_D_DZDN
      CUT_FACE(CFA%CUT_FACE_IND1)%ZZ_FACE(N,   CFA%CUT_FACE_IND2) = B1%ZZ_F(N)
      CUT_FACE(CFA%CUT_FACE_IND1)%TMP_FACE(    CFA%CUT_FACE_IND2) = B1%TMP_F
   ENDDO CFACE_LOOP

   ! Mesh Boundaries:
   ! Regular Faces:
   ! For Regular Faces connecting regular cells we use the WALL_CELL array to fill RHO_D_DZDN, in the same way as
   ! done in WALL_LOOP_2 of DIVERGENCE_PART_1 (divg.f90):
   ! IAXIS faces:
   X1AXIS = IAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_IAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
      IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE ! Already done on previous loops.
      CALL GET_BBREGFACE_RHO_D_DZDN
      ! NOTE: Boundary condition diffusive mass fluxes are already made realizable:
      CC_REGFACE_IAXIS_Z(IFACE)%RHO_D_DZDN(N) = RHO_D_DZDN ! Use single value of RHO_D_DZDN
   ENDDO

   ! JAXIS faces:
   X1AXIS = JAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_JAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
      IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      CALL GET_BBREGFACE_RHO_D_DZDN
      ! NOTE: Boundary condition diffusive mass fluxes are already made realizable:
      CC_REGFACE_JAXIS_Z(IFACE)%RHO_D_DZDN(N) = RHO_D_DZDN ! Use single value of RHO_D_DZDN
   ENDDO

   ! KAXIS faces:
   X1AXIS = KAXIS
   DO IFACE=1,MESHES(NM)%CC_NREGFACE_Z(X1AXIS)
      IW = CC_REGFACE_KAXIS_Z(IFACE)%IWC; IF(IW<1) CYCLE
      IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      CALL GET_BBREGFACE_RHO_D_DZDN
      ! NOTE: Boundary condition diffusive mass fluxes are already made realizable:
      CC_REGFACE_KAXIS_Z(IFACE)%RHO_D_DZDN(N) = RHO_D_DZDN ! Use single value of RHO_D_DZDN
   ENDDO

   ! Regular faces connecting gasphase-gasphase or gasphase- cut-cells:
   DO IFACE=1,MESHES(NM)%CC_NBBRCFACE_Z
      IW = MESHES(NM)%RC_FACE(IFACE)%IWC
      IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE
      I      = RC_FACE(IFACE)%IJK(IAXIS)
      J      = RC_FACE(IFACE)%IJK(JAXIS)
      K      = RC_FACE(IFACE)%IJK(KAXIS)
      X1AXIS = RC_FACE(IFACE)%IJK(KAXIS+1)
      CALL GET_BBRCFACE_RHO_D_DZDN
      RC_FACE(IFACE)%RHO_D_DZDN(N) = RHO_D_DZDN
      RC_FACE(IFACE)%ZZ_FACE(N)   = ZZ_FACE
      RC_FACE(IFACE)%TMP_FACE     = TMP_FACE
   ENDDO

   ! GASPHASE cut-faces:
   ! In case of Cut Faces and OPEN boundaries redefine the location of the guard cells with atmospheric conditions:
   DO ICF = 1,MESHES(NM)%N_CUTFACE_MESH
      IF ( CUT_FACE(ICF)%STATUS /= CC_GASPHASE .OR. MESHES(NM)%CUT_FACE(ICF)%IWC<1) CYCLE
      IW = MESHES(NM)%CUT_FACE(ICF)%IWC
      IF( WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY         .OR. &
          WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .OR. &
          WALL(IW)%BOUNDARY_TYPE==PERIODIC_BOUNDARY ) CYCLE

      I = CUT_FACE(ICF)%IJK(IAXIS)
      J = CUT_FACE(ICF)%IJK(JAXIS)
      K = CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)

      ! External boundary cut-cells of type OPEN_BOUNDARY:
      GASBOUND_IF : IF(WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) THEN
         ! Run over local cut-faces:
         DO IFACE=1,CUT_FACE(ICF)%NFACE
            X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
            IF (BOUNDARY_COORD(WALL(IW)%BC_INDEX)%IOR > 0) THEN
               IDX= 0.5_EB/(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F) ! Assumes DX twice the distance from WALL_CELL to
                                                                      ! internal cut-cell centroid.
            ELSE
               IDX= 0.5_EB/(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))
            ENDIF
            CCM1= 0.5_EB; CCP1= 0.5_EB
            SELECT CASE (X1AXIS)
            CASE(IAXIS)
               MUV(-1:0)       =     MU(I:I+1,J,K)
               MU_DNSV(-1:0)   = MU_DNS(I:I+1,J,K)
               TMPV(-1:0)      =    TMP(I:I+1,J,K)
               RHOPV(-1:0)     =   RHOP(I:I+1,J,K)
               ZZPV(-1:0)      =    ZZP(I:I+1,J,K,N)
            CASE(JAXIS)
               MUV(-1:0)       =     MU(I,J:J+1,K)
               MU_DNSV(-1:0)   = MU_DNS(I,J:J+1,K)
               TMPV(-1:0)      =    TMP(I,J:J+1,K)
               RHOPV(-1:0)     =   RHOP(I,J:J+1,K)
               ZZPV(-1:0)      =    ZZP(I,J:J+1,K,N)
            CASE(KAXIS)
               MUV(-1:0)       =     MU(I,J,K:K+1)
               MU_DNSV(-1:0)   = MU_DNS(I,J,K:K+1)
               TMPV(-1:0)      =    TMP(I,J,K:K+1)
               RHOPV(-1:0)     =   RHOP(I,J,K:K+1)
               ZZPV(-1:0)      =    ZZP(I,J,K:K+1,N)
            END SELECT
            ! Interpolate D_Z to the face, linear interpolation:
            DO ISIDE=-1,0
               SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
               CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
                  ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
                  JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
                  TMPV(ISIDE) = CUT_CELL(ICC)%TMP(JCC)
                  RHOPV(ISIDE)=        PRFCT *CUT_CELL(ICC)%RHO(JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%RHOS(JCC)
                  ZZPV(ISIDE) =        PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                                (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
               END SELECT
               CALL GET_CC_CELL_DIFFUSIVITY(RHOPV(ISIDE),D_Z_N,MUV(ISIDE),&
                                                MU_DNSV(ISIDE),TMPV(ISIDE),D_Z_TEMP(ISIDE))
            ENDDO

            ! One Term defined flux:
            DIFF_FACE = CCM1*RHOPV(-1)*D_Z_TEMP(-1) + CCP1*RHOPV(0)*D_Z_TEMP(0)
            CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE) = DIFF_FACE*IDX*(ZZPV(0) - ZZPV(-1) ) ! rho D_a Grad(Y_a)
            CUT_FACE(ICF)%ZZ_FACE(N,IFACE) = CCM1*ZZPV(-1) + CCP1*ZZPV(0) ! Linear interpolation of ZZ to the face.
            CUT_FACE(ICF)%TMP_FACE(IFACE)  = CCM1*TMPV(-1) + CCP1*TMPV(0) ! Linear interpolation of TMP to the face.

         ENDDO ! IFACE

      ELSE

         ! Other boundary conditions:
         DO IFACE=1,CUT_FACE(ICF)%NFACE
            CALL GET_BBCUTFACE_RHO_D_DZDN
            CUT_FACE(ICF)%RHO_D_DZDN(N,IFACE) = RHO_D_DZDN
            CUT_FACE(ICF)%ZZ_FACE(N,IFACE) = ZZ_FACE
            CUT_FACE(ICF)%TMP_FACE(IFACE)  = TMP_FACE
         ENDDO

      ENDIF GASBOUND_IF

   ENDDO ! ICF

   ! Finally INBOUNDARY cut-faces, compute RHO_D_DZDN using CFACES:
   ! TO DO.

   ! Finally EXIM faces -> we use RHO_D_DZDX,Y,Z previously defined on divg.f90:
   ! No need to do anything on this initial DIFFUSIVE_FLUX_LOOP, as consistency already enforced
   ! on divg.f90.

ENDDO DIFFUSIVE_FLUX_LOOP

DEALLOCATE(ZZ_GET,RHO_D_DZDN_GET,N_ZZ_MAX_V)

RETURN

CONTAINS

SUBROUTINE GET_BBREGFACE_RHO_D_DZDN

INTEGER :: IIG, JJG, KKG, IOR, N_ZZ_MAX
REAL(EB) :: RHO_D_DZDN_GET(1:N_TRACKED_SPECIES)
WC => WALL(IW)
B1 => BOUNDARY_PROP1(WC%B1_INDEX)
BC => BOUNDARY_COORD(WC%BC_INDEX)
IIG = BC%IIG
JJG = BC%JJG
KKG = BC%KKG
IOR = BC%IOR
N_ZZ_MAX = MAXLOC(B1%ZZ_F(1:N_TRACKED_SPECIES),1)
RHO_D_DZDN = 2._EB*B1%RHO_D_F(N)*(ZZP(IIG,JJG,KKG,N)-B1%ZZ_F(N))*B1%RDN
IF (N==N_ZZ_MAX) THEN
   RHO_D_DZDN_GET = 2._EB*B1%RHO_D_F(:)*(ZZP(IIG,JJG,KKG,:)-B1%ZZ_F(:))*B1%RDN
   RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(:))-RHO_D_DZDN)
ENDIF

IF (IOR < 0) RHO_D_DZDN = -RHO_D_DZDN ! This is to switch the sign of the spatial derivative in high side boundaries.

END SUBROUTINE GET_BBREGFACE_RHO_D_DZDN

SUBROUTINE GET_BBRCFACE_RHO_D_DZDN

INTEGER :: IIG, JJG, KKG, IOR, N_ZZ_MAX
REAL(EB) :: ZZ_G, ZZ_GV(1:N_TRACKED_SPECIES),RHO_D_DZDN_GET(1:N_TRACKED_SPECIES)

WC => WALL(IW)
B1 => BOUNDARY_PROP1(WC%B1_INDEX)
BC => BOUNDARY_COORD(WC%BC_INDEX)
IIG = BC%IIG
JJG = BC%JJG
KKG = BC%KKG
IOR = BC%IOR
! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
!                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
ISIDE = -1 + (SIGN(1,IOR)+1) / 2
SELECT CASE(RC_FACE(IFACE)%CELL_LIST(1,ISIDE+2))
CASE(CC_FTYPE_RGGAS) ! Regular cell.
   ZZ_G = ZZP(IIG,JJG,KKG,N)
   ZZ_GV(1:N_TRACKED_SPECIES)= ZZP(IIG,JJG,KKG,1:N_TRACKED_SPECIES)
CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
   ICC = RC_FACE(IFACE)%CELL_LIST(2,ISIDE+2)
   JCC = RC_FACE(IFACE)%CELL_LIST(3,ISIDE+2)
   ZZ_G =               PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
   ZZ_GV(1:N_TRACKED_SPECIES)= PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                        (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
END SELECT

SELECT CASE(X1AXIS)
    CASE(IAXIS)
       X1F= MESHES(NM)%X(I)
    CASE(JAXIS)
       X1F= MESHES(NM)%Y(J)
    CASE(KAXIS)
       X1F= MESHES(NM)%Z(K)
END SELECT

IF (IOR > 0) THEN !Cell or cutcell on high side of RC face:
   IDX = 1._EB / (RC_FACE(IFACE)%XCEN(X1AXIS,HIGH_IND)-X1F)
ELSE
   IDX = 1._EB / (X1F-RC_FACE(IFACE)%XCEN(X1AXIS, LOW_IND))
ENDIF

N_ZZ_MAX = MAXLOC(B1%ZZ_F(1:N_TRACKED_SPECIES),1)
RHO_D_DZDN = B1%RHO_D_F(N)*(ZZ_G-B1%ZZ_F(N))*IDX
IF (N==N_ZZ_MAX) THEN
   RHO_D_DZDN_GET = B1%RHO_D_F(:)*(ZZ_GV(:)-B1%ZZ_F(:))*IDX
   RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(:))-RHO_D_DZDN)
ENDIF

IF (IOR < 0) RHO_D_DZDN = -RHO_D_DZDN ! This is to switch the sign of the spatial derivative in high side boundaries.
DIFF_FACE = B1%RHO_D_F(N)/B1%RHO_F
ZZ_FACE   = B1%ZZ_F(N)
TMP_FACE  = B1%TMP_F

END SUBROUTINE GET_BBRCFACE_RHO_D_DZDN


SUBROUTINE GET_BBCUTFACE_RHO_D_DZDN

INTEGER :: IOR, N_ZZ_MAX
REAL(EB) :: ZZ_G, ZZ_GV(1:N_TRACKED_SPECIES),RHO_D_DZDN_GET(1:N_TRACKED_SPECIES)

WC => WALL(IW)
B1 => BOUNDARY_PROP1(WC%B1_INDEX)
IOR = BOUNDARY_COORD(WC%BC_INDEX)%IOR

X1F= CUT_FACE(ICF)%XYZCEN(X1AXIS,IFACE)
IF (IOR > 0) THEN
   IDX= 1._EB/(CUT_FACE(ICF)%XCENHIGH(X1AXIS,IFACE)-X1F)
ELSE
   IDX= 1._EB/(X1F-CUT_FACE(ICF)%XCENLOW(X1AXIS, IFACE))
ENDIF
! This expression is such that when sign of IOR is -1 -> use Low Side cell  -> ISIDE=-1,
!                              when sign of IOR is  1 -> use High Side cell -> ISIDE= 0 .
ISIDE = -1 + (SIGN(1,IOR)+1) / 2
SELECT CASE(CUT_FACE(ICF)%CELL_LIST(1,ISIDE+2,IFACE))
CASE(CC_FTYPE_CFGAS) ! Cut-cell -> use Temperature value from CUT_CELL data struct:
   ICC = CUT_FACE(ICF)%CELL_LIST(2,ISIDE+2,IFACE)
   JCC = CUT_FACE(ICF)%CELL_LIST(3,ISIDE+2,IFACE)
   ZZ_G =               PRFCT *CUT_CELL(ICC)%ZZ(N,JCC) + &
                 (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(N,JCC)
   ZZ_GV(1:N_TRACKED_SPECIES)= PRFCT *CUT_CELL(ICC)%ZZ(1:N_TRACKED_SPECIES,JCC) + &
                        (1._EB-PRFCT)*CUT_CELL(ICC)%ZZS(1:N_TRACKED_SPECIES,JCC)
END SELECT

N_ZZ_MAX = MAXLOC(B1%ZZ_F(1:N_TRACKED_SPECIES),1)
RHO_D_DZDN = B1%RHO_D_F(N)*(ZZ_G-B1%ZZ_F(N))*IDX
IF (N==N_ZZ_MAX) THEN
   RHO_D_DZDN_GET = B1%RHO_D_F(:)*(ZZ_GV(:)-B1%ZZ_F(:))*IDX
   RHO_D_DZDN = -(SUM(RHO_D_DZDN_GET(:))-RHO_D_DZDN)
ENDIF

IF (IOR < 0) RHO_D_DZDN = -RHO_D_DZDN ! This is to switch the sign of the spatial derivative in high side boundaries.
DIFF_FACE = B1%RHO_D_F(N)/B1%RHO_F
ZZ_FACE   = B1%ZZ_F(N)
TMP_FACE  = B1%TMP_F

END SUBROUTINE GET_BBCUTFACE_RHO_D_DZDN

END SUBROUTINE CC_DIFFUSIVE_MASS_FLUXES



! -------------------------- CC_CHECK_DIVERGENCE -----------------------------

SUBROUTINE CC_CHECK_DIVERGENCE(T,DT,PREDVEL)

USE MPI_F08

! This routine is to be used at the end of predictor or corrector:
REAL(EB),INTENT(IN) :: T,DT
LOGICAL, INTENT(IN) :: PREDVEL

! Local Variables:
INTEGER :: NM, I, J, K, ICC, NCELL, JCC, IPZ

REAL(EB):: PRFCT, DIV, RES, DIVVOL, DIV_JCC, VOL, DPCC, DIV2,TLOC,DTLOC

REAL(EB), POINTER, DIMENSION(:,:,:)  :: UP, VP, WP, DP
REAL(EB), ALLOCATABLE, DIMENSION(:)  :: RESMAXV, RESVOLMX
REAL(EB), ALLOCATABLE, DIMENSION(:,:):: DIVMNX, DIVVOLMNX, VOLMNX
INTEGER,  ALLOCATABLE, DIMENSION(:,:):: IJKRM, RESICJCMX
INTEGER,  ALLOCATABLE, DIMENSION(:,:,:):: IJKMNX    , DIVVOLIJKMNX    , DIVVOLICJCMNX     ,DIVICJCMNX
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:):: XYZMNX
INTEGER :: NMV(1), IERR
REAL(EB), ALLOCATABLE, DIMENSION(:)  :: RESMAXV_AUX, RESVOLMX_AUX
REAL(EB), ALLOCATABLE, DIMENSION(:,:):: DIVMNX_AUX, DIVVOLMNX_AUX, VOLMNX_AUX
INTEGER,  ALLOCATABLE, DIMENSION(:,:):: IJKRM_AUX, RESICJCMX_AUX
INTEGER,  ALLOCATABLE, DIMENSION(:,:,:):: IJKMNX_AUX, DIVVOLIJKMNX_AUX, DIVVOLICJCMNX_AUX ,DIVICJCMNX_AUX
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:):: XYZMNX_AUX
REAL(EB), POINTER, DIMENSION(:) :: D_PBAR_DT_P

! Allocate div Containers
ALLOCATE( RESMAXV(NMESHES), DIVMNX(LOW_IND:HIGH_IND,NMESHES), DIVVOLMNX(LOW_IND:HIGH_IND,NMESHES) )
ALLOCATE( IJKRM(MAX_DIM,NMESHES), IJKMNX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES), XYZMNX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES), &
          DIVVOLIJKMNX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES),  DIVVOLICJCMNX(2,LOW_IND:HIGH_IND,1:NMESHES), &
          DIVICJCMNX(2,LOW_IND:HIGH_IND,1:NMESHES) ,  VOLMNX(LOW_IND:HIGH_IND,1:NMESHES) )
ALLOCATE( RESICJCMX(1:2,1:NMESHES), RESVOLMX(1:NMESHES) )

IF(STORE_CUTCELL_DIVERGENCE) CCVELDIV = 1.E6_EB

! Initialize div containers
RESMAXV(1:NMESHES) = 0._EB
DIVMNX(LOW_IND:HIGH_IND,1:NMESHES)    = 0._EB
DIVVOLMNX(LOW_IND:HIGH_IND,1:NMESHES) = 0._EB
IJKRM(IAXIS:KAXIS,1:NMESHES)                         = 0
IJKMNX(IAXIS:KAXIS,LOW_IND:HIGH_IND,1:NMESHES)       = 0
DIVVOLIJKMNX(IAXIS:KAXIS,LOW_IND:HIGH_IND,1:NMESHES) = 0
DIVVOLICJCMNX(1:2,LOW_IND:HIGH_IND,1:NMESHES)        = 0
DIVICJCMNX(1:2,LOW_IND:HIGH_IND,1:NMESHES)           = 0
RESICJCMX(1:2,1:NMESHES)                             = 0
VOLMNX(LOW_IND:HIGH_IND,1:NMESHES)                   = 0._EB
RESVOLMX(1:NMESHES)                                  = 0._EB
XYZMNX(IAXIS:KAXIS,LOW_IND:HIGH_IND,1:NMESHES)       = 0._EB
TLOC = T
DTLOC= DT
! Meshes Loop:
MESHES_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   DIVMNX(HIGH_IND,NM)     = -10000._EB
   DIVMNX(LOW_IND ,NM)     =  10000._EB
   DIVVOLMNX(HIGH_IND,NM)  = -10000._EB
   DIVVOLMNX(LOW_IND ,NM)  =  10000._EB

   CALL POINT_TO_MESH(NM)

   IF (PREDVEL) THEN ! Take divergence from predicted velocities
      UP => US
      VP => VS
      WP => WS
      DP => DS ! Thermodynamic divergence
      D_PBAR_DT_P => D_PBAR_DT_S
      PRFCT= 1._EB
   ELSE ! Take divergence from final velocities
      UP => U
      VP => V
      WP => W
      DP => D !DDT
      D_PBAR_DT_P => D_PBAR_DT
      PRFCT= 0._EB
   ENDIF

   ! First Regular GASPHASE cells:
   DO K=1,KBAR
      DO J=1,JBAR
         LOOP1: DO I=1,IBAR
            IF ( CCVAR(I,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
            IF ( CELL(CELL_INDEX(I,J,K))%SOLID ) CYCLE
            ! 3D Cartesian divergence:
            DIV = (UP(I,J,K)-UP(I-1,J,K))*RDX(I) + &
                  (VP(I,J,K)-VP(I,J-1,K))*RDY(J) + &
                  (WP(I,J,K)-WP(I,J,K-1))*RDZ(K)
            RES = ABS(DIV-DP(I,J,K))
            IF (RES >= RESMAXV(NM)) THEN
               RESMAXV(NM) = RES
               IJKRM(IAXIS:KAXIS,NM)= (/ I,J,K /)
               RESVOLMX(NM) = DX(I)*DY(J)*DZ(K)
            ENDIF
            IF (DIV >= DIVMNX(HIGH_IND,NM)) THEN
               DIVMNX(HIGH_IND,NM) = DIV
               IJKMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ I,J,K /)
               XYZMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ XC(I),YC(J),ZC(K) /)
               VOLMNX(HIGH_IND,NM) = DX(I)*DY(J)*DZ(K)
            ENDIF
            IF (DIV < DIVMNX(LOW_IND ,NM)) THEN
               DIVMNX(LOW_IND ,NM) = DIV
               IJKMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ I,J,K /)
               XYZMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ XC(I),YC(J),ZC(K) /)
               VOLMNX(LOW_IND,NM) = DX(I)*DY(J)*DZ(K)
            ENDIF
            DIVVOL = DIV*DX(I)*DY(J)*DZ(K)
            IF (DIVVOL >= DIVVOLMNX(HIGH_IND,NM)) THEN
               DIVVOLMNX(HIGH_IND,NM) = DIVVOL
               DIVVOLIJKMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ I,J,K /)
            ENDIF
            IF (DIVVOL < DIVVOLMNX(LOW_IND ,NM)) THEN
               DIVVOLMNX(LOW_IND ,NM) = DIVVOL
               DIVVOLIJKMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ I,J,K /)
            ENDIF
         ENDDO LOOP1
      ENDDO
   ENDDO

   ! Then cut-cells:
   ICC_LOOP : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      NCELL  = CUT_CELL(ICC)%NCELL
      I      = CUT_CELL(ICC)%IJK(IAXIS)
      J      = CUT_CELL(ICC)%IJK(JAXIS)
      K      = CUT_CELL(ICC)%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      IPZ = PRESSURE_ZONE(I,J,K)
      DIVVOL = 0._EB
      DPCC   = 0._EB
      VOL    = 0._EB
      JCC_LOOP : DO JCC=1,NCELL
         VOL  = VOL + CUT_CELL(ICC)%VOLUME(JCC)
         CALL GET_VELOC_DIVERGENCE_CUTCELL(MESHES(NM),ICC,JCC, &
                                          PRFCT,DIV_JCC)
         DIVVOL = DIVVOL + DIV_JCC*CUT_CELL(ICC)%VOLUME(JCC)
         ! Thermodynamic divergence * vol:
         DPCC= DPCC + ((1._EB-PRFCT)*CUT_CELL(ICC)%D(JCC)+PRFCT*CUT_CELL(ICC)%DS(JCC))*CUT_CELL(ICC)%VOLUME(JCC)
      ENDDO JCC_LOOP

      DIV = DIVVOL / (DX(I)*DY(J)*DZ(K))
      RES = ABS(DIVVOL-DPCC)/(DX(I)*DY(J)*DZ(K))
      DIV2 = (UP(I,J,K)-UP(I-1,J,K))*RDX(I) + &
             (VP(I,J,K)-VP(I,J-1,K))*RDY(J) + &
             (WP(I,J,K)-WP(I,J,K-1))*RDZ(K)

      IF(STORE_CUTCELL_DIVERGENCE) CCVELDIV(I,J,K) = DIVVOL/VOL

      IF (RES >= RESMAXV(NM)) THEN
         RESMAXV(NM) = RES
         IJKRM(IAXIS:KAXIS,NM)= (/ I,J,K /)
         RESICJCMX(1:2,NM) = (/ ICC, NCELL /)
         RESVOLMX(NM) = VOL !CUT_CELL(ICC)%VOLUME(JCC)
      ENDIF
      IF (DIV >= DIVMNX(HIGH_IND,NM)) THEN
         DIVMNX(HIGH_IND,NM) = DIV
         IJKMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ I,J,K /)
         XYZMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ XC(I),YC(J),ZC(K) /)
         DIVICJCMNX(1:2,HIGH_IND,NM) = (/ ICC, NCELL /)
         VOLMNX(HIGH_IND,NM) = VOL !CUT_CELL(ICC)%VOLUME(JCC)
      ENDIF
      IF (DIV < DIVMNX(LOW_IND ,NM)) THEN
         DIVMNX(LOW_IND ,NM) = DIV
         IJKMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ I,J,K /)
         XYZMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ XC(I),YC(J),ZC(K) /)
         DIVICJCMNX(1:2,LOW_IND,NM) = (/ ICC, NCELL /)
         VOLMNX(LOW_IND,NM) = VOL !CUT_CELL(ICC)%VOLUME(JCC)
      ENDIF
      IF (DIVVOL >= DIVVOLMNX(HIGH_IND,NM)) THEN
         DIVVOLMNX(HIGH_IND,NM) = DIVVOL
         DIVVOLIJKMNX(IAXIS:KAXIS,HIGH_IND,NM) = (/ I,J,K /)
         DIVVOLICJCMNX(1:2,HIGH_IND,NM) = (/ ICC, NCELL/)
      ENDIF
      IF (DIVVOL < DIVVOLMNX(LOW_IND ,NM)) THEN
         DIVVOLMNX(LOW_IND ,NM) = DIVVOL
         DIVVOLIJKMNX(IAXIS:KAXIS,LOW_IND ,NM) = (/ I,J,K /)
         DIVVOLICJCMNX(1:2,LOW_IND,NM) = (/ ICC, NCELL /)
      ENDIF
   ENDDO ICC_LOOP

   IF(STORE_CARTESIAN_DIVERGENCE) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CELL(CELL_INDEX(I,J,K))%SOLID .OR. CCVAR(I,J,K,CC_CGSC)==CC_SOLID) CARTVELDIV(I,J,K) = 0._EB
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Assign max residual and divergence to corresponding location in MESHES(NM):
   RESMAX = RESMAXV(NM)
   IRM = IJKRM(IAXIS,NM)
   JRM = IJKRM(JAXIS,NM)
   KRM = IJKRM(KAXIS,NM)

   DIVMN = DIVMNX(LOW_IND ,NM)
   IMN = IJKMNX(IAXIS,LOW_IND ,NM)
   JMN = IJKMNX(JAXIS,LOW_IND ,NM)
   KMN = IJKMNX(KAXIS,LOW_IND ,NM)

   DIVMX = DIVMNX(HIGH_IND ,NM)
   IMX = IJKMNX(IAXIS,HIGH_IND ,NM)
   JMX = IJKMNX(JAXIS,HIGH_IND ,NM)
   KMX = IJKMNX(KAXIS,HIGH_IND ,NM)

ENDDO MESHES_LOOP

! Here All_Reduce SUM all mesh values to write if GET_CUTCELLS_VERBOSE:
DEBUG_CC_SCALAR_TRANSPORT_IF : IF (DEBUG_CC_SCALAR_TRANSPORT) THEN
   IF (GET_CUTCELLS_VERBOSE) THEN
      IF (N_MPI_PROCESSES>1) THEN
         ! Allocate aux div Containers
         ALLOCATE( RESMAXV_AUX(NMESHES), DIVMNX_AUX(LOW_IND:HIGH_IND,NMESHES), DIVVOLMNX_AUX(LOW_IND:HIGH_IND,NMESHES) )
         ALLOCATE( IJKRM_AUX(MAX_DIM,NMESHES), IJKMNX_AUX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES), &
                   XYZMNX_AUX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES),&
                   DIVVOLIJKMNX_AUX(MAX_DIM,LOW_IND:HIGH_IND,NMESHES),  DIVVOLICJCMNX_AUX(2,LOW_IND:HIGH_IND,1:NMESHES), &
                   DIVICJCMNX_AUX(2,LOW_IND:HIGH_IND,1:NMESHES) ,  VOLMNX_AUX(LOW_IND:HIGH_IND,1:NMESHES) )
         ALLOCATE( RESICJCMX_AUX(1:2,1:NMESHES), RESVOLMX_AUX(1:NMESHES) )
         RESMAXV_AUX(:)           = RESMAXV(:)
         DIVMNX_AUX(:,:)          = DIVMNX(:,:)
         DIVVOLMNX_AUX(:,:)       = DIVVOLMNX(:,:)
         IJKRM_AUX(:,:)           = IJKRM(:,:)
         IJKMNX_AUX(:,:,:)        = IJKMNX(:,:,:)
         XYZMNX_AUX(:,:,:)        = XYZMNX(:,:,:)
         DIVVOLIJKMNX_AUX(:,:,:)  = DIVVOLIJKMNX(:,:,:)
         DIVVOLICJCMNX_AUX(:,:,:) = DIVVOLICJCMNX(:,:,:)
         DIVICJCMNX_AUX(:,:,:)    = DIVICJCMNX(:,:,:)
         VOLMNX_AUX(:,:)          = VOLMNX(:,:)
         RESICJCMX_AUX(:,:)       = RESICJCMX(:,:)
         RESVOLMX_AUX(:)          = RESVOLMX(:)
         ! Reals:
         CALL MPI_ALLREDUCE(RESMAXV_AUX(1) , RESMAXV(1) ,   NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(DIVMNX_AUX(1,1), DIVMNX(1,1), 2*NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(DIVVOLMNX_AUX(1,1), DIVVOLMNX(1,1), 2*NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(VOLMNX_AUX(1,1), VOLMNX(1,1), 2*NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(RESVOLMX_AUX(1), RESVOLMX(1),   NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(XYZMNX_AUX(1,1,1), XYZMNX(1,1,1), MAX_DIM*2*NMESHES, MPI_DOUBLE_PRECISION, MPI_SUM, &
                            MPI_COMM_WORLD, IERR)
         ! Integers:
         CALL MPI_ALLREDUCE(IJKRM_AUX(1,1), IJKRM(1,1), MAX_DIM*NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(IJKMNX_AUX(1,1,1), IJKMNX(1,1,1), MAX_DIM*2*NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(DIVVOLIJKMNX_AUX(1,1,1), DIVVOLIJKMNX(1,1,1), MAX_DIM*2*NMESHES, MPI_INTEGER, MPI_SUM, &
                            MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(DIVVOLICJCMNX_AUX(1,1,1), DIVVOLICJCMNX(1,1,1), 2*2*NMESHES, MPI_INTEGER, MPI_SUM, &
                            MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(DIVICJCMNX_AUX(1,1,1), DIVICJCMNX(1,1,1), 2*2*NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
         CALL MPI_ALLREDUCE(RESICJCMX_AUX(1,1), RESICJCMX(1,1), 2*NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)

         DEALLOCATE(RESMAXV_AUX, DIVMNX_AUX, DIVVOLMNX_AUX, IJKRM_AUX, IJKMNX_AUX, DIVVOLIJKMNX_AUX, DIVVOLICJCMNX_AUX, &
                    DIVICJCMNX_AUX, VOLMNX_AUX, RESICJCMX_AUX, RESVOLMX_AUX, XYZMNX_AUX)
      ENDIF
      IF (MY_RANK==0) THEN
         WRITE(LU_ERR,*) ' '
         WRITE(LU_ERR,*) "N Step    =",ICYC," T, DT=",TLOC,DTLOC
         NMV(1)=MINLOC(DIVMNX(LOW_IND ,1:NMESHES),DIM=1)
         WRITE(LU_ERR,*) "Div Min   =",NMV(1),DIVMNX(LOW_IND ,NMV(1)),IJKMNX(IAXIS:KAXIS,LOW_IND ,NMV(1)),&
         XYZMNX(IAXIS:KAXIS,LOW_IND ,NMV(1)),&
         DIVICJCMNX(1:2,LOW_IND,NMV(1)),VOLMNX(LOW_IND,NMV(1))
         NMV(1)=MAXLOC(DIVMNX(HIGH_IND ,1:NMESHES),DIM=1)
         WRITE(LU_ERR,*) "Div Max   =",NMV(1),DIVMNX(HIGH_IND,NMV(1)),IJKMNX(IAXIS:KAXIS,HIGH_IND,NMV(1)),&
         XYZMNX(IAXIS:KAXIS,HIGH_IND,NMV(1)),&
         DIVICJCMNX(1:2,HIGH_IND,NMV(1)),VOLMNX(HIGH_IND,NMV(1))

         NMV(1)=MAXLOC(RESMAXV(1:NMESHES),DIM=1)
         WRITE(LU_ERR,*) "Res Max   =",NMV(1),RESMAXV(NMV(1)),IJKRM(IAXIS:KAXIS,NMV(1)),RESICJCMX(1:2,NMV(1)),RESVOLMX(NMV(1))

         NMV(1)=MINLOC(DIVVOLMNX(LOW_IND ,1:NMESHES),DIM=1)
         WRITE(LU_ERR,*) "DivVol Min=",NMV(1),DIVVOLMNX(LOW_IND ,NMV(1)),DIVVOLIJKMNX(IAXIS:KAXIS,LOW_IND ,NMV(1)),&
         DIVVOLICJCMNX(1:2,LOW_IND,NMV(1))
         NMV(1)=MAXLOC(DIVVOLMNX(HIGH_IND ,1:NMESHES),DIM=1)
         WRITE(LU_ERR,*) "DivVol Max=",NMV(1),DIVVOLMNX(HIGH_IND,NMV(1)),DIVVOLIJKMNX(IAXIS:KAXIS,HIGH_IND,NMV(1)),&
         DIVVOLICJCMNX(1:2,HIGH_IND,NMV(1))
      ENDIF
   ENDIF
ENDIF DEBUG_CC_SCALAR_TRANSPORT_IF

! DeAllocate div Containers
DEALLOCATE( RESMAXV, DIVMNX, DIVVOLMNX )
DEALLOCATE( IJKRM, IJKMNX, DIVVOLIJKMNX, DIVVOLICJCMNX, RESVOLMX, RESICJCMX )
DEALLOCATE( XYZMNX )
RETURN
END SUBROUTINE  CC_CHECK_DIVERGENCE

END MODULE CC_DIVERGENCE
