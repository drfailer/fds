!> \brief Pure computation kernels extracted from WALL_ROUTINES
!> These routines take TYPE(MESH_TYPE) as an explicit argument instead of relying on MESH_POINTERS.

MODULE WALL_KERNELS

USE PRECISION_PARAMETERS
USE TYPES
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC CALCULATE_RHO_D_F, CALC_DEPOSITION, PYROLYSIS

CONTAINS


!> \brief Calculate the diffusion coefficient, RHO*D, at the boundary
!> \param M Mesh data structure
!> \param B1 Pointer to BOUNDARY_PROP1 derived type variable
!> \param BC Pointer to BOUNDARY_COORD derived type variable
!> \param WALL_INDEX Optional WALL cell index
!> \param CFACE_INDEX Optional immersed boundary (CFACE) index

SUBROUTINE CALCULATE_RHO_D_F(M,B1,BC,WALL_INDEX,CFACE_INDEX)

INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,CFACE_INDEX
REAL(EB) :: MU_G
TYPE(MESH_TYPE), INTENT(INOUT) :: M
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER :: N,ITMP
REAL(EB) :: RSC_LOC

IF (PRESENT(WALL_INDEX)) THEN
   B1%RHO_G = M%RHO(BC%IIG,BC%JJG,BC%KKG)
   MU_G = M%MU(BC%IIG,BC%JJG,BC%KKG)
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   MU_G  = M%CFACE(CFACE_INDEX)%MU_G
ENDIF

SELECT CASE(SIM_MODE)
   CASE DEFAULT
      DO N=1,N_TRACKED_SPECIES
         B1%RHO_D_F(N) = MU_G*RSC_T*B1%RHO_F/B1%RHO_G
      ENDDO
   CASE (LES_MODE)
      ITMP = MIN(I_MAX_TEMP-1,NINT(B1%TMP_F))
      DO N=1,N_TRACKED_SPECIES
         RSC_LOC = RSC_T
         IF (SPECIES_MIXTURE(N)%SC_T_USER>TWENTY_EPSILON_EB) RSC_LOC=1._EB/SPECIES_MIXTURE(N)%SC_T_USER
         B1%RHO_D_F(N) = B1%RHO_F*( D_Z(ITMP,N) + (MU_G-M%MU_DNS(BC%IIG,BC%JJG,BC%KKG))/B1%RHO_G*RSC_LOC )
      ENDDO
   CASE (DNS_MODE)
      ITMP = MIN(I_MAX_TEMP-1,NINT(B1%TMP_F))
      DO N=1,N_TRACKED_SPECIES
         B1%RHO_D_F(N) = B1%RHO_F*D_Z(ITMP,N)
      ENDDO
END SELECT

END SUBROUTINE CALCULATE_RHO_D_F


!> \brief Calculate aerosol deposition onto a solid surface
!> \param M Mesh data structure
!> \param DT Current time step
!> \param BC Pointer to Boundary Coordinate derived type variable
!> \param B1 Pointer to Boundary Property derived type variable
!> \param B2 Pointer to Boundary Property 2 derived type variable
!> \param WALL_INDEX Optional WALL cell index
!> \param CFACE_INDEX Optional immersed boundary (CFACE) index

SUBROUTINE CALC_DEPOSITION(M,DT,BC,B1,B2,WALL_INDEX,CFACE_INDEX)

USE PHYSICAL_FUNCTIONS, ONLY: GET_VISCOSITY,GET_CONDUCTIVITY,CUNNINGHAM
USE GLOBAL_CONSTANTS, ONLY: K_BOLTZMANN,GRAVITATIONAL_DEPOSITION,TURBULENT_DEPOSITION,THERMOPHORETIC_DEPOSITION,GVEC
REAL(EB), INTENT(IN) :: DT
INTEGER, INTENT(IN), OPTIONAL :: WALL_INDEX,CFACE_INDEX
INTEGER:: N
TYPE(MESH_TYPE), INTENT(INOUT) :: M
TYPE(SPECIES_MIXTURE_TYPE), POINTER :: SM
TYPE(SPECIES_TYPE), POINTER :: SS
REAL(EB), PARAMETER :: CS=1.17_EB,CT=2.2_EB,CM=1.146_EB
REAL(EB), PARAMETER :: CM3=3._EB*CM,CS2=CS*2._EB,CT2=2._EB*CT
REAL(EB), PARAMETER :: ZZ_MIN_DEP=1.E-14_EB
REAL(EB) :: U_THERM,U_TURB,MU_G,Y_AEROSOL,ZZ_GET(1:MAX_SPECIES),YDEP,K_G,TMP_FILM,ALPHA,DTMPDX,&
            TAU_PLUS,U_GRAV,D_SOLID,MW_RATIO,KN,KN_FAC,RSUM_G,ZZ_G(1:MAX_SPECIES),TAU_PLUS_C,U_NORMAL
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(BOUNDARY_PROP2_TYPE), POINTER :: B2
TYPE(SURFACE_TYPE), POINTER :: SFD

IF (PRESENT(WALL_INDEX)) THEN
   SFD=>SURFACE(M%WALL(WALL_INDEX)%SURF_INDEX)
   RSUM_G = M%RSUM(BC%IIG,BC%JJG,BC%KKG)
   ZZ_G(1:N_TRACKED_SPECIES) = M%ZZ(BC%IIG,BC%JJG,BC%KKG,1:N_TRACKED_SPECIES)
ELSEIF (PRESENT(CFACE_INDEX)) THEN
   SFD=>SURFACE(M%CFACE(CFACE_INDEX)%SURF_INDEX)
   RSUM_G= M%CFACE(CFACE_INDEX)%RSUM_G
   ZZ_G(1:N_TRACKED_SPECIES) = B1%ZZ_G(1:N_TRACKED_SPECIES)
ENDIF

IF (ANY(SFD%LEAK_PATH>0)) THEN
   U_NORMAL = 0._EB
ELSE
   U_NORMAL = B1%U_NORMAL
ENDIF

SMIX_LOOP: DO N=1,N_TRACKED_SPECIES

   ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,ZZ_G(1:N_TRACKED_SPECIES))
   IF (ZZ_GET(N) < ZZ_MIN_DEP) CYCLE SMIX_LOOP
   SM => SPECIES_MIXTURE(N)
   IF (.NOT.SM%DEPOSITING) CYCLE SMIX_LOOP
   SS => SPECIES(SM%SINGLE_SPEC_INDEX)
   MW_RATIO = SPECIES_MIXTURE(N)%RCON/RSUM_G
   TMP_FILM = 0.5_EB*(B1%TMP_G+B1%TMP_F)
   CALL GET_VISCOSITY(ZZ_GET,MU_G,TMP_FILM)
   CALL GET_CONDUCTIVITY(ZZ_GET,K_G,TMP_FILM)
   ! Kn=2 lambda/d, lambda has sqrt(1/2). kn_fac has 2*sqrt(1/2)=sqrt(4/2)=sqrt(2)
   KN_FAC = MU_G*SQRT(2._EB*PI/(M%PBAR(BC%KKG,M%PRESSURE_ZONE(BC%IIG,BC%JJG,BC%KKG))*B1%RHO_G))
   ALPHA = K_G/SM%CONDUCTIVITY_SOLID
   DTMPDX = B1%HEAT_TRANS_COEF*(B1%TMP_G-B1%TMP_F)/K_G
   U_THERM = 0._EB
   U_TURB = 0._EB
   U_GRAV = 0._EB

   IF (THERMOPHORETIC_DEPOSITION) THEN
      KN = KN_FAC/SM%THERMOPHORETIC_DIAMETER
      U_THERM = CS2*(ALPHA+CT*KN)*CUNNINGHAM(KN)/((1._EB+CM3*KN)*(1+2*ALPHA+CT2*KN)) * MU_G/(B1%TMP_G*B1%RHO_G)*DTMPDX
   ENDIF
   IF (GRAVITATIONAL_DEPOSITION) THEN
      KN = KN_FAC/SM%MEAN_DIAMETER
      U_GRAV = - DOT_PRODUCT(GVEC,BC%NVEC)*CUNNINGHAM(KN)*SM%MEAN_DIAMETER**2*SM%DENSITY_SOLID/(18._EB*MU_G)
      U_GRAV = MAX(0._EB,U_GRAV)  ! Prevent negative settling velocity at downward facing surfaces
   ENDIF

   IF (TURBULENT_DEPOSITION) THEN
      KN = KN_FAC/SM%MEAN_DIAMETER
      TAU_PLUS_C = SM%DENSITY_SOLID*SM%MEAN_DIAMETER**2/18._EB
      TAU_PLUS = TAU_PLUS_C/MU_G**2*B2%U_TAU**2*B1%RHO_G
      IF (TAU_PLUS < 0.2_EB) THEN ! Diffusion regime
         D_SOLID = K_BOLTZMANN*B1%TMP_G*CUNNINGHAM(KN)/(3._EB*PI*MU_G*SM%MEAN_DIAMETER)
         U_TURB = B2%U_TAU * 0.086_EB*(MU_G/B1%RHO_G/D_SOLID)**(-0.7_EB)
      ELSEIF (TAU_PLUS >= 0.2_EB .AND. TAU_PLUS < 22.9_EB) THEN ! Diffusion-impaction regime
         U_TURB = B2%U_TAU * 3.5E-4_EB * TAU_PLUS**2
      ELSE ! Inertia regime
         U_TURB = B2%U_TAU * 0.17_EB
      ENDIF
   ENDIF
   B2%V_DEP = MAX(0._EB,U_THERM+U_TURB+U_GRAV+U_NORMAL)
   IF (B2%V_DEP <= TWENTY_EPSILON_EB) CYCLE SMIX_LOOP
   ZZ_GET = ZZ_GET * B1%RHO_G
   Y_AEROSOL = ZZ_GET(N)
   YDEP = Y_AEROSOL*MIN(1._EB,(B2%V_DEP)*DT*B1%RDN)
   ZZ_GET(N) = Y_AEROSOL - YDEP
   IF (SM%AWM_INDEX > 0) B1%AWM_AEROSOL(SM%AWM_INDEX)= B1%AWM_AEROSOL(SM%AWM_INDEX)+YDEP/B1%RDN
   IF (SS%AWM_INDEX > 0) B1%AWM_AEROSOL(SS%AWM_INDEX)= B1%AWM_AEROSOL(SS%AWM_INDEX)+YDEP/B1%RDN
   M%D_SOURCE(BC%IIG,BC%JJG,BC%KKG) = M%D_SOURCE(BC%IIG,BC%JJG,BC%KKG) - MW_RATIO*YDEP / B1%RHO_G / DT
   M%M_DOT_PPP(BC%IIG,BC%JJG,BC%KKG,N) = M%M_DOT_PPP(BC%IIG,BC%JJG,BC%KKG,N) - YDEP / DT

ENDDO SMIX_LOOP

END SUBROUTINE CALC_DEPOSITION


!> \brief Calculate the pyrolysis of solid and liquid materials
!> \param M Mesh data structure

SUBROUTINE PYROLYSIS(M,N_MATS,MATL_INDEX,SURF_INDEX,IIG,JJG,KKG,TMP_S,TMP_F,Y_O2_F,IOR,&
                     RHO_DOT_OUT,RHO_S,DEPTH,ASH_DEPTH,DX_S,DT_BC,&
                     M_DOT_G_PPP_ADJUST,M_DOT_G_PPP_ACTUAL,M_DOT_S_PPP,Q_DOT_S_PPP,Q_DOT_G_PPP,Q_DOT_O2_PPP,&
                     Q_DOT_PART,M_DOT_PART,T_BOIL_EFF,B_NUMBER,LAYER_INDEX,REMOVE_LAYER,ONE_D,B1,SOLID_CELL_INDEX,&
                     R_DROP,LPU,LPV,LPW)

USE PHYSICAL_FUNCTIONS, ONLY: GET_MASS_FRACTION,GET_VISCOSITY,GET_PARTICLE_ENTHALPY,GET_SPECIFIC_HEAT,&
                              GET_MASS_FRACTION_ALL,GET_EQUIL_DATA,GET_SENSIBLE_ENTHALPY,GET_Y_SURF,GET_FILM_PROPERTIES,&
                              RAYLEIGH_HEAT_FLUX_MODEL,RAYLEIGH_MASS_FLUX_MODEL
USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
TYPE(MESH_TYPE), INTENT(INOUT) :: M
INTEGER, INTENT(IN) :: N_MATS,SURF_INDEX,IIG,JJG,KKG,IOR,LAYER_INDEX
INTEGER, INTENT(IN), OPTIONAL :: SOLID_CELL_INDEX
LOGICAL, INTENT(IN) :: REMOVE_LAYER
REAL(EB), INTENT(OUT), DIMENSION(:,:) :: RHO_DOT_OUT(N_MATS)
REAL(EB), INTENT(IN) :: TMP_S,TMP_F,DT_BC,DEPTH,RHO_S(N_MATS),Y_O2_F,ASH_DEPTH
REAL(EB), INTENT(IN), OPTIONAL :: R_DROP,LPU,LPV,LPW
REAL(EB), INTENT(IN), DIMENSION(NWP_MAX) :: DX_S
REAL(EB), DIMENSION(:) :: ZZ_GET(1:MAX_SPECIES),Y_ALL(1:MAX_SPECIES)
REAL(EB), DIMENSION(:), INTENT(OUT) :: M_DOT_G_PPP_ADJUST(N_TRACKED_SPECIES),M_DOT_G_PPP_ACTUAL(N_TRACKED_SPECIES)
REAL(EB), DIMENSION(:), INTENT(OUT) :: M_DOT_S_PPP(MAX_MATERIALS),Q_DOT_PART(MAX_LPC),M_DOT_PART(MAX_LPC)
REAL(EB), INTENT(OUT) :: Q_DOT_S_PPP,Q_DOT_G_PPP,Q_DOT_O2_PPP,B_NUMBER
REAL(EB), INTENT(INOUT) :: T_BOIL_EFF(MAX_MATERIALS)
INTEGER, INTENT(IN), DIMENSION(:) :: MATL_INDEX(N_MATS)
INTEGER :: N,NN,NNN,J,NS,SMIX_INDEX(MAX_MATERIALS),NWP,NP,NP2,ITMP
TYPE(MATERIAL_TYPE), POINTER :: ML
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(BOUNDARY_ONE_D_TYPE), POINTER :: ONE_D
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
REAL(EB) :: REACTION_RATE,Y_O2,X_O2,MW(MAX_MATERIALS),Y_GAS(MAX_MATERIALS),Y_TMP(MAX_MATERIALS),Y_SV(MAX_MATERIALS),&
            X_SV(MAX_MATERIALS),X_L(MAX_MATERIALS),&
            D_FILM,H_MASS,RE_L,SHERWOOD,MFLUX,MU_FILM,SC_FILM,TMP_FILM,TMP_G,U2,V2,W2,VEL,&
            DR,R_S_0,R_S_1,H_R,H_R_B,H_S_B,H_S,LENGTH_SCALE,SUM_Y_GAS,SUM_Y_SV,NU_O2_CHAR,Y_O2_S,&
            SUM_Y_SV_SMIX(MAX_SPECIES),X_L_SUM,RHO_DOT_EXTRA,MFLUX_MAX,RHO_FILM,CP_FILM,PR_FILM,K_FILM,&
            RHO_DOT,RHO_DOT_REAC(MAX_REACTIONS),RHO_DOT_REAC_SUM,H_MASS_DNS
LOGICAL :: LIQUID(MAX_MATERIALS),SPEC_ID_ALREADY_USED(MAX_MATERIALS),DO_EVAPORATION

B_NUMBER = 0._EB
Q_DOT_S_PPP = 0._EB
Q_DOT_G_PPP = 0._EB
Q_DOT_O2_PPP = 0._EB
M_DOT_S_PPP = 0._EB
M_DOT_G_PPP_ADJUST = 0._EB
M_DOT_G_PPP_ACTUAL = 0._EB
M_DOT_PART = 0._EB
Q_DOT_PART = 0._EB
RHO_DOT_OUT = 0._EB
SF => SURFACE(SURF_INDEX)

! Determine if any liquids are present. If they are, determine if this is a the surface layer.

DO_EVAPORATION = .FALSE.
IF (ANY(MATERIAL(MATL_INDEX(:))%PYROLYSIS_MODEL==PYROLYSIS_LIQUID)) THEN
   IF (PRESENT(SOLID_CELL_INDEX)) THEN
      IF (SOLID_CELL_INDEX==1) DO_EVAPORATION = .TRUE.
   ENDIF
ENDIF

! If this is surface liquid layer, calculate the Spalding B number and other liquid-specific variables

IF_DO_EVAPORATION: IF (DO_EVAPORATION) THEN

   ! Calculate a sum needed to calculate the volume fraction of liquid components

   LIQUID  = .FALSE.
   X_L_SUM = 0._EB
   MATERIAL_LOOP_00: DO N=1,N_MATS
      ML => MATERIAL(MATL_INDEX(N))
      IF (ML%PYROLYSIS_MODEL/=PYROLYSIS_LIQUID) CYCLE MATERIAL_LOOP_00
      IF (RHO_S(N) < TWENTY_EPSILON_EB) CYCLE MATERIAL_LOOP_00
      LIQUID(N) = .TRUE.
      X_L_SUM = X_L_SUM + RHO_S(N)/ML%RHO_S
   ENDDO MATERIAL_LOOP_00

   IF (X_L_SUM < TWENTY_EPSILON_EB) RETURN

   Y_GAS = 0._EB
   SUM_Y_GAS = 0._EB
   SPEC_ID_ALREADY_USED = .FALSE.
   SMIX_INDEX = 0
   X_L = 0._EB
   X_SV = 0._EB

   MATERIAL_LOOP_0: DO N=1,N_MATS

      IF (.NOT.LIQUID(N)) CYCLE MATERIAL_LOOP_0
      ML => MATERIAL(MATL_INDEX(N))

      SMIX_INDEX(N) = MAXLOC(ML%NU_GAS(:,1),1)
      ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,M%ZZ(IIG,JJG,KKG,1:N_TRACKED_SPECIES))

      IF (ML%MW<0._EB) THEN  ! No molecular weight specified; assume the liquid component evaporates into the defined gas species
         MW(N) = SPECIES_MIXTURE(SMIX_INDEX(N))%MW
      ELSE  ! the user has specified a molecular weight for the liquid component
         MW(N) = ML%MW
      ENDIF

      ! Determine the mass fraction of evaporated MATL N in the first gas phase grid cell

      IF (SPECIES_MIXTURE(SMIX_INDEX(N))%SINGLE_SPEC_INDEX > 0) THEN
         CALL GET_MASS_FRACTION_ALL(ZZ_GET,Y_ALL)
         Y_GAS(N) = Y_ALL(SPECIES_MIXTURE(SMIX_INDEX(N))%SINGLE_SPEC_INDEX)
         IF (SPECIES_MIXTURE(SMIX_INDEX(N))%CONDENSATION_SMIX_INDEX > 0) &
               Y_GAS(N) = Y_GAS(N) - ZZ_GET(SPECIES_MIXTURE(SMIX_INDEX(N))%CONDENSATION_SMIX_INDEX)
      ELSE
         Y_GAS(N) = ZZ_GET(SMIX_INDEX(N))
      ENDIF

      ! Determine volume fraction of MATL N in the liquid and then the surface vapor layer

      T_BOIL_EFF(N) = ML%TMP_BOIL
      CALL GET_EQUIL_DATA(MW(N),TMP_F,M%PBAR(KKG,M%PRESSURE_ZONE(IIG,JJG,KKG)),H_R,H_R_B,T_BOIL_EFF(N),X_SV(N),ML%H_R(1,:))
      X_L(N)  = RHO_S(N)/(ML%RHO_S*X_L_SUM)  ! Volume fraction of MATL component N in the liquid
      X_SV(N) = X_L(N)*X_SV(N)               ! Volume fraction of MATL component N in the surface vapor based on Raoult's law

      ! Calculate sums to be used to compute B number
      IF (.NOT.SPEC_ID_ALREADY_USED(N)) SUM_Y_GAS = SUM_Y_GAS + Y_GAS(N)
      SPEC_ID_ALREADY_USED(N) = .TRUE.

   ENDDO MATERIAL_LOOP_0

   ! Convert mole fraction to mass fraction
   CALL GET_Y_SURF(N_MATS,ZZ_GET,X_SV,Y_SV,MW,SMIX_INDEX)

   ! Compute the Spalding B number

   SUM_Y_SV = SUM(Y_SV(1:N_MATS))
   SUM_Y_SV_SMIX = 0._EB
   MATERIAL_LOOP_1: DO N=1,N_MATS
      IF (.NOT.LIQUID(N)) CYCLE MATERIAL_LOOP_1
      SUM_Y_SV_SMIX(SMIX_INDEX(N)) = SUM_Y_SV_SMIX(SMIX_INDEX(N)) + Y_SV(N)
   ENDDO MATERIAL_LOOP_1

   IF (SUM_Y_SV<ONE_M_EPS) THEN
      B_NUMBER = MAX(0._EB,(SUM_Y_SV-SUM(Y_GAS(1:N_MATS)))/(1._EB-SUM_Y_SV))
   ELSE
      B_NUMBER = 1.E6_EB  ! Fictitiously high B number intended to push mass flux to its upper limit
   ENDIF

   ! Compute an effective gas phase mass fraction, Y_GAS, corresponding to each liquid component, N

   Y_TMP = 0._EB
   MATERIAL_LOOP_2: DO N=1,N_MATS
      IF (.NOT.LIQUID(N)) CYCLE MATERIAL_LOOP_2
      IF (SUM_Y_SV_SMIX(SMIX_INDEX(N))>TWENTY_EPSILON_EB) Y_TMP(N) = Y_SV(N)*Y_GAS(N)/SUM_Y_SV_SMIX(SMIX_INDEX(N))
   ENDDO MATERIAL_LOOP_2
   Y_GAS = Y_TMP

   CALL GET_FILM_PROPERTIES(N_MATS,SF%FILM_FACTOR,Y_SV,Y_GAS,SMIX_INDEX,TMP_F,M%TMP(IIG,JJG,KKG),ZZ_GET,&
                           M%PBAR(KKG,M%PRESSURE_ZONE(IIG,JJG,KKG)),TMP_FILM,MU_FILM,K_FILM,CP_FILM,D_FILM,&
                           RHO_FILM,PR_FILM,SC_FILM)

   ! Compute mass transfer coefficient

   H_MASS_IF: IF (SF%HM_FIXED>=0._EB) THEN

      H_MASS = SF%HM_FIXED

   ELSE H_MASS_IF

      SELECT CASE(ABS(IOR))
         CASE(0); H_MASS_DNS = 0._EB
         CASE(1); H_MASS_DNS = 2._EB*D_FILM*M%RDX(IIG)
         CASE(2); H_MASS_DNS = 2._EB*D_FILM*M%RDY(JJG)
         CASE(3); H_MASS_DNS = 2._EB*D_FILM*M%RDZ(KKG)
      END SELECT

      IF (SIM_MODE==DNS_MODE) THEN

         H_MASS = H_MASS_DNS

      ELSE

         IF (PRESENT(LPU) .AND. PRESENT(LPV) .AND. PRESENT(LPW)) THEN
            U2 = 0.5_EB*(M%U(IIG,JJG,KKG)+M%U(IIG-1,JJG,KKG))
            V2 = 0.5_EB*(M%V(IIG,JJG,KKG)+M%V(IIG,JJG-1,KKG))
            W2 = 0.5_EB*(M%W(IIG,JJG,KKG)+M%W(IIG,JJG,KKG-1))
            VEL = SQRT((U2-LPU)**2+(V2-LPV)**2+(W2-LPW)**2)
         ELSE
            VEL = SQRT(2._EB*M%KRES(IIG,JJG,KKG))
         ENDIF
         CALL GET_VISCOSITY(ZZ_GET,MU_FILM,TMP_FILM)
         IF (PRESENT(R_DROP)) THEN
            LENGTH_SCALE = 2._EB*R_DROP
         ELSE
            LENGTH_SCALE = SF%CONV_LENGTH
         ENDIF
         RE_L     = RHO_FILM*VEL*LENGTH_SCALE/MU_FILM
         SELECT CASE(SF%GEOMETRY)
            CASE DEFAULT         ; SHERWOOD = 0.037_EB*SC_FILM**ONTH*RE_L**0.8_EB
            CASE(SURF_SPHERICAL) ; SHERWOOD = 2._EB + 0.6_EB*SC_FILM**ONTH*SQRT(RE_L)
         END SELECT
         H_MASS   = MAX(H_MASS_DNS,SHERWOOD*D_FILM/LENGTH_SCALE)

      ENDIF

   ENDIF H_MASS_IF

ENDIF IF_DO_EVAPORATION


! Calculate reaction rates for liquids, solids and vegetation

MATERIAL_LOOP: DO N=1,N_MATS  ! Loop over all materials in the cell (alpha subscript)

   IF (RHO_S(N) < TWENTY_EPSILON_EB) CYCLE MATERIAL_LOOP  ! If component alpha density is zero, go on to the next material.
   ML => MATERIAL(MATL_INDEX(N))

   REACTION_LOOP_1: DO J=1,ML%N_REACTIONS  ! Loop over the reactions (beta subscript)

      SELECT CASE (ML%PYROLYSIS_MODEL)

         CASE (PYROLYSIS_LIQUID)

            ! Limit the burning rate to (200 kW/m2) / h_g

             MFLUX_MAX = 200.E3_EB/ML%H_R(J,INT(TMP_F))

            ! Calculate the mass flux of liquid component N at the surface if this is a surface cell.

            IF (DO_EVAPORATION) THEN
               IF (B_NUMBER>TWENTY_EPSILON_EB) THEN
                  MFLUX = MAX(0._EB,MIN(MFLUX_MAX,H_MASS*RHO_FILM*LOG(1._EB+B_NUMBER)*(Y_SV(N) + (Y_SV(N)-Y_GAS(N))/B_NUMBER)))
               ELSE
                  MFLUX = 0._EB
               ENDIF
            ELSE
               MFLUX = 0._EB
            ENDIF

            IF (DX_S(SOLID_CELL_INDEX)>TWENTY_EPSILON_EB) THEN

               ! If the liquid temperature (TMP_S) is greater than the boiling temperature of the current liquid component
               ! ((T_BOIL_EFF(N)), calculate the additional mass loss rate of this component (RHO_DOT_EXTRA) necessary to bring
               ! the liquid temperature back to the boiling temperature.

               RHO_DOT_EXTRA = 0._EB
               IF (TMP_S>T_BOIL_EFF(N)) THEN
                  ITMP = MIN(I_MAX_TEMP,INT(TMP_S))
                  H_S = ML%H(ITMP) + (TMP_S-REAL(ITMP,EB))*(ML%H(ITMP+1)-ML%H(ITMP))
                  ITMP = INT(T_BOIL_EFF(N))
                  H_S = H_S - (ML%H(ITMP) + (T_BOIL_EFF(N)-REAL(ITMP,EB))*(ML%H(ITMP+1)-ML%H(ITMP)))
                  H_S = H_S * RHO_S(N)
                  H_R = ML%H_R(1,NINT(T_BOIL_EFF(N)))
                  RHO_DOT_EXTRA = H_S/(H_R*DT_BC)  ! kg/m3/s
               ENDIF

               ! Calculate the mass loss rate per unit volume of this liquid component (RHO_DOT)

               SELECT CASE(SF%GEOMETRY)
                  CASE DEFAULT
                     MFLUX = MIN(MFLUX_MAX,MFLUX + RHO_DOT_EXTRA*DX_S(SOLID_CELL_INDEX))
                     RHO_DOT_REAC(J) =MFLUX/DX_S(SOLID_CELL_INDEX) ! kg/m3/s
                  CASE(SURF_SPHERICAL)
                     NWP = SUM(ONE_D%N_LAYER_CELLS(1:ONE_D%N_LAYERS))
                     R_S_0 = SF%INNER_RADIUS + ONE_D%X(NWP) - ONE_D%X(0)
                     R_S_1 = SF%INNER_RADIUS + ONE_D%X(NWP) - ONE_D%X(1)
                     DR = (R_S_0**3-R_S_1**3)/(3._EB*R_S_0**2)
                     MFLUX = MIN(MFLUX_MAX,MFLUX + RHO_DOT_EXTRA*DR)
                     RHO_DOT_REAC(J) = MFLUX/DR
               END SELECT

            ENDIF

         CASE (PYROLYSIS_SOLID)

            ! Reaction rate in 1/s (Tech Guide: r_alpha_beta)

            REACTION_RATE = ML%A(J)*(RHO_S(N))**ML%N_S(J)*EXP(-ML%E(J)/(R0*TMP_S))

            ! power term

            IF (ABS(ML%N_T(J))>=TWENTY_EPSILON_EB) REACTION_RATE = REACTION_RATE * TMP_S**ML%N_T(J)

            ! Oxidation reaction?

            IF ( ML%N_O2(J)>0._EB .AND. O2_INDEX>0 ) THEN
               ! Calculate oxygen volume fraction at the surface
               X_O2 = SPECIES(O2_INDEX)%RCON*Y_O2_F/M%RSUM(IIG,JJG,KKG)
               ! Calculate oxygen concentration inside the material, assuming decay function
               X_O2 = X_O2 * EXP(-MAX(0._EB,DEPTH-ASH_DEPTH)/(TWENTY_EPSILON_EB+ML%GAS_DIFFUSION_DEPTH(J)))
               REACTION_RATE = REACTION_RATE * X_O2**ML%N_O2(J)
            ENDIF
            REACTION_RATE = MIN(REACTION_RATE,ML%MAX_REACTION_RATE(J))  ! User-specified limit
            RHO_DOT_REAC(J) = REACTION_RATE  ! Tech Guide: rho_s(0)*r_alpha,beta kg/m3/s

         CASE (PYROLYSIS_SURFACE_OXIDATION)

            ! Reaction rate in kg/m2/s
            REACTION_RATE = ML%A(J)*EXP(-ML%E(J)/(R0*TMP_S))

            ! Estimate surface oxygen concentration from mass transport
            TMP_FILM = (TMP_F+M%TMP(IIG,JJG,KKG))/2._EB
            ! Get oxygen mass fraction
            ZZ_GET(1:N_TRACKED_SPECIES) = MAX(0._EB,B1%ZZ_G(1:N_TRACKED_SPECIES))
            CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
            CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_FILM,TMP_FILM)
            ! Mass transfer coefficient
            H_MASS = B1%HEAT_TRANS_COEF/CP_FILM
            ! Mass stoichiometric coefficient for oxygen
            NU_O2_CHAR = ML%NU_GAS_M(O2_INDEX,J)
            Y_O2_S = 0._EB
            IF (H_MASS>0._EB) THEN
               Y_O2_S = (SQRT(4._EB*REACTION_RATE/H_MASS*Y_O2+(REACTION_RATE*NU_O2_CHAR/H_MASS)**2._EB + &
                  2._EB*REACTION_RATE*NU_O2_CHAR/H_MASS+1._EB)-REACTION_RATE*NU_O2_CHAR/H_MASS-1) / &
                  (2._EB*REACTION_RATE/H_MASS)
            ENDIF

            ! Compute LENGTH_SCALE: 1/(surface-to-volume ratio)
            IF (SF%BOUNDARY_FUEL_MODEL) THEN
               LENGTH_SCALE = 1._EB/(SF%SURFACE_VOLUME_RATIO(LAYER_INDEX)*SF%PACKING_RATIO(LAYER_INDEX))
            ELSE
               LENGTH_SCALE = SF%INNER_RADIUS + SUM(ONE_D%LAYER_THICKNESS(1:ONE_D%N_LAYERS))
               SELECT CASE(SF%GEOMETRY)
                  CASE(SURF_SPHERICAL)
                     LENGTH_SCALE = LENGTH_SCALE/3._EB
                  CASE DEFAULT
                     LENGTH_SCALE = LENGTH_SCALE/2._EB
               END SELECT
            ENDIF

            REACTION_RATE = Y_O2_S/LENGTH_SCALE*REACTION_RATE
            REACTION_RATE = MIN(REACTION_RATE,ML%MAX_REACTION_RATE(J))  ! User-specified limit
            RHO_DOT_REAC(J) = MAX(REACTION_RATE,0._EB)  ! Tech Guide: rho_s(0)*r_alpha,beta kg/m3/s

      END SELECT
   ENDDO REACTION_LOOP_1

   RHO_DOT_REAC_SUM = SUM(RHO_DOT_REAC(1:ML%N_REACTIONS))
   IF (RHO_DOT_REAC_SUM < TWENTY_EPSILON_EB .AND. .NOT. REMOVE_LAYER) CYCLE MATERIAL_LOOP

   IF (REMOVE_LAYER) THEN
      RHO_DOT = RHO_S(N)/DT_BC
      ! If layer is being removed but zero reaction rate, apportion mass loss equally over reactions
      IF (RHO_DOT_REAC_SUM < TWENTY_EPSILON_EB) RHO_DOT_REAC_SUM = RHO_DOT/ML%N_REACTIONS
   ELSE
      RHO_DOT = MIN(RHO_DOT_REAC_SUM,RHO_S(N)/DT_BC)
   ENDIF

   RHO_DOT_REAC(1:ML%N_REACTIONS) = RHO_DOT_REAC(1:ML%N_REACTIONS)*RHO_DOT/RHO_DOT_REAC_SUM

   ! Optional limiting of mass loss rate based on specified fuel burnout time

   IF (SF%MINIMUM_BURNOUT_TIME<1.E5_EB) RHO_DOT = MIN(RHO_DOT,SF%LAYER_DENSITY(LAYER_INDEX)/SF%MINIMUM_BURNOUT_TIME)

   REACTION_LOOP_2:DO J=1,ML%N_REACTIONS

      RHO_DOT_OUT(N) = RHO_DOT_OUT(N) + RHO_DOT_REAC(J)

      DO NN=1,ML%N_RESIDUE(J) ! Get residue production (alpha' represents the other materials)
         NNN = FINDLOC(MATL_INDEX,ML%RESIDUE_MATL_INDEX(NN,J),1)
         RHO_DOT_OUT(NNN) = RHO_DOT_OUT(NNN) - ML%NU_RESIDUE(NN,J)*RHO_DOT_REAC(J)
         M_DOT_S_PPP(NNN) = M_DOT_S_PPP(NNN) + ML%NU_RESIDUE(NN,J)*RHO_DOT_REAC(J) ! (m_dot_alpha')'''
      ENDDO

      ! Optional variable heat of reaction

      ITMP = MIN(I_MAX_TEMP,NINT(TMP_S))
      H_R = ML%H_R(J,ITMP)

      ! Calculate various energy and mass source terms

      Q_DOT_S_PPP    = Q_DOT_S_PPP    - RHO_DOT_REAC(J)*H_R ! Tech Guide: q_dot_s'''
      M_DOT_S_PPP(N) = M_DOT_S_PPP(N) - RHO_DOT_REAC(J)     ! m_dot_alpha''' = -rho_s(0) * sum_beta r_alpha,beta
      TMP_G = M%TMP(IIG,JJG,KKG)
      ! Tech Guide: m_dot_gamma'''

      M_DOT_G_PPP_ACTUAL(:) = M_DOT_G_PPP_ACTUAL(:) + ML%NU_GAS(:,J)*RHO_DOT_REAC(J)
      M_DOT_G_PPP_ADJUST(:) = M_DOT_G_PPP_ADJUST(:) + ML%NU_GAS(:,J)*RHO_DOT_REAC(J) * ML%ADJUST_BURN_RATE(:,J)
      DO NS=1,N_SPECIES
         IF (ML%NU_GAS_P(NS,J) <= 0._EB .AND. ML%NU_GAS_M(NS,J) <= 0._EB) CYCLE
         CALL INTERPOLATE1D_UNIFORM(0,SPECIES(NS)%H_G,TMP_S,H_S_B)
         CALL INTERPOLATE1D_UNIFORM(0,SPECIES(NS)%H_G,TMP_G,H_S)
         IF (ML%NU_GAS_P(NS,J) > 0._EB) &
            Q_DOT_G_PPP = Q_DOT_G_PPP + ML%ADJUST_BURN_RATE_P(NS,J)*ML%NU_GAS_P(NS,J)*RHO_DOT_REAC(J)*(H_S-H_S_B)
         IF (ML%NU_GAS_M(NS,J) > 0._EB) &
            Q_DOT_G_PPP = Q_DOT_G_PPP -                             ML%NU_GAS_M(NS,J)*RHO_DOT_REAC(J)*(H_S-H_S_B)
      ENDDO

      IF (ANY(ML%NU_LPC(:,J)>0._EB)) THEN
         DO NP=1,ML%N_LPC(J)
            IF (ML%NU_LPC(NP,J)<=0._EB) CYCLE
            NP2 = ML%LPC_INDEX(NP,J)
            IF (SF%MATL_PART_INDEX(NP2)==NP) THEN
               M_DOT_PART(NP2)=ML%NU_LPC(NP,J)*RHO_DOT_REAC(J)
               Q_DOT_PART(NP2)=GET_PARTICLE_ENTHALPY(NP,TMP_S)*M_DOT_PART(NP2)
            ENDIF
         ENDDO
      ENDIF

      ! If there is char oxidation, save the HRR per unit volume generated

      IF (ML%N_O2(J)>0._EB) Q_DOT_O2_PPP = Q_DOT_O2_PPP - RHO_DOT_REAC(J)*H_R

   ENDDO REACTION_LOOP_2

ENDDO MATERIAL_LOOP

END SUBROUTINE PYROLYSIS


END MODULE WALL_KERNELS
