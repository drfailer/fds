!  +++++++++++++++++++++++ CC_EXCHANGE ++++++++++++++++++++++++++

! MPI exchange routines for cut-cell / immersed-boundary data.

MODULE CC_EXCHANGE

USE CC_SCALARS_DATA
USE CC_SCALARS, ONLY: CC_H_INTERP, CC_RHO0W_INTERP
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: MESH_CC_EXCHANGE

CONTAINS



! ------------------------------- MESH_CC_EXCHANGE ---------------------------------

SUBROUTINE MESH_CC_EXCHANGE(CODE)

USE MPI_F08

INTEGER, INTENT(IN) :: CODE

! Local Variables:
INTEGER :: NM,NOM,RNODE,SNODE,IERR
INTEGER :: II1,JJ1,KK1,NCELL,ICC,ICC1,JCC1,NQT2,JCC,LL,NN
INTEGER :: I,J,K,II,JJ,KK,IFC,ICF,X1AXIS,ICF1,JCF
TYPE (MESH_TYPE), POINTER :: M,M1
TYPE (OMESH_TYPE), POINTER :: M2,M3
REAL(EB), POINTER, DIMENSION(:,:,:) :: UP,UP2,VP,VP2,WP,WP2
LOGICAL, SAVE :: INITIALIZE_CC_SCALARS_FORC=.TRUE.

INTEGER :: EP,INPE,INT_NPE_LO,INT_NPE_HI,VIND,ICELL,IEDGE,IFEP,IW,IIO,JJO,KKO
! For solid phase only return. All variables exchanged currently here are gas-phase.
IF (SOLID_PHASE_ONLY) RETURN
! In case of initialization code from main return.
! Initialization of cut-cell communications needs to be done later in the main.f90 sequence and will be done using
! INITIALIZE_CC_SCALARS/VELOCITY logicals.
IF (CODE == 0 .OR. CODE==2 .OR. CODE>6) RETURN
! No need to do mesh exchange within pressure iteration scheme here, when no IBM forcing, or call to fill GLMAT H ghost cells.
IF (.NOT.CC_MATVEC_DEFINED) RETURN
IF (CODE == 3 .AND. CALL_FROM_GLMAT_SETUP) RETURN

! First Allocate and setup persistent send-receives for scalars:
INITIALIZE_CC_SCALARS_FORC_COND : IF (INITIALIZE_CC_SCALARS_FORC) THEN

   ! Allocate REQ11, for scalar transport quantities, reduced cycling conditionals:
   N_REQ11=0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NICC_S(1)==0 .AND. M3%NICC_R(1)==0) CYCLE
         N_REQ11 = N_REQ11+1
      ENDDO
   ENDDO
   ALLOCATE(REQ11(N_REQ11*4)); N_REQ11=0


   ! Allocate REQ112: Exchange cut-face data in block boundaries.
   N_REQ112=0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NICF_S(1)==0 .AND. M3%NICF_R(1)==0 .AND. M3%NLKF_S==0 .AND. M3%NLKF_R==0) CYCLE
         N_REQ112 = N_REQ112+1
      ENDDO
   ENDDO
   ALLOCATE(REQ112(N_REQ112*4)); N_REQ112=0

   ! Allocate REQ12: Dual use, IBM forcing or TAU,OMG computation vars.
   N_REQ12=0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_S(1)==0 .AND. M3%NFCC_R(1)==0) CYCLE
         N_REQ12 = N_REQ12+1
      ENDDO
   ENDDO
   ALLOCATE(REQ12(N_REQ12*4)); N_REQ12=0

   ! Allocate REQ13, for end of step H and RHO_0*W interpolation (cell) quantities:
   N_REQ13=0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_S(2)==0 .AND. M3%NFCC_R(2)==0) CYCLE
         N_REQ13 = N_REQ13+1
      ENDDO
   ENDDO
   ALLOCATE(REQ13(N_REQ13*4)); N_REQ13=0


   ! 1. Receives:
   MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      RNODE = PROCESS(NM)

      ! REQ11:
      OTHER_MESH_LOOP_11: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NICC_R(1)==0) CYCLE OTHER_MESH_LOOP_11
         SNODE = PROCESS(NOM)
         IF (M3%NICC_R(1)>0) THEN
            ! Cell centered variables on cut-cells:
            ALLOCATE(M3%REAL_RECV_PKG11(M3%NICC_R(2)*(4+N_TOTAL_SCALARS)))
            IF (RNODE/=SNODE) THEN
               N_REQ11 = N_REQ11 + 1
               CALL MPI_RECV_INIT(M3%REAL_RECV_PKG11(1),SIZE(M3%REAL_RECV_PKG11),MPI_DOUBLE_PRECISION, &
                                  SNODE,NOM,MPI_COMM_WORLD,REQ11(N_REQ11),IERR)
            ENDIF
         ENDIF
      ENDDO OTHER_MESH_LOOP_11

      ! REQ112:
      OTHER_MESH_LOOP_112: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NICF_R(1)==0 .AND. M3%NLKF_R==0) CYCLE OTHER_MESH_LOOP_112
         SNODE = PROCESS(NOM)
         ! Cut-face centered variables VEL/VELS, F, FB, ICG Hi-1,Hi:
         ALLOCATE(M3%REAL_RECV_PKG112(M3%NICF_R(2) * 4 + M3%NLKF_R * 3))
         IF (RNODE/=SNODE) THEN
            N_REQ112 = N_REQ112 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG112(1),SIZE(M3%REAL_RECV_PKG112),MPI_DOUBLE_PRECISION, &
                               SNODE,NOM,MPI_COMM_WORLD,REQ112(N_REQ112),IERR)
         ENDIF
      ENDDO OTHER_MESH_LOOP_112

      ! REQ12:
      OTHER_MESH_LOOP_12: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_R(1)==0) CYCLE OTHER_MESH_LOOP_12
         SNODE = PROCESS(NOM)
         ! Face centered variables Ux1, Fvx1, dHdx1:
         ALLOCATE(M3%REAL_RECV_PKG12(M3%NFCC_R(1) * 2))
         IF (RNODE/=SNODE) THEN
            N_REQ12 = N_REQ12 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG12(1),SIZE(M3%REAL_RECV_PKG12),MPI_DOUBLE_PRECISION, &
                               SNODE,NOM,MPI_COMM_WORLD,REQ12(N_REQ12),IERR)
         ENDIF
      ENDDO OTHER_MESH_LOOP_12

      ! REQ13:
      OTHER_MESH_LOOP_13: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_R(2)==0) CYCLE OTHER_MESH_LOOP_13
         SNODE = PROCESS(NOM)
         ! Cell centered variables:
         ALLOCATE(M3%REAL_RECV_PKG13(M3%NFCC_R(2)*(NQT2C+N_TRACKED_SPECIES)))
         IF (RNODE/=SNODE) THEN
            N_REQ13 = N_REQ13 + 1
            CALL MPI_RECV_INIT(M3%REAL_RECV_PKG13(1),SIZE(M3%REAL_RECV_PKG13),MPI_DOUBLE_PRECISION, &
                               SNODE,NOM,MPI_COMM_WORLD,REQ13(N_REQ13),IERR)
         ENDIF
      ENDDO OTHER_MESH_LOOP_13

   ENDDO MESH_LOOP_1

   ! 2. Sends:
   SENDING_MESH_LOOP_1: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      RNODE = PROCESS(NM)
      M =>MESHES(NM)

      ! REQ11:
      RECEIVING_MESH_LOOP_11: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         SNODE = PROCESS(NOM)
         ! Initialize persistent send requests
         IF (M3%NICC_S(1)>0 .AND. RNODE/=SNODE) THEN
            ALLOCATE(M3%REAL_SEND_PKG11(M3%NICC_S(2)*(4+N_TOTAL_SCALARS)))
            N_REQ11 = N_REQ11 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG11(1),SIZE(M3%REAL_SEND_PKG11),MPI_DOUBLE_PRECISION, &
                               SNODE,NM,MPI_COMM_WORLD,REQ11(N_REQ11),IERR)
         ENDIF
      ENDDO RECEIVING_MESH_LOOP_11

      ! REQ112:
      RECEIVING_MESH_LOOP_112: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NICF_S(1)==0 .AND. M3%NLKF_S==0)  CYCLE RECEIVING_MESH_LOOP_112
         SNODE = PROCESS(NOM)
         ! Initialize persistent send requests
         IF ((M3%NICF_S(1)>0 .OR. M3%NLKF_S>0) .AND. RNODE/=SNODE) THEN
            ALLOCATE(M3%REAL_SEND_PKG112(M3%NICF_S(2) * 4 + M3%NLKF_S * 3))
            N_REQ112 = N_REQ112 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG112(1),SIZE(M3%REAL_SEND_PKG112),MPI_DOUBLE_PRECISION, &
                               SNODE,NM,MPI_COMM_WORLD,REQ112(N_REQ112),IERR)
         ENDIF
      ENDDO RECEIVING_MESH_LOOP_112

      ! REQ12:
      RECEIVING_MESH_LOOP_12: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_S(1)==0)  CYCLE RECEIVING_MESH_LOOP_12
         SNODE = PROCESS(NOM)
         ! Initialize persistent send requests
         IF (M3%NFCC_S(1)>0 .AND. RNODE/=SNODE) THEN
            ALLOCATE(M3%REAL_SEND_PKG12(M3%NFCC_S(1) * 2))
            N_REQ12 = N_REQ12 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG12(1),SIZE(M3%REAL_SEND_PKG12),MPI_DOUBLE_PRECISION, &
                               SNODE,NM,MPI_COMM_WORLD,REQ12(N_REQ12),IERR)
         ENDIF
      ENDDO RECEIVING_MESH_LOOP_12

      ! REQ13:
      RECEIVING_MESH_LOOP_13: DO NOM=1,NMESHES
         M3=>MESHES(NM)%OMESH(NOM)
         IF (M3%NFCC_S(2)==0)  CYCLE RECEIVING_MESH_LOOP_13
         SNODE = PROCESS(NOM)
         ! Initialize persistent send requests
         IF (M3%NFCC_S(2)>0 .AND. RNODE/=SNODE) THEN
            ALLOCATE(M3%REAL_SEND_PKG13(M3%NFCC_S(2)*(NQT2C+N_TRACKED_SPECIES)))
            N_REQ13 = N_REQ13 + 1
            CALL MPI_SEND_INIT(M3%REAL_SEND_PKG13(1),SIZE(M3%REAL_SEND_PKG13),MPI_DOUBLE_PRECISION, &
                               SNODE,NM,MPI_COMM_WORLD,REQ13(N_REQ13),IERR)
         ENDIF
      ENDDO RECEIVING_MESH_LOOP_13

   ENDDO SENDING_MESH_LOOP_1

   INITIALIZE_CC_SCALARS_FORC = .FALSE.

ENDIF INITIALIZE_CC_SCALARS_FORC_COND


! Exchange Scalars in cut-cells:
SENDING_MESH_LOOP_2: DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M =>MESHES(NM)
   RECEIVING_MESH_LOOP_2: DO NOM=1,NMESHES

      M1=>MESHES(NOM)
      M3=>MESHES(NM)%OMESH(NOM)

      SNODE = PROCESS(NOM)
      RNODE = PROCESS(NM)

      ! Exchange of density and species mass fractions following the PREDICTOR update

      IF (CODE==1 .AND. M3%NICC_S(1)>0) THEN
         NQT2 = 4+N_TOTAL_SCALARS
         LL = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG11: DO ICC1=1,M3%NICC_S(1)
               ICC=M3%ICC_UNKZ_CC_S(ICC1)
               NCELL=M%CUT_CELL(ICC)%NCELL
               II1=M%CUT_CELL(ICC)%IJK(IAXIS)
               JJ1=M%CUT_CELL(ICC)%IJK(JAXIS)
               KK1=M%CUT_CELL(ICC)%IJK(KAXIS)
               DO JCC=1,NCELL
                  LL = LL + 1
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+1) = M%CUT_CELL(ICC)%RHOS(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+2) = M%CUT_CELL(ICC)%TMP(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%RSUM(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+4) = M%CUT_CELL(ICC)%D(JCC)
                  DO NN=1,N_TOTAL_SCALARS
                     M3%REAL_SEND_PKG11(NQT2*(LL-1)+4+NN) = M%CUT_CELL(ICC)%ZZS(NN,JCC)
                  ENDDO
               ENDDO
            ENDDO PACK_REAL_SEND_PKG11
         ENDIF
      ENDIF

      ! Information for cell centered variables:
      IF (CODE==1 .AND. M3%NFCC_S(2)>0) THEN
         NQT2 = NQT2C+N_TOTAL_SCALARS
         LL = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG213 : DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+1) = M%HS(I,J,K)                            ! Prev H in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%W(I,J,K-1)+M%W(I,J,K))       ! Wcen^n in I,J,K.
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+4) = M%RHOS(I,J,K)                          ! RHO^*
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+5) = M%TMP(I,J,K)                           ! TMP^*
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+6) = M%RSUM(I,J,K)                          ! RSUM^*
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+7) = M%MU(I,J,K)                            ! MU^n
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+8) = M%MU_DNS(I,J,K)                        ! MU_DNS^n
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+NQT2C)= M%RHO(I,J,K)*(M%HS(I,J,K)-M%KRES(I,J,K)) ! Previous substep pressure.
               DO NN=1,N_TOTAL_SCALARS
                  M3%REAL_SEND_PKG13(NQT2*(LL-1)+NQT2C+NN)= M%ZZS(I,J,K,NN)
               ENDDO
            ENDDO PACK_REAL_SEND_PKG213
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG213: DO IFEP=1,M2%NFEP_R(2)
               ICF = M2%IFEP_R_2( LOW_IND,IFEP)
               INPE= M2%IFEP_R_2(HIGH_IND,IFEP)
               LL  = M1%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_CC_S(LL)
               J     = M3%JJO_CC_S(LL)
               K     = M3%KKO_CC_S(LL)
               M1%CUT_FACE(ICF)%INT_CVARS(    INT_H_IND,INPE)= M%HS(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(  INT_RHO_IND,INPE)= M%RHOS(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(  INT_TMP_IND,INPE)= M%TMP(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS( INT_RSUM_IND,INPE)= M%RSUM(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(   INT_MU_IND,INPE)= M%MU(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(INT_MUDNS_IND,INPE)= M%MU_DNS(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(    INT_P_IND,INPE)= M%RHO(I,J,K)*(M%HS(I,J,K)-M%KRES(I,J,K))
               DO NN=1,N_TOTAL_SCALARS
                  M1%CUT_FACE(ICF)%INT_CVARS(INT_P_IND+NN,INPE)=M%ZZS(I,J,K,NN)
               ENDDO
            ENDDO PACK_REAL_RECV_PKG213
         ENDIF
      ENDIF

      ! Exchange velocity, momentum rhs and previous substep dH/Dx1 for cut-faces, in PREDICTOR, IBM forcing:
      IF (CODE==5 .AND. PREDICTOR .AND. M3%NICF_S(1)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG112A: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  LL = LL + 1
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+1) = CF%FN(JCF)
                  ICC=CF%CELL_LIST(2,LOW_IND,JCF); JCC=CF%CELL_LIST(3,LOW_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+2) = M%CUT_CELL(ICC)%H(JCC) ! H_LO
                  ICC=CF%CELL_LIST(2,HIGH_IND,JCF); JCC=CF%CELL_LIST(3,HIGH_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%H(JCC) ! H_HI
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112A
         ELSE
            PACK_REAL_SEND_PKG112A2: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  CF%FN_OMESH(JCF) = CF%FN(JCF)
                  ! No need to copy H_LO, H_HI
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112A2
         ENDIF
      ENDIF

      ! Exchange Velocity at end of PREDICTOR: To be used in RCEDGEs estimation of OMEGA and TAU at next substep.

      IF (CODE==3 .AND. M3%NICF_S(1)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG112A3: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  LL = LL + 1
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+1) = CF%VELS(JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+2) = CF%VEL_LNK(JCF)
                  ICC=CF%CELL_LIST(2,LOW_IND,JCF); JCC=CF%CELL_LIST(3,LOW_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%H(JCC) ! H_LO
                  ICC=CF%CELL_LIST(2,HIGH_IND,JCF); JCC=CF%CELL_LIST(3,HIGH_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+4) = M%CUT_CELL(ICC)%H(JCC) ! H_HI
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112A3
         ELSE
            PACK_REAL_SEND_PKG112A4: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  CF%VELS_OMESH(JCF)    = CF%VELS(JCF)
                  CF%VEL_LNK_OMESH(JCF) = CF%VEL_LNK(JCF)
                  ! No need to copy H_LO and H_HI.
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112A4
         ENDIF
      ENDIF
      IF (CODE==3 .AND. M3%NFCC_S(1)>0) THEN
         NQT2 = 1
         LL = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG121 : DO IFC=1,M3%NFCC_S(1)
               I     = M3%IIO_FC_S(IFC)
               J     = M3%JJO_FC_S(IFC)
               K     = M3%KKO_FC_S(IFC)
               X1AXIS= M3%AXS_FC_S(IFC)
               LL = LL + 1
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%US(I,J,K)                           ! U^* in x face I,J,K
               CASE(JAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%VS(I,J,K)                           ! V^* in y face I,J,K
               CASE(KAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%WS(I,J,K)                           ! W^* in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_SEND_PKG121
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG121: DO IFEP=1,M2%NFEP_R(1)
               ICF = M2%IFEP_R_1( LOW_IND,IFEP)
               INPE= M2%IFEP_R_1(HIGH_IND,IFEP)
               LL  = M1%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VELS_IND,INPE) = M%US(I,J,K)               ! U^* in x face I,J,K
               CASE(JAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VELS_IND,INPE) = M%VS(I,J,K)               ! V^* in y face I,J,K
               CASE(KAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VELS_IND,INPE) = M%WS(I,J,K)               ! W^* in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG121
            ! Second Loop cut-edges:
            PACK_REAL_RECV_PKG121E: DO IFEP=1,M2%NFEP_R(3)
               IEDGE= M2%IFEP_R_3( LOW_IND,IFEP)
               INPE = M2%IFEP_R_3(HIGH_IND,IFEP)
               LL   = M1%CC_RCEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%US(I,J,K)               ! U^* in x face I,J,K
               CASE(JAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%VS(I,J,K)               ! V^* in y face I,J,K
               CASE(KAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%WS(I,J,K)               ! W^* in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG121E
            PACK_REAL_RECV_PKG121EIB: DO IFEP=1,M2%NFEP_R(4)
               IEDGE= M2%IFEP_R_4( LOW_IND,IFEP)
               INPE = M2%IFEP_R_4(HIGH_IND,IFEP)
               LL   = M1%CC_IBEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%US(I,J,K)               ! U^* in x face I,J,K
               CASE(JAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%VS(I,J,K)               ! V^* in y face I,J,K
               CASE(KAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%WS(I,J,K)               ! W^* in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG121EIB
         ENDIF
      ENDIF

      ! Exchange of density and species mass fractions following the CORRECTOR update

      IF (CODE==4 .AND. M3%NICC_S(1)>0) THEN
         NQT2 = 4+N_TOTAL_SCALARS
         LL = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG111: DO ICC1=1,M3%NICC_S(1)
               ICC=M3%ICC_UNKZ_CC_S(ICC1)
               NCELL=M%CUT_CELL(ICC)%NCELL
               II1=M%CUT_CELL(ICC)%IJK(IAXIS)
               JJ1=M%CUT_CELL(ICC)%IJK(JAXIS)
               KK1=M%CUT_CELL(ICC)%IJK(KAXIS)
               DO JCC=1,NCELL
                  LL = LL + 1
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+1) = M%CUT_CELL(ICC)%RHO(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+2) = M%CUT_CELL(ICC)%TMP(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%RSUM(JCC)
                  M3%REAL_SEND_PKG11(NQT2*(LL-1)+4) = M%CUT_CELL(ICC)%DS(JCC)
                  DO NN=1,N_TOTAL_SCALARS
                     M3%REAL_SEND_PKG11(NQT2*(LL-1)+4+NN) = M%CUT_CELL(ICC)%ZZ(NN,JCC)
                  ENDDO
               ENDDO
            ENDDO PACK_REAL_SEND_PKG111
         ENDIF
      ENDIF

      ! Information for cell centered variables:
      IF (CODE==4 .AND. M3%NFCC_S(2)>0) THEN
         NQT2 = NQT2C+N_TOTAL_SCALARS
         LL = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG313 : DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+1) = M%H(I,J,K)                             ! Prev H in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%WS(I,J,K-1)+M%WS(I,J,K))     ! Wcen^* in I,J,K.
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+4) = M%RHO(I,J,K)                           ! RHO^n+1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+5) = M%TMP(I,J,K)                           ! TMP^n+1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+6) = M%RSUM(I,J,K)                          ! RSUM^n+1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+7) = M%MU(I,J,K)                            ! MU^*
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+8) = M%MU_DNS(I,J,K)                        ! MU_DNS^*
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+NQT2C)= M%RHOS(I,J,K)*(M%H(I,J,K)-M%KRES(I,J,K)) ! Previous substep pressure.
               DO NN=1,N_TOTAL_SCALARS
                  M3%REAL_SEND_PKG13(NQT2*(LL-1)+NQT2C+NN)= M%ZZ(I,J,K,NN)
               ENDDO
            ENDDO PACK_REAL_SEND_PKG313
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG313: DO IFEP=1,M2%NFEP_R(2)
               ICF = M2%IFEP_R_2( LOW_IND,IFEP)
               INPE= M2%IFEP_R_2(HIGH_IND,IFEP)
               LL  = M1%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_CC_S(LL)
               J     = M3%JJO_CC_S(LL)
               K     = M3%KKO_CC_S(LL)
               M1%CUT_FACE(ICF)%INT_CVARS(    INT_H_IND,INPE)= M%H(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(  INT_RHO_IND,INPE)= M%RHO(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(  INT_TMP_IND,INPE)= M%TMP(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS( INT_RSUM_IND,INPE)= M%RSUM(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(   INT_MU_IND,INPE)= M%MU(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(INT_MUDNS_IND,INPE)= M%MU_DNS(I,J,K)
               M1%CUT_FACE(ICF)%INT_CVARS(    INT_P_IND,INPE)= M%RHOS(I,J,K)*(M%H(I,J,K)-M%KRES(I,J,K))
               DO NN=1,N_TOTAL_SCALARS
                  M1%CUT_FACE(ICF)%INT_CVARS(INT_P_IND+NN,INPE)=M%ZZ(I,J,K,NN)
               ENDDO
            ENDDO PACK_REAL_RECV_PKG313
         ENDIF
      ENDIF


      ! Exchange velocity, momentum rhs and previous substep dH/Dx1 for cut-faces, in CORRECTOR, IBM forcing:
      IF (CODE==5 .AND. CORRECTOR .AND. M3%NICF_S(1)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG112B: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  LL = LL + 1
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+1) = CF%FN(JCF)
                  ICC=CF%CELL_LIST(2,LOW_IND,JCF); JCC=CF%CELL_LIST(3,LOW_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+2) = M%CUT_CELL(ICC)%HS(JCC) ! HS_LO
                  ICC=CF%CELL_LIST(2,HIGH_IND,JCF); JCC=CF%CELL_LIST(3,HIGH_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%HS(JCC) ! HS_HI
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112B
         ELSE
            PACK_REAL_SEND_PKG112B2: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  CF%FN_OMESH(JCF) = CF%FN(JCF)
                  ! No need to copy H_LO and H_HI.
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112B2
         ENDIF
      ENDIF

      ! Exchange Velocity and Pressure at end of CORRECTOR: To be used in RCEDGEs estimation of OMEGA and TAU next substep.

      IF (CODE==6 .AND. M3%NICF_S(1)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG112B3: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  LL = LL + 1
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+1) = CF%VEL(JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+2) = CF%VEL_LNK(JCF)
                  ICC=CF%CELL_LIST(2,LOW_IND,JCF); JCC=CF%CELL_LIST(3,LOW_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+3) = M%CUT_CELL(ICC)%HS(JCC) ! H_LO
                  ICC=CF%CELL_LIST(2,HIGH_IND,JCF); JCC=CF%CELL_LIST(3,HIGH_IND,JCF)
                  M3%REAL_SEND_PKG112(NQT2*(LL-1)+4) = M%CUT_CELL(ICC)%HS(JCC) ! H_HI
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112B3
         ELSE
            PACK_REAL_SEND_PKG112B4: DO ICF1=1,M3%NICF_S(1)
               ICF=M3%ICF_UFFB_CF_S(ICF1); CF => M%CUT_FACE(ICF)
               DO JCF=1,CF%NFACE
                  CF%VEL_OMESH(JCF)     = CF%VEL(JCF)
                  CF%VEL_LNK_OMESH(JCF) = CF%VEL_LNK(JCF)
                  ! No need to copy H_LO and H_HI.
               ENDDO
            ENDDO PACK_REAL_SEND_PKG112B4
         ENDIF
      ENDIF
      IF ((CODE==3 .OR. CODE==6) .AND. M3%NLKF_S>0) THEN
         IF (CODE==3) THEN
            UP => M%US ; VP => M%VS ; WP => M%WS
         ELSE
            UP => M%U  ; VP => M%V  ; WP => M%W
         ENDIF
         IF (RNODE/=SNODE) THEN
            LL = 4 * M3%NICF_S(2)
            DO KK=M3%K_MIN_S,M3%K_MAX_S
               DO JJ=M3%J_MIN_S,M3%J_MAX_S
                  DO II=M3%I_MIN_S,M3%I_MAX_S
                     ! U linked Velocity:
                     M3%REAL_SEND_PKG112(LL+1) = UP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS)>0) THEN ! Regular Face
                        M3%REAL_SEND_PKG112(LL+1) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,IAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,IAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS)>0) &
                        M3%REAL_SEND_PKG112(LL+1) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                     ! V linked Velocity:
                     M3%REAL_SEND_PKG112(LL+2) = VP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS)>0) THEN ! Regular Face
                        M3%REAL_SEND_PKG112(LL+2) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,JAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,JAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS)>0) &
                        M3%REAL_SEND_PKG112(LL+2) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                     ! W linked velocity:
                     M3%REAL_SEND_PKG112(LL+3) = WP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS)>0) THEN ! Regular Face
                        M3%REAL_SEND_PKG112(LL+3) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,KAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,KAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS)>0) &
                        M3%REAL_SEND_PKG112(LL+3) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                     LL = LL+3
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            UP2 => M2%U_LNK ; VP2 => M2%V_LNK ; WP2 => M2%W_LNK
            DO KK=M3%K_MIN_S,M3%K_MAX_S
               DO JJ=M3%J_MIN_S,M3%J_MAX_S
                  DO II=M3%I_MIN_S,M3%I_MAX_S
                     ! U linked Velocity:
                     UP2(II,JJ,KK) = UP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS)>0) THEN ! Regular Face
                        UP2(II,JJ,KK) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,IAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,IAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,IAXIS)>0) &
                        UP2(II,JJ,KK) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                     ! V linked Velocity:
                     VP2(II,JJ,KK) = VP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS)>0) THEN ! Regular Face
                        VP2(II,JJ,KK) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,JAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,JAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,JAXIS)>0) &
                        VP2(II,JJ,KK) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                     ! W linked velocity:
                     WP2(II,JJ,KK) = WP(II,JJ,KK)
                     IF (M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS)>0) THEN ! Regular Face
                        WP2(II,JJ,KK) = M%UN_LNK(M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS))
                     ELSEIF(M%FCVAR(II,JJ,KK,CC_IDRC,KAXIS)>0) THEN ! RC Face.
                        ICF=M%FCVAR(II,JJ,KK,CC_IDRC,KAXIS); IF(M%FCVAR(II,JJ,KK,CC_UNKF,KAXIS)>0) &
                        WP2(II,JJ,KK) = M%UN_LNK(M%RC_FACE(ICF)%UNKF)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF (CODE==6 .AND. M3%NFCC_S(1)>0) THEN
         NQT2 = 1
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG122 : DO IFC=1,M3%NFCC_S(1)
               I     = M3%IIO_FC_S(IFC)
               J     = M3%JJO_FC_S(IFC)
               K     = M3%KKO_FC_S(IFC)
               X1AXIS= M3%AXS_FC_S(IFC)
               LL = LL + 1
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%U(I,J,K)                            ! U^n+1 in x face I,J,K
               CASE(JAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%V(I,J,K)                            ! V^n+1 in y face I,J,K
               CASE(KAXIS)
                  M3%REAL_SEND_PKG12(NQT2*(LL-1)+1) = M%W(I,J,K)                            ! W^n+1 in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_SEND_PKG122
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG122: DO IFEP=1,M2%NFEP_R(1)
               ICF = M2%IFEP_R_1( LOW_IND,IFEP)
               INPE= M2%IFEP_R_1(HIGH_IND,IFEP)
               LL  = M1%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VEL_IND,INPE) = M%U(I,J,K)                ! U^n+1 in x face I,J,K
               CASE(JAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VEL_IND,INPE) = M%V(I,J,K)                ! V^n+1 in y face I,J,K
               CASE(KAXIS)
                  M1%CUT_FACE(ICF)%INT_FVARS( INT_VEL_IND,INPE) = M%W(I,J,K)                ! W^n+1 in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG122
            ! Second Loop cut-edges:
            PACK_REAL_RECV_PKG122E: DO IFEP=1,M2%NFEP_R(3)
               IEDGE= M2%IFEP_R_3( LOW_IND,IFEP)
               INPE = M2%IFEP_R_3(HIGH_IND,IFEP)
               LL   = M1%CC_RCEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%U(I,J,K)             ! U^n+1 in x face I,J,K
               CASE(JAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%V(I,J,K)             ! V^n+1 in y face I,J,K
               CASE(KAXIS)
                  M1%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%W(I,J,K)             ! W^n+1 in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG122E
            PACK_REAL_RECV_PKG122EIB: DO IFEP=1,M2%NFEP_R(4)
               IEDGE= M2%IFEP_R_4( LOW_IND,IFEP)
               INPE = M2%IFEP_R_4(HIGH_IND,IFEP)
               LL   = M1%CC_IBEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               I     = M3%IIO_FC_S(LL)
               J     = M3%JJO_FC_S(LL)
               K     = M3%KKO_FC_S(LL)
               X1AXIS= M3%AXS_FC_S(LL)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%U(I,J,K)             ! U^n+1 in x face I,J,K
               CASE(JAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%V(I,J,K)             ! V^n+1 in y face I,J,K
               CASE(KAXIS)
                  M1%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE) = M%W(I,J,K)             ! W^n+1 in z face I,J,K
               END SELECT
            ENDDO PACK_REAL_RECV_PKG122EIB
         ENDIF
      ENDIF

      ! Exchange H, RHO_0 and W velocity averaged to cell center, at PREDICTOR end of step:

      IF (CODE==3 .AND. M3%NFCC_S(2)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG13 : DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+1) = M%H(I,J,K)                             ! H^n in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%WS(I,J,K-1)+M%WS(I,J,K))     ! Wcen^* in I,J,K.
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+4) = M%MU(I,J,K)                            ! MU in I,J,K
            ENDDO PACK_REAL_SEND_PKG13
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG13: DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+1) = M%H(I,J,K)                             ! H^n in cell I,J,K
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%WS(I,J,K-1)+M%WS(I,J,K))     ! Wcen^* in I,J,K.
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+4) = M%MU(I,J,K)                            ! MU in I,J,K
            ENDDO PACK_REAL_RECV_PKG13
         ENDIF
      ENDIF

      ! Exchange H, RHO_0 and W velocity averaged to cell center, at CORRECTOR end of step:

      IF (CODE==6 .AND. M3%NFCC_S(2)>0) THEN
         NQT2 = 4
         LL   = 0
         IF (RNODE/=SNODE) THEN
            PACK_REAL_SEND_PKG113 : DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+1) = M%HS(I,J,K)                            ! H^* in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%W(I,J,K-1)+M%W(I,J,K))       ! Wcen  in I,J,K.
               M3%REAL_SEND_PKG13(NQT2*(LL-1)+4) = M%MU(I,J,K)                            ! MU in I,J,K
            ENDDO PACK_REAL_SEND_PKG113
         ELSE
            M2=>MESHES(NOM)%OMESH(NM)
            PACK_REAL_RECV_PKG113: DO ICC=1,M3%NFCC_S(2)
               I     = M3%IIO_CC_S(ICC)
               J     = M3%JJO_CC_S(ICC)
               K     = M3%KKO_CC_S(ICC)
               LL = LL + 1
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+1) = M%HS(I,J,K)                            ! H^* in cell I,J,K
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+2) = M%RHO_0(K)                             ! RHO_0 in cell I,J,K
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+3) = 0.5_EB*(M%W(I,J,K-1)+M%W(I,J,K))       ! Wcen  in I,J,K.
               M2%REAL_RECV_PKG13(NQT2*(LL-1)+4) = M%MU(I,J,K)                            ! MU in I,J,K
            ENDDO PACK_REAL_RECV_PKG113
         ENDIF
      ENDIF

   ENDDO RECEIVING_MESH_LOOP_2
ENDDO SENDING_MESH_LOOP_2

! Exchange Scalars:
IF (N_MPI_PROCESSES>1 .AND. (CODE==1.OR.CODE==4) .AND. N_REQ11>0) THEN
   CALL MPI_STARTALL(N_REQ11,REQ11(1:N_REQ11),IERR)
   CALL CC_TIMEOUT('REQ11',N_REQ11,REQ11(1:N_REQ11))
ENDIF

! Exchange FVX,FVY,FVZ data for gas cut-faces:
IF (N_MPI_PROCESSES>1 .AND. ANY(CODE==(/3,5,6/)) .AND. N_REQ112>0) THEN
   CALL MPI_STARTALL(N_REQ112,REQ112(1:N_REQ112),IERR)
   CALL CC_TIMEOUT('REQ112',N_REQ112,REQ112(1:N_REQ112))
ENDIF

! Exchange End of Step velocity data for RCEDGEs:
IF (N_MPI_PROCESSES>1 .AND. (CODE==3.OR.CODE==6) .AND. N_REQ12>0) THEN
   CALL MPI_STARTALL(N_REQ12,REQ12(1:N_REQ12),IERR)
   CALL CC_TIMEOUT('REQ12',N_REQ12,REQ12(1:N_REQ12))
ENDIF

! Exchange scalar data End of step cell-centered data:
IF (N_MPI_PROCESSES>1 .AND. (CODE==1.OR.CODE==4.OR.CODE==3.OR.CODE==6) .AND. N_REQ13>0) THEN
   CALL MPI_STARTALL(N_REQ13,REQ13(1:N_REQ13),IERR)
   CALL CC_TIMEOUT('REQ13',N_REQ13,REQ13(1:N_REQ13))
ENDIF

! Receive the information sent above into the appropriate arrays.

RECV_MESH_LOOP: DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M =>MESHES(NOM)

   SEND_MESH_LOOP: DO NM=1,NMESHES

      M2=>MESHES(NOM)%OMESH(NM)

      RNODE = PROCESS(NOM)
      SNODE = PROCESS(NM)

      RNODE_SNODE_IF: IF (RNODE/=SNODE) THEN

         ! Unpack densities and species mass fractions following PREDICTOR exchange

         IF (CODE==1 .AND. M2%NICC_R(1)>0) THEN
            NQT2 = 4+N_TOTAL_SCALARS
            LL   = 0
            ! Copy-cut cell scalar quantities from MESHES(NOM)%OMESH(NM) cells to MESHES(NM) (i.e. other mesh) cut-cells:
            ! Use External wall cell loop:
            EXTERNAL_WALL_LOOP_1 : DO IW=1,M%N_EXTERNAL_WALL_CELLS
               WC=>M%WALL(IW)
               EWC=>M%EXTERNAL_WALL(IW)
               BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
               IF (.NOT.(WC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY)) CYCLE EXTERNAL_WALL_LOOP_1
               IF (EWC%NOM/=NM) CYCLE EXTERNAL_WALL_LOOP_1
               IF (M%CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1
               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                       ICC   = MESHES(NM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
                       IF (ICC > 0) THEN
                          DO JCC=1,MESHES(NM)%CUT_CELL(ICC)%NCELL
                             LL = LL + 1
                             MESHES(NM)%CUT_CELL(ICC)%RHOS(JCC) = M2%REAL_RECV_PKG11(NQT2*(LL-1)+1)
                             MESHES(NM)%CUT_CELL(ICC)%TMP(JCC)  = M2%REAL_RECV_PKG11(NQT2*(LL-1)+2)
                             MESHES(NM)%CUT_CELL(ICC)%RSUM(JCC) = M2%REAL_RECV_PKG11(NQT2*(LL-1)+3)
                             MESHES(NM)%CUT_CELL(ICC)%D(JCC)    = M2%REAL_RECV_PKG11(NQT2*(LL-1)+4)
                             DO NN=1,N_TOTAL_SCALARS
                                MESHES(NM)%CUT_CELL(ICC)%ZZS(NN,JCC) = M2%REAL_RECV_PKG11(NQT2*(LL-1)+4+NN)
                             ENDDO
                          ENDDO
                       ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO EXTERNAL_WALL_LOOP_1
         ENDIF

         IF((CODE==1 .OR. CODE==4) .AND. M2%NFCC_R(2)>0) THEN
            NQT2 = NQT2C+N_TOTAL_SCALARS
            DO IFEP=1,M2%NFEP_R(2)
               ICF = M2%IFEP_R_2( LOW_IND,IFEP)
               INPE= M2%IFEP_R_2(HIGH_IND,IFEP)
               LL  = M%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               M%CUT_FACE(ICF)%INT_CVARS(    INT_H_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+1)
               M%CUT_FACE(ICF)%INT_CVARS(  INT_RHO_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+4)
               M%CUT_FACE(ICF)%INT_CVARS(  INT_TMP_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+5)
               M%CUT_FACE(ICF)%INT_CVARS( INT_RSUM_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+6)
               M%CUT_FACE(ICF)%INT_CVARS(   INT_MU_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+7)
               M%CUT_FACE(ICF)%INT_CVARS(INT_MUDNS_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+8)
               M%CUT_FACE(ICF)%INT_CVARS(    INT_P_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+NQT2C)
               DO NN=1,N_TOTAL_SCALARS
                  M%CUT_FACE(ICF)%INT_CVARS(INT_P_IND+NN,INPE)=M2%REAL_RECV_PKG13(NQT2*(LL-1)+NQT2C+NN)
               ENDDO
            ENDDO
         ENDIF

         ! Unpack velocity, momentum rhs and previous substep dH/Dx1 for cut-faces, in PREDICTOR or CORRECTOR, IBM forcing:
         IF (CODE==5  .AND. M2%NICF_R(1)>0) THEN
            NQT2 = 4
            LL   = 0
            DO ICF1=1,M2%NICF_R(1)
               ICF = M2%ICF_UFFB_CF_R(ICF1)
               IF (ICF > 0) THEN
                  CF => MESHES(NM)%CUT_FACE(ICF)
                  DO JCF=1,CF%NFACE
                     LL = LL + 1
                     CF%FN_OMESH(JCF)   = M2%REAL_RECV_PKG112(NQT2*(LL-1)+1)
                     ICC =CF%CELL_LIST(2, LOW_IND,JCF); JCC =CF%CELL_LIST(3, LOW_IND,JCF)
                     ICC1=CF%CELL_LIST(2,HIGH_IND,JCF); JCC1=CF%CELL_LIST(3,HIGH_IND,JCF)
                     IF(PREDICTOR) THEN
                        MESHES(NM)%CUT_CELL(ICC )%H(JCC )  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+2)
                        MESHES(NM)%CUT_CELL(ICC1)%H(JCC1)  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+3)
                     ELSE
                        MESHES(NM)%CUT_CELL(ICC )%HS(JCC ) = M2%REAL_RECV_PKG112(NQT2*(LL-1)+2)
                        MESHES(NM)%CUT_CELL(ICC1)%HS(JCC1) = M2%REAL_RECV_PKG112(NQT2*(LL-1)+3)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         IF (CODE==3 .AND. M2%NFCC_R(1)>0) THEN
            NQT2 = 1
            ! First loop cut-faces:
            DO IFEP=1,M2%NFEP_R(1) ! Gasphase and Boundary cut-faces:
               ICF = M2%IFEP_R_1( LOW_IND,IFEP);
               INPE= M2%IFEP_R_1(HIGH_IND,IFEP)
               LL  = M%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               M%CUT_FACE(ICF)%INT_FVARS(INT_VELS_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^*
            ENDDO
            ! Second Loop cut-edges:
            DO IFEP=1,M2%NFEP_R(3)
               IEDGE= M2%IFEP_R_3( LOW_IND,IFEP)
               INPE = M2%IFEP_R_3(HIGH_IND,IFEP)
               LL   = M%CC_RCEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               M%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^*, added to INT_VEL_IND pos.
            ENDDO
            DO IFEP=1,M2%NFEP_R(4)
               IEDGE= M2%IFEP_R_4( LOW_IND,IFEP)
               INPE = M2%IFEP_R_4(HIGH_IND,IFEP)
               LL   = M%CC_IBEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               M%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^*, added to INT_VEL_IND pos.
            ENDDO
         ENDIF
         IF(CODE==3 .AND. M2%NICF_R(1)>0) THEN
            NQT2 = 4
            LL   = 0
            DO ICF1=1,M2%NICF_R(1)
               ICF = M2%ICF_UFFB_CF_R(ICF1)
               IF (ICF > 0) THEN
                  CF => MESHES(NM)%CUT_FACE(ICF)
                  DO JCF=1,CF%NFACE
                     LL = LL + 1
                     CF%VELS_OMESH(JCF)   = M2%REAL_RECV_PKG112(NQT2*(LL-1)+1)
                     CF%VEL_LNK_OMESH(JCF)= M2%REAL_RECV_PKG112(NQT2*(LL-1)+2)
                     ICC =CF%CELL_LIST(2, LOW_IND,JCF); JCC =CF%CELL_LIST(3, LOW_IND,JCF)
                     ICC1=CF%CELL_LIST(2,HIGH_IND,JCF); JCC1=CF%CELL_LIST(3,HIGH_IND,JCF)
                     MESHES(NM)%CUT_CELL(ICC )%H(JCC )  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+3)
                     MESHES(NM)%CUT_CELL(ICC1)%H(JCC1)  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+4)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         ! Unpack densities and species mass fractions following CORRECTOR exchange

         IF (CODE==4 .AND. M2%NICC_R(1)>0) THEN
            NQT2 = 4+N_TOTAL_SCALARS
            LL   = 0
            ! Copy-cut cell scalar quantities from MESHES(NOM)%OMESH(NM) cells to MESHES(NM) (i.e. other mesh) cut-cells:
            ! Use External wall cell loop:
            EXTERNAL_WALL_LOOP_2 : DO IW=1,M%N_EXTERNAL_WALL_CELLS
               WC=>M%WALL(IW)
               IF (.NOT.(WC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY)) CYCLE EXTERNAL_WALL_LOOP_2
               BC=>M%BOUNDARY_COORD(WC%BC_INDEX)
               EWC=>M%EXTERNAL_WALL(IW)
               IF (EWC%NOM/=NM) CYCLE EXTERNAL_WALL_LOOP_2
               IF (M%CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) &
                  CYCLE EXTERNAL_WALL_LOOP_2
               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                       ICC   = MESHES(NM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
                       IF (ICC > 0) THEN
                          DO JCC=1,MESHES(NM)%CUT_CELL(ICC)%NCELL
                             LL = LL + 1
                             MESHES(NM)%CUT_CELL(ICC)%RHO(JCC)  = M2%REAL_RECV_PKG11(NQT2*(LL-1)+1)
                             MESHES(NM)%CUT_CELL(ICC)%TMP(JCC)  = M2%REAL_RECV_PKG11(NQT2*(LL-1)+2)
                             MESHES(NM)%CUT_CELL(ICC)%RSUM(JCC) = M2%REAL_RECV_PKG11(NQT2*(LL-1)+3)
                             MESHES(NM)%CUT_CELL(ICC)%DS(JCC)   = M2%REAL_RECV_PKG11(NQT2*(LL-1)+4)
                             DO NN=1,N_TOTAL_SCALARS
                                MESHES(NM)%CUT_CELL(ICC)%ZZ(NN,JCC) = M2%REAL_RECV_PKG11(NQT2*(LL-1)+4+NN)
                             ENDDO
                          ENDDO
                       ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO EXTERNAL_WALL_LOOP_2
         ENDIF

         IF (CODE==6 .AND. M2%NFCC_R(1)>0) THEN
            NQT2 = 1
            ! First loop cut-faces:
            DO IFEP=1,M2%NFEP_R(1) ! Gasphase and Boundary cut-faces:
               ICF = M2%IFEP_R_1( LOW_IND,IFEP);
               INPE= M2%IFEP_R_1(HIGH_IND,IFEP)
               LL  = M%CUT_FACE(ICF)%INT_NOMIND(HIGH_IND,INPE)
               M%CUT_FACE(ICF)%INT_FVARS(INT_VEL_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^n+1
            ENDDO
            ! Second Loop cut-edges:
            DO IFEP=1,M2%NFEP_R(3)
               IEDGE= M2%IFEP_R_3( LOW_IND,IFEP)
               INPE = M2%IFEP_R_3(HIGH_IND,IFEP)
               LL   = M%CC_RCEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               M%CC_RCEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^n+1
            ENDDO
            DO IFEP=1,M2%NFEP_R(4)
               IEDGE= M2%IFEP_R_4( LOW_IND,IFEP)
               INPE = M2%IFEP_R_4(HIGH_IND,IFEP)
               LL   = M%CC_IBEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
               M%CC_IBEDGE(IEDGE)%INT_FVARS(INT_VEL_IND,INPE)= M2%REAL_RECV_PKG12(NQT2*(LL-1)+1) ! Vel^n+1
            ENDDO
         ENDIF
         IF(CODE==6 .AND. M2%NICF_R(1)>0) THEN
            NQT2 = 4
            LL   = 0
            DO ICF1=1,M2%NICF_R(1)
               ICF = M2%ICF_UFFB_CF_R(ICF1)
               IF (ICF > 0) THEN
                  CF => MESHES(NM)%CUT_FACE(ICF)
                  DO JCF=1,CF%NFACE
                     LL = LL + 1
                     CF%VEL_OMESH(JCF)    = M2%REAL_RECV_PKG112(NQT2*(LL-1)+1)
                     CF%VEL_LNK_OMESH(JCF)= M2%REAL_RECV_PKG112(NQT2*(LL-1)+2)
                     ICC =CF%CELL_LIST(2, LOW_IND,JCF); JCC =CF%CELL_LIST(3, LOW_IND,JCF)
                     ICC1=CF%CELL_LIST(2,HIGH_IND,JCF); JCC1=CF%CELL_LIST(3,HIGH_IND,JCF)
                     MESHES(NM)%CUT_CELL(ICC )%HS(JCC )  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+3)
                     MESHES(NM)%CUT_CELL(ICC1)%HS(JCC1)  = M2%REAL_RECV_PKG112(NQT2*(LL-1)+4)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         IF ((CODE==3.OR.CODE==6) .AND. M2%NLKF_R>0) THEN
            UP2 => M2%U_LNK ; VP2 => M2%V_LNK ; WP2 => M2%W_LNK
            LL = 4 * M2%NICF_R(2)
            DO KK=M2%K_MIN_R,M2%K_MAX_R
               DO JJ=M2%J_MIN_R,M2%J_MAX_R
                  DO II=M2%I_MIN_R,M2%I_MAX_R
                     UP2(II,JJ,KK) = M2%REAL_RECV_PKG112(LL+1)
                     VP2(II,JJ,KK) = M2%REAL_RECV_PKG112(LL+2)
                     WP2(II,JJ,KK) = M2%REAL_RECV_PKG112(LL+3)
                     LL = LL+3
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ENDIF RNODE_SNODE_IF

      ! Unpack H, RHO_0 and W velocity averaged to cell center, at PREDICTOR or CORRECTOR end of step:

      IF ( (CODE==3 .OR. CODE==6) .AND. M2%NFCC_R(2)>0) THEN
         NQT2 = 4
         ! First loop cut-cells:
         VIND = 0
         DO ICC=1,M%N_CUTCELL_MESH
            DO ICELL=0,M%CUT_CELL(ICC)%NCELL
               DO EP=1,INT_N_EXT_PTS  ! External point for cell ICELL
                  INT_NPE_LO = M%CUT_CELL(ICC)%INT_NPE(LOW_IND,VIND,EP,ICELL)
                  INT_NPE_HI = M%CUT_CELL(ICC)%INT_NPE(HIGH_IND,VIND,EP,ICELL)
                  DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                     IF (M%CUT_CELL(ICC)%INT_NOMIND( LOW_IND,INPE) /= NM) CYCLE
                     LL     = M%CUT_CELL(ICC)%INT_NOMIND(HIGH_IND,INPE)
                     M%CUT_CELL(ICC)%INT_CCVARS(   INT_H_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+1) ! H^n, or H^s
                     M%CUT_CELL(ICC)%INT_CCVARS(INT_RHO0_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+2) ! RHO_0
                     M%CUT_CELL(ICC)%INT_CCVARS(INT_WCEN_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+3) ! Wcen^*, or Wcen^n+1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         ! Then Loop IBEDGES:
         DO IFEP=1,M2%NFEP_R(5)
            IEDGE= M2%IFEP_R_5( LOW_IND,IFEP)
            INPE = M2%IFEP_R_5(HIGH_IND,IFEP)
            LL   = M%CC_IBEDGE(IEDGE)%INT_NOMIND(HIGH_IND,INPE)
            M%CC_IBEDGE(IEDGE)%INT_CVARS(INT_MU_IND,INPE)= M2%REAL_RECV_PKG13(NQT2*(LL-1)+4) ! MU.
         ENDDO
      ENDIF

   ENDDO SEND_MESH_LOOP
ENDDO RECV_MESH_LOOP

IF(CODE==1 .OR. CODE==4) THEN
  CALL FILL_GCCUTCELL_SPECIES
ELSEIF(CODE==3 .OR. CODE==6) THEN
   CALL CC_H_INTERP
   CALL CC_RHO0W_INTERP
ENDIF

RETURN

CONTAINS

SUBROUTINE FILL_GCCUTCELL_SPECIES


REAL(EB):: PRFCT, VCELL, RHO_CC, TMP_CC, RSUM_CC, D_CC, ZZ_CC(1:N_TOTAL_SCALARS), VOL
TYPE (OMESH_TYPE), POINTER :: OM
TYPE(CC_CUTCELL_TYPE), POINTER :: OCC
INTEGER :: NM,NOM,NN,ICC,JCC,IW,IIO,JJO,KKO

! Here inject OMESH cut-cell info obtained in MESH_CC_EXCHANGE into ghost-cell cc containers:
PRFCT = 0._EB; IF (PREDICTOR) PRFCT = 1._EB
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   EXTERNAL_WALL_LOOP : DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      IF (.NOT.(WC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY)) CYCLE EXTERNAL_WALL_LOOP
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      EWC=>EXTERNAL_WALL(IW)
      IF (CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP
      ! Do volume average to a cell container for ghost cell II,JJ,KK:
      NOM = EWC%NOM
      OM  => MESHES(NM)%OMESH(NOM)
      RHO_CC=0._EB; TMP_CC=0._EB; RSUM_CC=0._EB; D_CC=0._EB; ZZ_CC(1:N_TOTAL_SCALARS)=0._EB; VOL=0._EB
      DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
         DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
            DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
              IF (MESHES(NOM)%CELL(MESHES(NOM)%CELL_INDEX(IIO,JJO,KKO))%SOLID) CYCLE
              ICC   = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
              IF (ICC > 0) THEN ! Cut-cells:
                 OCC => MESHES(NOM)%CUT_CELL(ICC)
                 DO JCC=1,OCC%NCELL
                    VCELL   = OCC%VOLUME(JCC)
                    RHO_CC  = RHO_CC  + (PRFCT *OCC%RHOS(JCC) + (1._EB-PRFCT)*OCC%RHO(JCC))*VCELL
                    TMP_CC  = TMP_CC  + OCC%TMP(JCC)*VCELL
                    RSUM_CC = RSUM_CC + OCC%RSUM(JCC)*VCELL
                    D_CC    = D_CC    + (PRFCT *OCC%D(JCC) + (1._EB-PRFCT)*OCC%DS(JCC))*VCELL
                    DO NN=1,N_TOTAL_SCALARS
                       ZZ_CC(NN) = ZZ_CC(NN) + (PRFCT *OCC%ZZS(NN,JCC) + (1._EB-PRFCT)*OCC%ZZ(NN,JCC))*VCELL
                    ENDDO
                    VOL     = VOL  + VCELL
                 ENDDO
              ELSEIF(MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_CGSC) == CC_GASPHASE) THEN ! Regular cell:
                 VCELL = MESHES(NOM)%DX(IIO)*MESHES(NOM)%DY(JJO)*MESHES(NOM)%DZ(KKO)
                 RHO_CC  = RHO_CC + (PRFCT*RHOS(BC%II,BC%JJ,BC%KK) + (1._EB-PRFCT)*RHO(BC%II,BC%JJ,BC%KK))*VCELL
                 TMP_CC  = TMP_CC + TMP(BC%II,BC%JJ,BC%KK)*VCELL
                 D_CC    = D_CC   + (PRFCT*D(BC%II,BC%JJ,BC%KK) + (1._EB-PRFCT)*DS(BC%II,BC%JJ,BC%KK))*VCELL
                 DO NN=1,N_TOTAL_SCALARS
                 ZZ_CC(NN)=ZZ_CC(NN)+(PRFCT*ZZS(BC%II,BC%JJ,BC%KK,NN)+(1._EB-PRFCT)*ZZ(BC%II,BC%JJ,BC%KK,NN))*VCELL
                 ENDDO
                 VOL     = VOL + VCELL
              ENDIF
            ENDDO
         ENDDO
      ENDDO
      ! Add volume averaged variables into ghost cut-cell:
      ICC   = CCVAR(BC%II,BC%JJ,BC%KK,CC_IDCC)
      IF (PREDICTOR) THEN
         DO JCC=1,CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%RHOS(JCC) = RHO_CC/VOL
            CUT_CELL(ICC)%TMP(JCC)  = TMP_CC/VOL
            CUT_CELL(ICC)%RSUM(JCC) = RSUM_CC/VOL
            CUT_CELL(ICC)%D(JCC)    = D_CC/VOL
            CUT_CELL(ICC)%ZZS(1:N_TOTAL_SCALARS,JCC) = ZZ_CC(1:N_TOTAL_SCALARS)/VOL
         ENDDO
      ELSE
         DO JCC=1,CUT_CELL(ICC)%NCELL
            CUT_CELL(ICC)%RHO(JCC)   = RHO_CC/VOL
            CUT_CELL(ICC)%TMP(JCC)   = TMP_CC/VOL
            !CUT_CELL(ICC)%RSUM(JCC) = RSUM_CC/VOL
            CUT_CELL(ICC)%DS(JCC)    = D_CC/VOL
            CUT_CELL(ICC)%ZZ(1:N_TOTAL_SCALARS,JCC) = ZZ_CC(1:N_TOTAL_SCALARS)/VOL
         ENDDO
     ENDIF
   ENDDO EXTERNAL_WALL_LOOP
ENDDO MESH_LOOP


END SUBROUTINE FILL_GCCUTCELL_SPECIES


SUBROUTINE CC_TIMEOUT(RNAME,NR,RR)

REAL(EB) :: START_TIME,WAIT_TIME
INTEGER :: NR
TYPE (MPI_REQUEST), DIMENSION(:) :: RR
LOGICAL :: FLAG
CHARACTER(*) :: RNAME

IF (.NOT.PROFILING) THEN

   START_TIME = MPI_WTIME()
   FLAG = .FALSE.
   DO WHILE(.NOT.FLAG)
      CALL MPI_TESTALL(NR,RR(1:NR),FLAG,MPI_STATUSES_IGNORE,IERR)
      WAIT_TIME = MPI_WTIME() - START_TIME
      IF (WAIT_TIME>MPI_TIMEOUT) THEN
         WRITE(LU_ERR,'(A,A,A,I6,A,A)') 'CC_TIMEOUT Error: ',TRIM(RNAME),' timed out for MPI process ',MY_RANK
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      ENDIF
   ENDDO
ELSE

   CALL MPI_WAITALL(NR,RR(1:NR),MPI_STATUSES_IGNORE,IERR)

ENDIF

END SUBROUTINE CC_TIMEOUT

END SUBROUTINE MESH_CC_EXCHANGE

END MODULE CC_EXCHANGE
