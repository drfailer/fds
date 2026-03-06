!  +++++++++++++++++++++++ CC_INIT ++++++++++++++++++++++++++

! Initialization, setup, and finalization routines for the
! cut-cell / immersed-boundary method.

MODULE CC_INIT

USE CC_SCALARS_DATA
USE CC_SCALARS, ONLY: GET_LINKED_VELOCITIES, CC_VELOCITY_CUTFACES, &
                      CC_H_INTERP, CC_RHO0W_INTERP, &
                      CC_VELOCITY_BC, CC_CHECK_DIVERGENCE
USE CC_EXCHANGE, ONLY: MESH_CC_EXCHANGE
USE COMPLEX_GEOMETRY
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: CURRENT_TIME, GET_FILE_NUMBER
USE MATH_FUNCTIONS, ONLY: GET_SCALAR_FACE_VALUE

IMPLICIT NONE (TYPE,EXTERNAL)

PRIVATE

PUBLIC :: CC_SET_DATA, INIT_CUTCELL_DATA, FINISH_CC, &
          GET_CRTCFCC_INT_STENCILS

CONTAINS


! --------------------------- CC_EXCHANGE_UNPACKING_ARRAYS --------------------------

SUBROUTINE CC_EXCHANGE_UNPACKING_ARRAYS()

! Local Variables:
INTEGER :: NM,NOM,NOOM,IFEP,ICD_SGN
TYPE (MESH_TYPE), POINTER :: M
TYPE (OMESH_TYPE), POINTER :: M2
INTEGER :: EP,INPE,INT_NPE_LO,INT_NPE_HI,VIND,IEDGE

RECV_MESH_LOOP: DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M =>MESHES(NOM)

   SEND_MESH_LOOP: DO NM=1,NMESHES

      M2=>MESHES(NOM)%OMESH(NM)

      ! Boundary and gasphase cut-faces and rcedges, face centered variables for interpolation:
      CF_FC_IF : IF(M2%NFCC_R(1)>0) THEN
         ! RCEDGES:
         ! Count:
         DO IEDGE=1,M%CC_NRCEDGE
            DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
               DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
                  INT_NPE_LO = M%CC_RCEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,0)
                  INT_NPE_HI = M%CC_RCEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0)
                  DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                     NOOM   = M%CC_RCEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                     M2%NFEP_R(3) = M2%NFEP_R(3) + 1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF (M2%NFEP_R(3) > 0) THEN
            ! Allocate:
            IF (ALLOCATED(M2%IFEP_R_3)) DEALLOCATE(M2%IFEP_R_3)
            ALLOCATE(M2%IFEP_R_3(LOW_IND:HIGH_IND,M2%NFEP_R(3))); M2%IFEP_R_3 = CC_UNDEFINED
            ! Add index entries:
            IFEP = 0
            DO IEDGE=1,M%CC_NRCEDGE
               DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
                  DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
                     INT_NPE_LO = M%CC_RCEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,0)
                     INT_NPE_HI = M%CC_RCEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0)
                     DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                        NOOM   = M%CC_RCEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                        IFEP = IFEP + 1
                        M2%IFEP_R_3( LOW_IND:HIGH_IND,IFEP) = (/ IEDGE, INPE /)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         ! Then IBEDGES:
         ! Count:
         DO IEDGE=1,M%CC_NIBEDGE
            DO ICD_SGN=-2,2
               IF(ICD_SGN==0) CYCLE
               DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
                  DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
                    INT_NPE_LO = M%CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
                    INT_NPE_HI = M%CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN)
                    DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                       NOOM   = M%CC_IBEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                       M2%NFEP_R(4) = M2%NFEP_R(4) + 1
                    ENDDO
                 ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF (M2%NFEP_R(4) > 0) THEN
            ! Allocate:
            IF (ALLOCATED(M2%IFEP_R_4)) DEALLOCATE(M2%IFEP_R_4)
            ALLOCATE(M2%IFEP_R_4(LOW_IND:HIGH_IND,M2%NFEP_R(4))); M2%IFEP_R_4 = CC_UNDEFINED
            ! Add index entries:
            IFEP = 0
            DO IEDGE=1,M%CC_NIBEDGE
               DO ICD_SGN=-2,2
                  IF(ICD_SGN==0) CYCLE
                  DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
                     DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
                        INT_NPE_LO = M%CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
                        INT_NPE_HI = M%CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN)
                        DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                           NOOM   = M%CC_IBEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                           IFEP = IFEP + 1
                           M2%IFEP_R_4( LOW_IND:HIGH_IND,IFEP) = (/ IEDGE, INPE /)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF CF_FC_IF

      ! Boundary cut-faces, cell centered variables for interpolation:
      BNDCF_CC_IF : IF(M2%NFCC_R(2)>0) THEN
         VIND = 0 ! Cell centered variables.
         ! Case of IBEDGES:
         ! Count:
         DO IEDGE=1,M%CC_NIBEDGE
            DO ICD_SGN=-2,2
               IF(ICD_SGN==0) CYCLE
               DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
                  INT_NPE_LO = M%CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
                  INT_NPE_HI = M%CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN)
                  DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                     NOOM   = M%CC_IBEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                     M2%NFEP_R(5) = M2%NFEP_R(5) + 1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF (M2%NFEP_R(5)>0) THEN
            ! Allocate:
            IF (ALLOCATED(M2%IFEP_R_5)) DEALLOCATE(M2%IFEP_R_5)
            ALLOCATE(M2%IFEP_R_5(LOW_IND:HIGH_IND,M2%NFEP_R(5))); M2%IFEP_R_5 = CC_UNDEFINED
            ! Add index entries:
            IFEP = 0
            DO IEDGE=1,M%CC_NIBEDGE
               DO ICD_SGN=-2,2
                  IF(ICD_SGN==0) CYCLE
                  DO EP=1,INT_N_EXT_PTS  ! External point for face IEDGE
                     INT_NPE_LO = M%CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
                     INT_NPE_HI = M%CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN)
                     DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                        NOOM   = M%CC_IBEDGE(IEDGE)%INT_NOMIND( LOW_IND,INPE); IF (NOOM /= NM) CYCLE
                        IFEP = IFEP + 1
                        M2%IFEP_R_5( LOW_IND:HIGH_IND,IFEP) = (/ IEDGE, INPE /)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ENDIF BNDCF_CC_IF
   ENDDO SEND_MESH_LOOP
ENDDO RECV_MESH_LOOP

RETURN
END SUBROUTINE CC_EXCHANGE_UNPACKING_ARRAYS



! -------------------------------- CC_SET_DATA ----------------------------------

SUBROUTINE CC_SET_DATA(FIRST_CALL)

USE MEMORY_FUNCTIONS, ONLY: EXCHANGE_GEOMETRY_INFO

LOGICAL, INTENT(IN) :: FIRST_CALL

! Local Variables:
INTEGER :: NM,ICALL
REAL(EB):: LX,LY,LZ,MAX_DIST
REAL(EB):: TNOW,TNOW2,TDEL,MIN_XS(1:3),MAX_XF(1:3)

INTEGER :: ICF, IG
CHARACTER(80) :: FN_CCTIME
CHARACTER(200)::TCFORM

INTEGER :: ICC,JCC,I,J,K,ICC2,JCC2,JCF,IFACE,FTYPE,IFC2,IFACE2,ICF2,JCF2
REAL(EB):: ACRT
TYPE(CC_CUTCELL_TYPE), POINTER :: CC2
TYPE(CC_CUTFACE_TYPE), POINTER :: CF2

TNOW2 = CURRENT_TIME()

SET_CUTCELLS_CALL_IF : IF(FIRST_CALL) THEN

! Plane by plane Evaluation of stesses for IBEDGES, a la OBSTS.
CC_ONLY_IBEDGES_FLAG=.FALSE.
THRES_FCT_EP = -1._EB

IF (N_GEOMETRY==0 .AND. .NOT.(PERIODIC_TEST==103 .OR. PERIODIC_TEST==11 .OR. PERIODIC_TEST==7)) THEN
   IF (MY_RANK==0) THEN
      WRITE(LU_ERR,*) ' '
      WRITE(LU_ERR,*) 'CCIBM Setup Error : &MISC CC_IBM=.TRUE., but no &GEOM namelist defined on input file.'
      WRITE(LU_ERR,*) ' '
   ENDIF
   STOP_STATUS = SETUP_STOP
   RETURN
ENDIF

! Defined relative GEOMEPS:
! Find largest domain distance to define relative epsilon:
MIN_XS(1:3) = (/ MESHES(1)%XS, MESHES(1)%YS, MESHES(1)%ZS /)
MAX_XF(1:3) = (/ MESHES(1)%XF, MESHES(1)%YF, MESHES(1)%ZF /)
DO NM=2,NMESHES
   MIN_XS(1) = MIN(0._EB,MIN_XS(1),MESHES(NM)%XS)
   MIN_XS(2) = MIN(0._EB,MIN_XS(2),MESHES(NM)%YS)
   MIN_XS(3) = MIN(0._EB,MIN_XS(3),MESHES(NM)%ZS)
   MAX_XF(1) = MAX(0._EB,MAX_XF(1),MESHES(NM)%XF)
   MAX_XF(2) = MAX(0._EB,MAX_XF(2),MESHES(NM)%YF)
   MAX_XF(3) = MAX(0._EB,MAX_XF(3),MESHES(NM)%ZF)
ENDDO
LX = MAX_XF(1) - MIN_XS(1)
LY = MAX_XF(2) - MIN_XS(2)
LZ = MAX_XF(3) - MIN_XS(3)
MAX_DIST=MAX(LX,LY,LZ)

! Now test against GEOMETRY size:
DO IG=1,N_GEOMETRY
   LX = GEOMETRY(IG)%GEOM_BOX(HIGH_IND,IAXIS) - GEOMETRY(IG)%GEOM_BOX( LOW_IND,IAXIS)
   LY = GEOMETRY(IG)%GEOM_BOX(HIGH_IND,JAXIS) - GEOMETRY(IG)%GEOM_BOX( LOW_IND,JAXIS)
   LZ = GEOMETRY(IG)%GEOM_BOX(HIGH_IND,KAXIS) - GEOMETRY(IG)%GEOM_BOX( LOW_IND,KAXIS)
   MAX_DIST=MAX(MAX_DIST,LX,LY,LZ)
ENDDO

! Set relative epsilon for cut-cell definition:
MAX_DIST= MAX(1._EB,MAX_DIST)
GEOMEPS = GEOMEPS*MAX_DIST

IF(MY_RANK==0 .AND. GET_CUTCELLS_VERBOSE) THEN
   WRITE(LU_ERR,*) 'GEOMETRY intersection computation THRESHOLD GEOMEPS=',GEOMEPS
ENDIF

! Set CCVOL_LINK an epsilon higher than defined value to have all cells/faces around defined value linked.
CCVOL_LINK = CCVOL_LINK + GEOMEPS

IF (PERIODIC_TEST == 105) THEN ! Set cc-guard to zero, i.e. do not compute guard-cell cut-cells, for timings.
   NGUARD = 2
   CCGUARD= NGUARD-2
ENDIF

TNOW = CURRENT_TIME()
CALL SET_CUTCELLS_3D                    ! Defines CUT_CELL data for each mesh.
IF (STOP_STATUS==SETUP_STOP) RETURN

TDEL = CURRENT_TIME() - TNOW

IF (PERIODIC_TEST == 105) THEN ! Cut-cell definition timings test.
    IF(MY_RANK==0) WRITE(LU_ERR,*) ' '
    ICALL = 1
    IF(MY_RANK==0) WRITE(LU_ERR,*) 'CALL number ',ICALL,' to SET_CUTCELLS_3D finished. Max Time=',TDEL,' sec.'
    DO ICALL=2,N_SET_CUTCELLS_3D_CALLS
       TNOW = CURRENT_TIME()
       CALL SET_CUTCELLS_3D                    ! Defines CUT_CELL data for each mesh, average timings.
       TDEL = CURRENT_TIME() - TNOW
       IF(MY_RANK==0) WRITE(LU_ERR,*) 'CALL number ',ICALL,' to SET_CUTCELLS_3D finished. Max Time=',TDEL,' sec.'
    ENDDO
    WRITE_SET_CUTCELLS_TIMINGS = .TRUE.
ENDIF

! Write out SET_CUTCELLS_3D loop time:
IF (WRITE_SET_CUTCELLS_TIMINGS) THEN

   ! Total number of cut-cells and faces computed does not consider guard-cells:
   N_CUTCELLS_PROC     = 0
   N_INB_CUTFACES_PROC = 0
   N_REG_CUTFACES_PROC = 0
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      ! Cut-cells:
      N_CUTCELLS_PROC = N_CUTCELLS_PROC + MESHES(NM)%N_CUTCELL_MESH
      ! Cut-faces:
      DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
         SELECT CASE(CUT_FACE(ICF)%STATUS)
         CASE(CC_GASPHASE)
            N_REG_CUTFACES_PROC = N_REG_CUTFACES_PROC + CUT_FACE(ICF)%NFACE
         CASE(CC_INBOUNDARY)
            N_INB_CUTFACES_PROC = N_INB_CUTFACES_PROC + CUT_FACE(ICF)%NFACE
         END SELECT
      ENDDO
   ENDDO

   ! Write xxx_cc_cpu_0001.csv
   ! This csv file contains the following fields (14):
   ! N_CUTCELLS, N_INB_CUTFACES, N_REG_CUTFACES, SET_CUTCELLS_TIME, GET_BODINT_PLANE_TIME, GET_X2_INTERSECTIONS_TIME, &
   ! GET_X2_VERTVAR_TIME, GET_CARTEDGE_CUTEDGES_TIME, GET_BODX2X3_INTERSECTIONS_TIME, GET_CARTFACE_CUTEDGES_TIME, &
   ! GET_CARTCELL_CUTEDGES_TIME, GET_CARTFACE_CUTFACES_TIME, GET_CARTCELL_CUTFACES_TIME, GET_CARTCELL_CUTCELLS_TIME
   WRITE(FN_CCTIME,'(A,A,I3.3,A)') TRIM(CHID),'_cc_cpu_',MY_RANK,'.csv'
   OPEN(333,FILE=TRIM(FN_CCTIME),STATUS='UNKNOWN')
   WRITE(333,'(A,A,A,A)') "N_CUTCELLS, N_INB_CUTFACES, N_REG_CUTFACES, SET_CUTCELLS_TIME, GET_BODINT_PLANE_TIME, ",   &
                          "GET_X2_INTERSECTIONS_TIME, GET_X2_VERTVAR_TIME, GET_CARTEDGE_CUTEDGES_TIME, ",             &
                          "GET_BODX2X3_INTERSECTIONS_TIME, GET_CARTFACE_CUTEDGES_TIME, GET_CARTCELL_CUTEDGES_TIME, ", &
                          "GET_CARTFACE_CUTFACES_TIME, GET_CARTCELL_CUTFACES_TIME, GET_CARTCELL_CUTCELLS_TIME"
   WRITE(TCFORM,'(23A)')  "(I6,',',I6,',',I6,',',",FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",           &
                          FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",FMT_R,",',',",FMT_R,")"
   WRITE(333,TCFORM) N_CUTCELLS_PROC,N_INB_CUTFACES_PROC,N_REG_CUTFACES_PROC, &
                     T_CC_USED(SET_CUTCELLS_TIME_INDEX:GET_CARTCELL_CUTCELLS_TIME_INDEX)/ &
                     REAL(N_SET_CUTCELLS_3D_CALLS,EB)
   CLOSE(333)


   IF (MY_RANK == 0) THEN
      WRITE(LU_ERR,*) ' '
      WRITE(LU_ERR,*) 'Spheres NVERTS,NFACES',GEOMETRY(1)%N_VERTS,GEOMETRY(1)%N_FACES
      WRITE(LU_ERR,*) 'SET_CUTCELLS_3D loop time by process ',MY_RANK,' =',T_CC_USED(SET_CUTCELLS_TIME_INDEX), &
                      ' sec., cut-cells=',N_CUTCELLS_PROC,', cut-faces=',N_INB_CUTFACES_PROC,N_REG_CUTFACES_PROC
   ENDIF
ENDIF

! Redefine interpolated external wall_cells inside Geoms: We assume them SOLID_BOUNDARY
CALL BLOCK_CC_SOLID_EXTWALLCELLS(FIRST_CALL)

ELSE SET_CUTCELLS_CALL_IF

IF (GET_CUTCELLS_VERBOSE .AND. MY_RANK==0) THEN
   CALL CPU_TIME(TNOW)
ENDIF

! Redefine wall_cells inside Geoms: This is done before EDGE info as edges with WALL_CELL type NULL_BOUNDARY will be taken
! care of by GEOM edges. Note EDGE_INDEX will be reassigned the IBEDGE position in OMEGA, TAU arrays for velocity flux to
! be computed correctly.

CALL BLOCK_CC_SOLID_EXTWALLCELLS(FIRST_CALL)

! Reallocate and populate FDS edge and cell topology variables

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (.NOT.CC_ONLY_IBEDGES_FLAG) CALL GET_REGULAR_CUT_EDGES_BC(NM)
   CALL GET_SOLID_CUTCELL_EDGES_BC(NM)
ENDDO

CALL EXCHANGE_GEOMETRY_INFO

CALL GET_CRTCFCC_INT_STENCILS ! Computes interpolation stencils for face and cell centers.
IF (GET_CUTCELLS_VERBOSE .AND. MY_RANK==0) THEN
   CALL CPU_TIME(TDEL)
   WRITE(LU_ERR,'(A,F8.3,A)') ' Executed GET_CRTCFCC_INT_STENCILS. Time taken : ',TDEL-TNOW,' sec.'
ENDIF
CALL SET_CC_MATVEC_DATA              ! Defines data for discretization matrix-vectors.
IF (GET_CUTCELLS_VERBOSE .AND. MY_RANK==0) THEN
   CALL CPU_TIME(TNOW)
   WRITE(LU_ERR,'(A,F8.3,A)') ' Executing SET_CC_MATVEC_DATA. Time taken : ',TNOW-TDEL,' sec.'
ENDIF
CALL SET_CFACES_P1_RDN               ! Set inverse DXN for CFACES, uses cell linking information.
IF (GET_CUTCELLS_VERBOSE .AND. MY_RANK==0) THEN
   CALL CPU_TIME(TDEL)
   WRITE(LU_ERR,'(A,F8.3,A)') ' Executing SET_CFACES_P1_RDN. Time taken : ',TDEL-TNOW,' sec.'
ENDIF

! Give information for a particular cell:
IF (DEBUG_MATVEC_DATA) THEN
   NM = 1; I = 17; J = 24; K = 19
   CALL POINT_TO_MESH(NM)
   IF (CCVAR(I,J,K,CC_IDCC)<1) THEN
      WRITE(LU_ERR,*) 'No cut-cells in cell NM,I,J,K=',NM,I,J,K
   ELSE
      ICC = CCVAR(I,J,K,CC_IDCC); CC=>CUT_CELL(ICC)
      WRITE(LU_ERR,*) 'Cut Cell in cell NM,I,J,K : ICC,NCELL=',NM,I,J,K,':',ICC,CC%NCELL
      DO JCC=1,CC%NCELL
         ! Linking Info on cut-cell ICC,JCC
         WRITE(LU_ERR,*) 'JCC,ALPHA_CC, CC%IJK_LINK(:,JCC) : CC%LINK_LEV(JCC), UNKZ=',&
         JCC,CC%VOLUME(JCC)/(DX(I)*DY(J)*DZ(K)),CC%IJK_LINK(:,JCC),':',CC%LINK_LEV(JCC),CC%UNKZ(JCC)
         IF (CC%IJK_LINK(1,JCC)==CC_CUTCFE) THEN
            ICC2=CCVAR(CC%IJK_LINK(2,JCC),CC%IJK_LINK(3,JCC),CC%IJK_LINK(4,JCC),CC_IDCC);
            JCC2=CC%IJK_LINK(5,JCC)
            CC2 => CUT_CELL(ICC2)
            ! Linking Info on Parent if it is a cut-face:
            WRITE(LU_ERR,*) 'Parent CC2 I,J,K,JCC2 : ALPHA_CC2, CC2%IJK_LINK(:,JCC2) : CC2%LINK_LEV(JCC2) UNKZ=',&
            CC%IJK_LINK(2:5,JCC),':',CC2%VOLUME(JCC2)/(DX(I)*DY(J)*DZ(K)),CC2%IJK_LINK(:,JCC2),':',&
            CC2%LINK_LEV(JCC2),CC2%UNKZ(JCC2)
         ENDIF
         ! Linking info on cut-faces for ICC,JCC:
         DO JCF=2,CC%CCELEM(1,JCC)+1
            IFACE = CC%CCELEM(JCF,JCC)
            FTYPE = CC%FACE_LIST(1,IFACE)
            IF (FTYPE==CC_FTYPE_CFGAS) THEN
               IFC2    = CC%FACE_LIST(4,IFACE)
               IFACE2  = CC%FACE_LIST(5,IFACE)
               CF => CUT_FACE(IFC2)
               IF(CF%IJK(KAXIS+1)==IAXIS) ACRT = DY(J)*DZ(K)
               IF(CF%IJK(KAXIS+1)==JAXIS) ACRT = DX(I)*DZ(K)
               IF(CF%IJK(KAXIS+1)==KAXIS) ACRT = DX(I)*DY(J)
               WRITE(LU_ERR,*) 'CC Cut-face JFC,ICF,JCF,ALPHA_CF : CF%LINK_LEV(JCF),CF%UNKF(JCF)=',&
               JCF,IFC2,IFACE2,CF%AREA(IFACE2)/ACRT,':',CF%LINK_LEV(IFACE2),CF%UNKF(IFACE2)
               WRITE(LU_ERR,*) 'CUT-FACES with UNKF(JCF) : ',CF%UNKF(IFACE2)
               DO ICF2=1,MESHES(NM)%N_CUTFACE_MESH
                  CF2=>CUT_FACE(ICF2); IF(CF2%STATUS/=CC_GASPHASE) CYCLE
                  DO JCF2=1,CF2%NFACE
                     IF(CF2%UNKF(JCF2)==CF%UNKF(IFACE2)) THEN
                        WRITE(LU_ERR,*) 'CF with UNKF(JCF) ICF2,JCF2,I2,J2,K2,AX2,ALPHA_2,LINK_LEV2=',&
                        ICF2,JCF2,CF2%IJK(1:4),CF2%AREA(JCF2)/ACRT,':',CF2%LINK_LEV(JCF2),CF2%UNKF(JCF2)
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDIF
ENDIF



IF(GET_CUTCELLS_VERBOSE) CLOSE(LU_SETCC)

! Set flag that specifies cut-cell data as defined:
CC_MATVEC_DEFINED=.TRUE.

ENDIF SET_CUTCELLS_CALL_IF

T_USED(14) = T_USED(14) + CURRENT_TIME() - TNOW2

IF (TIME_CC_IBM) T_CC_USED(CC_SET_DATA_TIME_INDEX) = T_CC_USED(CC_SET_DATA_TIME_INDEX) + CURRENT_TIME() - TNOW2
RETURN

CONTAINS

! ------------------------ SET_CFACES_P1_RDN ---------------------------------

SUBROUTINE SET_CFACES_P1_RDN

! Local Variables:
INTEGER :: ICF, IFACE
INTEGER :: ICC, JCC, I, J, K
INTEGER :: IFACE_CELL, ICF_CELL, IROW
INTEGER :: IROW_DELTA
REAL(EB):: AREAI
REAL(EB), ALLOCATABLE, DIMENSION(:) :: DXN_UNKZ_LOC, AREA_UNKZ_LOC, VOL_UNKZ_LOC
INTEGER, PARAMETER :: RDN_METHOD = 2 ! 1: VOL/SOLIDAREA, 2: bbox projection
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: XYZMIN_UNKZ, XYZMAX_UNKZ
REAL(EB) :: DELTA_X, DELTA_Y, DELTA_Z, DXN, NVEC_ABS(MAX_DIM)
TYPE(CFACE_TYPE), POINTER :: CFA
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

! ALLOCATE local arrays
ALLOCATE(DXN_UNKZ_LOC(1:NUNKZ_LOCAL));  DXN_UNKZ_LOC(:)  = 0._EB
ALLOCATE(AREA_UNKZ_LOC(1:NUNKZ_LOCAL)); AREA_UNKZ_LOC(:) = 0._EB
ALLOCATE(VOL_UNKZ_LOC(1:NUNKZ_LOCAL));  VOL_UNKZ_LOC(:)  = 0._EB
IF (RDN_METHOD==1) THEN
   ! Main Loop:
   MESH_LOOP_01 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

      CALL POINT_TO_MESH(NM)

      ! Do a volume weighted average of distance to wall from linked cells, if one of them is a regular cell use 1/2 the
      ! distance of corner to corner sqrt(DX^2+DY^2+DZ^2).
      ! 1. Regular GASPHASE cells within the cc-region:
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,CC_UNKZ) <= 0 ) CYCLE ! Drop if regular gas cell has not been assigned unknown number.
               IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START) ! All row indexes must refer to ind_loc.
               VOL_UNKZ_LOC(IROW) = VOL_UNKZ_LOC(IROW) + (DX(I)*DY(J)*DZ(K))
            ENDDO
         ENDDO
      ENDDO
      ! 2. Cut-cells:
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
         IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
         DO JCC=1,CC%NCELL
            IROW = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
            ! Mean INBOUNDARY cut-face distance to this cut-cell center, projected to cut-face normal:
            AREAI = 0._EB
            DO ICF_CELL=1,CC%CCELEM(1,JCC)
               IFACE_CELL = CC%CCELEM(ICF_CELL+1,JCC)
               IF (CC%FACE_LIST(1,IFACE_CELL) /= CC_FTYPE_CFINB) CYCLE
               ! Indexes of INBOUNDARY cutface on CUT_FACE:
               ICF   = CC%FACE_LIST(4,IFACE_CELL)
               IFACE = CC%FACE_LIST(5,IFACE_CELL)
               ! Area sum:
               AREAI = AREAI + CUT_FACE(ICF)%AREA(IFACE)
            ENDDO
            AREA_UNKZ_LOC(IROW) = AREA_UNKZ_LOC(IROW) + AREAI
            VOL_UNKZ_LOC(IROW) = VOL_UNKZ_LOC(IROW) + CC%VOLUME(JCC)
         ENDDO
      ENDDO

   ENDDO MESH_LOOP_01

   ! Compute average DXN of all linked cells:
   DXN_UNKZ_LOC = VOL_UNKZ_LOC / (AREA_UNKZ_LOC + TWO_EPSILON_EB)
ELSE
   ALLOCATE(XYZMIN_UNKZ(IAXIS:KAXIS,1:NUNKZ_LOCAL))
   ALLOCATE(XYZMAX_UNKZ(IAXIS:KAXIS,1:NUNKZ_LOCAL))
   CALL GET_LINKED_CELLS_BBOX(XYZMIN_UNKZ,XYZMAX_UNKZ)
   IF (ALLOCATED(DELTA_UNKZ)) THEN
      IF (SIZE(DELTA_UNKZ,DIM=2) /= NUNKZ_LOCAL) DEALLOCATE(DELTA_UNKZ)
   ENDIF
   IF (.NOT.ALLOCATED(DELTA_UNKZ)) ALLOCATE(DELTA_UNKZ(IAXIS:KAXIS,1:NUNKZ_LOCAL))
   DO IROW_DELTA=1,NUNKZ_LOCAL
      DELTA_UNKZ(IAXIS,IROW_DELTA) = MAX(0._EB,XYZMAX_UNKZ(IAXIS,IROW_DELTA)-XYZMIN_UNKZ(IAXIS,IROW_DELTA))
      DELTA_UNKZ(JAXIS,IROW_DELTA) = MAX(0._EB,XYZMAX_UNKZ(JAXIS,IROW_DELTA)-XYZMIN_UNKZ(JAXIS,IROW_DELTA))
      DELTA_UNKZ(KAXIS,IROW_DELTA) = MAX(0._EB,XYZMAX_UNKZ(KAXIS,IROW_DELTA)-XYZMIN_UNKZ(KAXIS,IROW_DELTA))
   ENDDO
ENDIF

! Finally Define B1%RDN:
MESH_LOOP_02 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
      CF => CUT_FACE(ICF); IF(CF%STATUS /= CC_INBOUNDARY) CYCLE
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO IFACE=1,CF%NFACE
         IF (CF%CELL_LIST(1,LOW_IND,IFACE) /= CC_FTYPE_CFGAS) CYCLE
         ICC = CF%CELL_LIST(2,LOW_IND,IFACE)
         JCC = CF%CELL_LIST(3,LOW_IND,IFACE)
         IROW = CUT_CELL(ICC)%UNKZ(JCC) - UNKZ_IND(NM_START)
         SELECT CASE (RDN_METHOD)
         CASE (1)
            BOUNDARY_PROP1(CFACE(CF%CFACE_INDEX(IFACE))%B1_INDEX)%RDN = 1._EB/DXN_UNKZ_LOC(IROW)
         CASE (2)
            CFA => CFACE(CF%CFACE_INDEX(IFACE))
            BC  => BOUNDARY_COORD(CFA%BC_INDEX)
            DELTA_X = XYZMAX_UNKZ(IAXIS,IROW) - XYZMIN_UNKZ(IAXIS,IROW)
            DELTA_Y = XYZMAX_UNKZ(JAXIS,IROW) - XYZMIN_UNKZ(JAXIS,IROW)
            DELTA_Z = XYZMAX_UNKZ(KAXIS,IROW) - XYZMIN_UNKZ(KAXIS,IROW)
            NVEC_ABS(IAXIS:KAXIS) = ABS(BC%NVEC(IAXIS:KAXIS))
            DXN = DELTA_X*NVEC_ABS(IAXIS) + DELTA_Y*NVEC_ABS(JAXIS) + DELTA_Z*NVEC_ABS(KAXIS)
            BOUNDARY_PROP1(CFA%B1_INDEX)%RDN = 1._EB/(DXN + TWO_EPSILON_EB)
         END SELECT
      ENDDO
   ENDDO
ENDDO MESH_LOOP_02
IF (ALLOCATED(XYZMIN_UNKZ)) DEALLOCATE(XYZMIN_UNKZ, XYZMAX_UNKZ)
DEALLOCATE(DXN_UNKZ_LOC, VOL_UNKZ_LOC, AREA_UNKZ_LOC)

RETURN
END SUBROUTINE SET_CFACES_P1_RDN

! ---------------------- GET_LINKED_CELLS_BBOX --------------------------------

SUBROUTINE GET_LINKED_CELLS_BBOX(XYZMIN_UNKZ,XYZMAX_UNKZ)

! Local Variables:
INTEGER :: I, J, K, ICC, JCC, IROW
INTEGER :: IFACE, IFC, FTYPE, X1AXIS, LOWHIGH, ILH, IFC2, IFACE2, NVFACE, IPT, IVERT
REAL(EB) :: XLO, XHI, YLO, YHI, ZLO, ZHI
REAL(EB), INTENT(OUT) :: XYZMIN_UNKZ(IAXIS:KAXIS,1:NUNKZ_LOCAL)
REAL(EB), INTENT(OUT) :: XYZMAX_UNKZ(IAXIS:KAXIS,1:NUNKZ_LOCAL)

XYZMIN_UNKZ =  HUGE(1._EB); XYZMAX_UNKZ = -HUGE(1._EB)
MESH_LOOP_01 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0 ) CYCLE
            IROW = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            XLO = X(I-1); XHI = X(I)
            YLO = Y(J-1); YHI = Y(J)
            ZLO = Z(K-1); ZHI = Z(K)
            XYZMIN_UNKZ(IAXIS,IROW) = MIN(XYZMIN_UNKZ(IAXIS,IROW),XLO)
            XYZMIN_UNKZ(JAXIS,IROW) = MIN(XYZMIN_UNKZ(JAXIS,IROW),YLO)
            XYZMIN_UNKZ(KAXIS,IROW) = MIN(XYZMIN_UNKZ(KAXIS,IROW),ZLO)
            XYZMAX_UNKZ(IAXIS,IROW) = MAX(XYZMAX_UNKZ(IAXIS,IROW),XHI)
            XYZMAX_UNKZ(JAXIS,IROW) = MAX(XYZMAX_UNKZ(JAXIS,IROW),YHI)
            XYZMAX_UNKZ(KAXIS,IROW) = MAX(XYZMAX_UNKZ(KAXIS,IROW),ZHI)
         ENDDO
      ENDDO
   ENDDO

   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         XLO =  HUGE(1._EB); XHI = -HUGE(1._EB)
         YLO =  HUGE(1._EB); YHI = -HUGE(1._EB)
         ZLO =  HUGE(1._EB); ZHI = -HUGE(1._EB)
         DO IFC=1,CC%CCELEM(1,JCC)
            IFACE = CC%CCELEM(IFC+1,JCC)
            FTYPE = CC%FACE_LIST(1,IFACE)
            SELECT CASE (FTYPE)
            CASE (CC_FTYPE_CFGAS,CC_FTYPE_CFINB)
               IFC2   = CC%FACE_LIST(4,IFACE)
               IFACE2 = CC%FACE_LIST(5,IFACE)
               NVFACE = CUT_FACE(IFC2)%CFELEM(1,IFACE2)
               DO IPT=1,NVFACE
                  IVERT = CUT_FACE(IFC2)%CFELEM(IPT+1,IFACE2)
                  XLO = MIN(XLO,CUT_FACE(IFC2)%XYZVERT(IAXIS,IVERT))
                  XHI = MAX(XHI,CUT_FACE(IFC2)%XYZVERT(IAXIS,IVERT))
                  YLO = MIN(YLO,CUT_FACE(IFC2)%XYZVERT(JAXIS,IVERT))
                  YHI = MAX(YHI,CUT_FACE(IFC2)%XYZVERT(JAXIS,IVERT))
                  ZLO = MIN(ZLO,CUT_FACE(IFC2)%XYZVERT(KAXIS,IVERT))
                  ZHI = MAX(ZHI,CUT_FACE(IFC2)%XYZVERT(KAXIS,IVERT))
               ENDDO
            CASE (CC_FTYPE_RCGAS)
               LOWHIGH = CC%FACE_LIST(2,IFACE)
               X1AXIS  = CC%FACE_LIST(3,IFACE)
               ILH     = LOWHIGH - 1
               SELECT CASE (X1AXIS)
               CASE (IAXIS)
                  XLO = MIN(XLO,X(I-1+ILH)); XHI = MAX(XHI,X(I-1+ILH))
                  YLO = MIN(YLO,Y(J-1));     YHI = MAX(YHI,Y(J))
                  ZLO = MIN(ZLO,Z(K-1));     ZHI = MAX(ZHI,Z(K))
               CASE (JAXIS)
                  XLO = MIN(XLO,X(I-1));     XHI = MAX(XHI,X(I))
                  YLO = MIN(YLO,Y(J-1+ILH)); YHI = MAX(YHI,Y(J-1+ILH))
                  ZLO = MIN(ZLO,Z(K-1));     ZHI = MAX(ZHI,Z(K))
               CASE (KAXIS)
                  XLO = MIN(XLO,X(I-1));     XHI = MAX(XHI,X(I))
                  YLO = MIN(YLO,Y(J-1));     YHI = MAX(YHI,Y(J))
                  ZLO = MIN(ZLO,Z(K-1+ILH)); ZHI = MAX(ZHI,Z(K-1+ILH))
               END SELECT
            END SELECT
         ENDDO
         IF (XLO > XHI) THEN
            XLO = X(I-1); XHI = X(I)
            YLO = Y(J-1); YHI = Y(J)
            ZLO = Z(K-1); ZHI = Z(K)
         ENDIF
         XYZMIN_UNKZ(IAXIS,IROW) = MIN(XYZMIN_UNKZ(IAXIS,IROW),XLO)
         XYZMIN_UNKZ(JAXIS,IROW) = MIN(XYZMIN_UNKZ(JAXIS,IROW),YLO)
         XYZMIN_UNKZ(KAXIS,IROW) = MIN(XYZMIN_UNKZ(KAXIS,IROW),ZLO)
         XYZMAX_UNKZ(IAXIS,IROW) = MAX(XYZMAX_UNKZ(IAXIS,IROW),XHI)
         XYZMAX_UNKZ(JAXIS,IROW) = MAX(XYZMAX_UNKZ(JAXIS,IROW),YHI)
         XYZMAX_UNKZ(KAXIS,IROW) = MAX(XYZMAX_UNKZ(KAXIS,IROW),ZHI)
      ENDDO
   ENDDO
ENDDO MESH_LOOP_01

RETURN
END SUBROUTINE GET_LINKED_CELLS_BBOX

END SUBROUTINE CC_SET_DATA


! ----------------------------- INIT_CUTCELL_DATA -------------------------------

SUBROUTINE INIT_CUTCELL_DATA(T,DT,FIRST_CALL)

! This routine assumes INITIALIZE_MESH_VARIABLES_1(DT,NM) has already been called.

USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

REAL(EB), INTENT(IN) :: T, DT
LOGICAL, INTENT(IN)  :: FIRST_CALL

! Local Variables:
INTEGER :: NM,I,J,K,N,ICC,JCC,X1AXIS,NFACE,ICF
REAL(EB) TMP_CC,RHO_CC,AREAT,VEL_CF !,Z_CC,TMP_0_CC,P_0_CC
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_CC
INTEGER :: IW,IROW_LOC

REAL(EB) :: TNOW

TNOW = CURRENT_TIME()

ALLOCATE( ZZ_CC(1:N_TOTAL_SCALARS) )

! Loop Meshes:
! Get Z location of linked cells centroids, then P_0_CV for the control volumes:
ZCEN_CV(:) = 0._EB; RZ_Z = 0._EB
MESH_LOOP_0 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
            IROW_LOC     = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
            RZ_Z(IROW_LOC)= RZ_Z(IROW_LOC) + DX(I)*DY(J)*DZ(K)
            ZCEN_CV(IROW_LOC) = ZCEN_CV(IROW_LOC) + ZC(K)*DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IF (CELL(CELL_INDEX(CC%IJK(IAXIS),CC%IJK(JAXIS),CC%IJK(KAXIS)))%SOLID) CYCLE
      DO JCC=1,CC%NCELL
         IROW_LOC     = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         RZ_Z(IROW_LOC)= RZ_Z(IROW_LOC) + CC%VOLUME(JCC)
         ZCEN_CV(IROW_LOC) = ZCEN_CV(IROW_LOC) + CC%XYZCEN(KAXIS,JCC)*CC%VOLUME(JCC)
      ENDDO
   ENDDO
ENDDO MESH_LOOP_0
DO IROW_LOC=1,NUNKZ_LOCAL
   ZCEN_CV(IROW_LOC) = ZCEN_CV(IROW_LOC) / RZ_Z(IROW_LOC)
ENDDO

P_0_CV(:)  = P_INF
TMP_0_CV(:)= TMPA
IF (STRATIFICATION) THEN
   DO IROW_LOC=1,NUNKZ_LOCAL
      P_0_CV(IROW_LOC)  = EVALUATE_RAMP(ZCEN_CV(IROW_LOC),I_RAMP_P0_Z)
      TMP_0_CV(IROW_LOC)= TMPA*EVALUATE_RAMP(ZCEN_CV(IROW_LOC),I_RAMP_TMP0_Z)
   ENDDO
ENDIF
DO IROW_LOC=1,NUNKZ_LOCAL
   RHO_0_CV(IROW_LOC) = P_0_CV(IROW_LOC)/(TMP_0_CV(IROW_LOC)*RSUM0)
ENDDO

MESH_LOOP_1 : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! Default initialization:
   IF(.NOT.RESTART .AND. PERIODIC_TEST/=11 .AND. PERIODIC_TEST/=7) THEN
      DO K=1,KBAR ! Linked cells get initialized to control volume values of TMP_0 and RHO_0
         DO J=1,JBAR
            DO I=1,IBAR
               IF(CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE
               IROW_LOC     = CCVAR(I,J,K,CC_UNKZ) - UNKZ_IND(NM_START)
               TMP(I,J,K) = TMP_0_CV(IROW_LOC)
               RHO(I,J,K) = RHO_0_CV(IROW_LOC)
               RHOS(I,J,K)= RHO_0_CV(IROW_LOC)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Cut-cells inherit underlying Cartesian cell values of rho,T,Z, etc.:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC);  I=CC%IJK(IAXIS); J=CC%IJK(JAXIS);  K=CC%IJK(KAXIS); IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      TMP_CC = TMP(I,J,K)
      RHO_CC = RHO(I,J,K)
      ZZ_CC(1:N_TOTAL_SCALARS) = ZZ(I,J,K,1:N_TOTAL_SCALARS)
      DO JCC=1,CC%NCELL
         IROW_LOC     = CC%UNKZ(JCC) - UNKZ_IND(NM_START)
         CC%RHO_0(JCC)= RHO_0_CV(IROW_LOC)
         IF (.NOT.RESTART .AND. PERIODIC_TEST/=11 .AND. PERIODIC_TEST/=7) THEN
           TMP_CC = TMP_0_CV(IROW_LOC)
           RHO_CC = RHO_0_CV(IROW_LOC)
         ENDIF

         CC%TMP(JCC)  = TMP_CC
         CC%RHO(JCC)  = RHO_CC
         CC%RHOS(JCC) = RHO_CC

         CC%ZZ(1:N_TOTAL_SCALARS,JCC) = ZZ_CC(1:N_TOTAL_SCALARS)
         DO N=1,N_TRACKED_SPECIES
            CC%ZZS(N,JCC) = SPECIES_MIXTURE(N)%ZZ0
         ENDDO
         CC%MIX_TIME(JCC) = DT
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_CC(1:N_TRACKED_SPECIES),CC%RSUM(JCC))
         CC%D(JCC)        = 0._EB
         CC%DS(JCC)       = 0._EB
         CC%DVOL(JCC)     = 0._EB
         CC%D_SOURCE(JCC) = 0._EB
         CC%Q(JCC)        = 0._EB
         CC%QR(JCC)       = 0._EB
         CC%M_DOT_PPP(:,JCC) = 0._EB
         IF(RESTART) THEN
            CC%D(JCC)        = D(I,J,K)/CC%ALPHA_CC
            IF (ALLOCATED(MESHES(NM)%D_SOURCE)) CC%D_SOURCE(JCC) = D_SOURCE(I,J,K)/CC%ALPHA_CC
            CC%Q(JCC)        = Q(I,J,K)/CC%ALPHA_CC
            CC%QR(JCC)       = QR(I,J,K) ! Not needed for radiation.
            IF (ALLOCATED(MESHES(NM)%M_DOT_PPP)) CC%M_DOT_PPP(1:N_TRACKED_SPECIES,JCC) = &
                                                    M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES)/CC%ALPHA_CC
         ENDIF
      ENDDO
      IF (ALLOCATED(MESHES(NM)%D_SOURCE)) THEN
         D_SOURCE(I,J,K) = 0._EB
         M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES) = 0._EB
      ENDIF
   ENDDO

   ! Init guardcell cut-cells:
   DO ICC=MESHES(NM)%N_CUTCELL_MESH+1,MESHES(NM)%N_CUTCELL_MESH+MESHES(NM)%N_GCCUTCELL_MESH
      CC => CUT_CELL(ICC);  I = CC%IJK(IAXIS); J = CC%IJK(JAXIS);  K = CC%IJK(KAXIS)
      IF(I < 0 .OR. I > IBP1) CYCLE
      IF(J < 0 .OR. J > JBP1) CYCLE
      IF(K < 0 .OR. K > KBP1) CYCLE
      TMP_CC = TMP(I,J,K)
      RHO_CC = RHO(I,J,K)
      ZZ_CC(1:N_TOTAL_SCALARS) = ZZ(I,J,K,1:N_TOTAL_SCALARS)
      DO JCC=1,CC%NCELL
         CC%RHO_0(JCC)= RHO_0(K)
         CC%TMP(JCC)  = TMP_CC
         CC%RHO(JCC)  = RHO_CC
         CC%RHOS(JCC) = RHO_CC
         CC%ZZ(1:N_TOTAL_SCALARS,JCC) = ZZ_CC(1:N_TOTAL_SCALARS)
         DO N=1,N_TRACKED_SPECIES
            CC%ZZS(N,JCC) = SPECIES_MIXTURE(N)%ZZ0
         ENDDO
         CC%MIX_TIME(JCC) = DT
         CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_CC(1:N_TRACKED_SPECIES),CC%RSUM(JCC))
         CC%D(JCC)        = 0._EB
         CC%DS(JCC)       = 0._EB
         CC%DVOL(JCC)     = 0._EB
         CC%D_SOURCE(JCC) = 0._EB
         CC%Q(JCC)        = 0._EB
         CC%QR(JCC)       = 0._EB
         CC%M_DOT_PPP(:,JCC) = 0._EB
         IF(RESTART) THEN
            CC%D(JCC)        = D(I,J,K)/CC%ALPHA_CC
            IF (ALLOCATED(MESHES(NM)%D_SOURCE)) CC%D_SOURCE(JCC) = D_SOURCE(I,J,K)/CC%ALPHA_CC
            CC%Q(JCC)        = Q(I,J,K)/CC%ALPHA_CC
            CC%QR(JCC)       = QR(I,J,K) ! Not needed for radiation.
            IF (ALLOCATED(MESHES(NM)%M_DOT_PPP)) CC%M_DOT_PPP(1:N_TRACKED_SPECIES,JCC) = &
                                                    M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES)/CC%ALPHA_CC
         ENDIF
      ENDDO
      IF (ALLOCATED(MESHES(NM)%D_SOURCE)) THEN
         D_SOURCE(I,J,K) = 0._EB
         M_DOT_PPP(I,J,K,1:N_TRACKED_SPECIES) = 0._EB
      ENDIF
   ENDDO

   ! Gasphase Cut-faces inherit underlying Cartesian face values of Velocity (flux matched or not):
   PERIODIC_TEST_COND : IF (PERIODIC_TEST /= 21 .AND. PERIODIC_TEST /= 22 .AND. PERIODIC_TEST /= 24) THEN

      ! First GASPHASe cut-faces:
      CUTFACE_LOOP : DO ICF=1,MESHES(NM)%N_CUTFACE_MESH+MESHES(NM)%N_GCCUTFACE_MESH
         NFACE  = CUT_FACE(ICF)%NFACE
         IF (CUT_FACE(ICF)%STATUS /= CC_GASPHASE) CYCLE
         I      = CUT_FACE(ICF)%IJK(IAXIS); IF(I<0 .OR. I>IBP1) CYCLE
         J      = CUT_FACE(ICF)%IJK(JAXIS); IF(J<0 .OR. J>JBP1) CYCLE
         K      = CUT_FACE(ICF)%IJK(KAXIS); IF(K<0 .OR. K>KBP1) CYCLE
         X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)

         AREAT  = SUM( CUT_FACE(ICF)%AREA(1:NFACE) )

         ! Flux matched U to cut-face centroids, they all get same velocity:
         IF(RESTART) THEN
            SELECT CASE(X1AXIS)
            CASE(IAXIS); VEL_CF = (DY(J)*DZ(K))/(AREAT+TWENTY_EPSILON_EB) * U(I,J,K)
            CASE(JAXIS); VEL_CF = (DX(I)*DZ(K))/(AREAT+TWENTY_EPSILON_EB) * V(I,J,K)
            CASE(KAXIS); VEL_CF = (DX(I)*DY(J))/(AREAT+TWENTY_EPSILON_EB) * W(I,J,K)
            END SELECT
         ELSE
            SELECT CASE(X1AXIS)
            CASE(IAXIS); VEL_CF = U(I,J,K)
            CASE(JAXIS); VEL_CF = V(I,J,K)
            CASE(KAXIS); VEL_CF = W(I,J,K)
            END SELECT
         ENDIF

         CUT_FACE(ICF)%VEL(1:NFACE)  = VEL_CF
         CUT_FACE(ICF)%VELS(1:NFACE) = VEL_CF
         CUT_FACE(ICF)%VEL_CF        = VEL_CF
      ENDDO CUTFACE_LOOP

      ! Push cut-face velocities down to cartesian velocities:
      DO ICF=1,MESHES(NM)%N_CUTFACE_MESH+MESHES(NM)%N_GCCUTFACE_MESH
         NFACE  = CUT_FACE(ICF)%NFACE
         IF (CUT_FACE(ICF)%STATUS /= CC_GASPHASE) CYCLE
         I      = CUT_FACE(ICF)%IJK(IAXIS); IF(I<0 .OR. I>IBP1) CYCLE
         J      = CUT_FACE(ICF)%IJK(JAXIS); IF(J<0 .OR. J>JBP1) CYCLE
         K      = CUT_FACE(ICF)%IJK(KAXIS); IF(K<0 .OR. K>KBP1) CYCLE
         X1AXIS = CUT_FACE(ICF)%IJK(KAXIS+1)
         VEL_CF = DOT_PRODUCT(CUT_FACE(ICF)%VEL(1:NFACE),CUT_FACE(ICF)%AREA(1:NFACE))
         SELECT CASE(X1AXIS)
         CASE(IAXIS); U(I,J,K) = VEL_CF/(DY(J)*DZ(K)); CUT_FACE(ICF)%VEL_CRT=U(I,J,K)
         CASE(JAXIS); V(I,J,K) = VEL_CF/(DX(I)*DZ(K)); CUT_FACE(ICF)%VEL_CRT=V(I,J,K)
         CASE(KAXIS); W(I,J,K) = VEL_CF/(DX(I)*DY(J)); CUT_FACE(ICF)%VEL_CRT=W(I,J,K)
         END SELECT
      ENDDO

      ! INBOUNDARY cut-faces are initialized with 0._EB velocity, that will be changed in
      ! CFACE_PREDICT_NORMAL_VELOCITY.

   ENDIF PERIODIC_TEST_COND

   ! Force velocities in CC_SOLID faces to zero
   WHERE(FCVAR(0:IBP1,0:JBP1,0:KBP1,CC_FGSC,IAXIS)==CC_SOLID) U(0:IBP1,0:JBP1,0:KBP1) = 0._EB
   WHERE(FCVAR(0:IBP1,0:JBP1,0:KBP1,CC_FGSC,JAXIS)==CC_SOLID) V(0:IBP1,0:JBP1,0:KBP1) = 0._EB
   WHERE(FCVAR(0:IBP1,0:JBP1,0:KBP1,CC_FGSC,KAXIS)==CC_SOLID) W(0:IBP1,0:JBP1,0:KBP1) = 0._EB
   US = U; VS = V; WS = W

   ! External mesh CFACEs initialize P1 BCs:
   DO ICF=1,N_EXTERNAL_CFACE_CELLS+N_INTWALL_CFACE_CELLS
      CFA  => CFACE(ICF)
      IW = CUT_FACE(CFA%CUT_FACE_IND1)%IWC
      CALL INIT_CFACE_CELL(NM,CFA%CUT_FACE_IND1,CFA%CUT_FACE_IND2,ICF,CFA%SURF_INDEX,INTEGER_THREE,&
                           IS_INB=.FALSE.,IW=IW)
   ENDDO

   ! Geometry boundary CFACES initialize P1 BCs:
   IF(.NOT.RESTART) THEN ! Only if not restarting, otherwise the Boundary P1 vars are read from restart file.
      DO ICF=INTERNAL_CFACE_CELLS_LB+1,INTERNAL_CFACE_CELLS_LB+N_INTERNAL_CFACE_CELLS
         CFA  => CFACE(ICF)
         CALL INIT_CFACE_CELL(NM,CFA%CUT_FACE_IND1,CFA%CUT_FACE_IND2,ICF,CFA%SURF_INDEX,INTEGER_THREE,IS_INB=.TRUE.)
      ENDDO
   ENDIF

ENDDO MESH_LOOP_1

DEALLOCATE( ZZ_CC)
IF(((.NOT.RESTART .AND. FIRST_CALL) .OR. (RESTART .AND. .NOT.FIRST_CALL)) .AND. N_GEOMETRY>0) DEALLOCATE(FDS_AREA_GEOM)

! Populate Linked velocity arrays:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   CALL GET_LINKED_VELOCITIES(NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.,CMP_FLG=.TRUE.)
ENDDO

! Flux match Cartesian face velocity back to cut-faces:
CALL CC_VELOCITY_CUTFACES(APPLY_TO_ESTIMATED_VARIABLES=.FALSE.)


CALL MESH_CC_EXCHANGE(1)
CALL MESH_CC_EXCHANGE(4)
CALL MESH_CC_EXCHANGE(6)

CALL CC_H_INTERP
CALL CC_RHO0W_INTERP

! Recompute initial CC_VELOCITY_BC with RESTART data if needed
IF (RESTART .AND. .NOT. FIRST_CALL) THEN
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      DRAG_UVWMAX = 0._EB
      CALL CC_VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=.FALSE.,DO_IBEDGES=.TRUE.)
   ENDDO
ENDIF

! Check divergence of initial velocity field:
IF(GET_CUTCELLS_VERBOSE) CALL CC_CHECK_DIVERGENCE(T,DT,.FALSE.)

T_USED(14) = T_USED(14) + CURRENT_TIME() - TNOW
IF (TIME_CC_IBM) T_CC_USED(INIT_CUTCELL_DATA_TIME_INDEX) = T_CC_USED(INIT_CUTCELL_DATA_TIME_INDEX) + CURRENT_TIME() - TNOW
RETURN
END SUBROUTINE INIT_CUTCELL_DATA



! -------------------------------- FINISH_CC ---------------------------------

SUBROUTINE FINISH_CC

USE MPI_F08

! Local variables:
INTEGER :: I, TLB, TUB, LU_TCC, IERR
CHARACTER(MESSAGE_LENGTH) :: CC_CPU_FILE
REAL(EB), ALLOCATABLE, DIMENSION(:) :: T_CC_USED_MIN, T_CC_USED_MAX, T_CC_USED_MEA
CHARACTER(30) :: FRMT

IF (TIME_CC_IBM) THEN
   TLB = LBOUND(T_CC_USED,DIM=1)
   TUB = UBOUND(T_CC_USED,DIM=1)
   ALLOCATE(T_CC_USED_MIN(TLB:TUB),T_CC_USED_MAX(TLB:TUB),T_CC_USED_MEA(TLB:TUB))
   IF (N_MPI_PROCESSES > 1) THEN
      CALL MPI_ALLREDUCE(T_CC_USED(TLB) , T_CC_USED_MIN(TLB) , TUB-TLB+1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLREDUCE(T_CC_USED(TLB) , T_CC_USED_MAX(TLB) , TUB-TLB+1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLREDUCE(T_CC_USED(TLB) , T_CC_USED_MEA(TLB) , TUB-TLB+1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
      T_CC_USED_MEA = T_CC_USED_MEA / REAL(N_MPI_PROCESSES,EB)
   ELSE
      T_CC_USED_MIN(TLB:TUB) = T_CC_USED(TLB:TUB)
      T_CC_USED_MAX(TLB:TUB) = T_CC_USED(TLB:TUB)
      T_CC_USED_MEA(TLB:TUB) = T_CC_USED(TLB:TUB)
   ENDIF
   IF (MY_RANK==0) THEN
      WRITE(CC_CPU_FILE,'(A,A)') TRIM(CHID),'_cc_cpu.csv'
      LU_TCC = GET_FILE_NUMBER()
      OPEN(LU_TCC,FILE=TRIM(CC_CPU_FILE),STATUS='UNKNOWN')
      WRITE(LU_TCC,'(A,A,A)') 'CCCOMPUTE_RADIATION, CC_DENSITY, CC_VELOCITY_FLUX, CC_COMPUTE_VISCOSITY, ',&
                              'CC_INTERP_FACE_VEL, CC_DIVERGENCE_PART_1, CC_END_STEP, CC_TARGET_VELOCITY, ',&
                              'CC_NO_FLUX, CC_COMPUTE_VELOCITY_ERROR, MESH_CC_EXCHANGE (s)'
      WRITE(FRMT,'(A,I2.2,A)') '(',MESH_CC_EXCHANGE_TIME_INDEX-CCCOMPUTE_RADIATION_TIME_INDEX+1,'(",",ES10.3))'
      WRITE(LU_TCC,FRMT) (T_CC_USED_MIN(I),I=CCCOMPUTE_RADIATION_TIME_INDEX,MESH_CC_EXCHANGE_TIME_INDEX)
      WRITE(LU_TCC,FRMT) (T_CC_USED_MAX(I),I=CCCOMPUTE_RADIATION_TIME_INDEX,MESH_CC_EXCHANGE_TIME_INDEX)
      WRITE(LU_TCC,FRMT) (T_CC_USED_MEA(I),I=CCCOMPUTE_RADIATION_TIME_INDEX,MESH_CC_EXCHANGE_TIME_INDEX)
      CLOSE(LU_TCC)
   ENDIF
   DEALLOCATE(T_CC_USED_MIN,T_CC_USED_MAX,T_CC_USED_MEA)
ENDIF

! Release Requests:
DO I=1,N_REQ11  ; CALL MPI_REQUEST_FREE(REQ11(I) ,IERR) ; ENDDO
DO I=1,N_REQ12  ; CALL MPI_REQUEST_FREE(REQ12(I) ,IERR) ; ENDDO
DO I=1,N_REQ13  ; CALL MPI_REQUEST_FREE(REQ13(I) ,IERR) ; ENDDO

RETURN
END SUBROUTINE FINISH_CC



! --------------------------- GET_CRTCFCC_INT_STENCILS -----------------------------

SUBROUTINE GET_CRTCFCC_INT_STENCILS

USE GEOMETRY_FUNCTIONS, ONLY : SEARCH_OTHER_MESHES

! Local variables:
INTEGER :: NM
INTEGER :: X1AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:)   :: IJKCELL
INTEGER :: I,J,K,NCELL,ICC,JCC,IJK(MAX_DIM),ICF,IFACE
INTEGER :: ISTR, IEND, JSTR, JEND, KSTR, KEND
LOGICAL :: FOUND_POINT, INSEG, FOUNDPT
REAL(EB):: XYZ(MAX_DIM),XYZ_PP(MAX_DIM),XYZ_IP(MAX_DIM),DV(MAX_DIM),NVEC(MAX_DIM)
REAL(EB):: P0(MAX_DIM),P1(MAX_DIM),DIST,DISTANCE,DIR_FCT,NORM_DV,LASTDOTNVEC,DOTNVEC
INTEGER :: FOUND_INBFC(1:3), BODTRI(1:2)
INTEGER :: CCFC,NFC_CC,ICFC,INBFC,INBFC_LOC,IFCPT,IFCPT_LOC,IBOD,IWSEL,ICELL
INTEGER :: TESTVAR

INTEGER :: IW,II,JJ,KK,IGC
REAL(EB) :: MIN_DIST_VEL

! OMESH related arrays:
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: IJKFACE2
INTEGER :: IIO,JJO,KKO,NOM
LOGICAL :: FLGX,FLGY,FLGZ,INNM

INTEGER, ALLOCATABLE, DIMENSION(:) :: IIO_FC_R_AUX,JJO_FC_R_AUX,KKO_FC_R_AUX,AXS_FC_R_AUX
INTEGER, ALLOCATABLE, DIMENSION(:) :: IIO_CC_R_AUX,JJO_CC_R_AUX,KKO_CC_R_AUX
INTEGER :: SIZE_REC

INTEGER, PARAMETER :: DELTA_FC = 200

INTEGER :: VIND,EP,INPE,INT_NPE_LO,INT_NPE_HI,NPE_LIST_START,NPE_LIST_COUNT,SZ_1,SZ_2,NPE_COUNT,IEDGE
INTEGER,  ALLOCATABLE, DIMENSION(:,:,:,:) :: INT_NPE
INTEGER,  ALLOCATABLE, DIMENSION(:,:)     :: INT_IJK, INT_NOMIND
REAL(EB), ALLOCATABLE, DIMENSION(:)       :: INT_COEF
REAL(EB), ALLOCATABLE, DIMENSION(:,:)     :: INT_NOUT, INT_DCOEF
REAL(EB) :: DELN,INT_XN(0:INT_N_EXT_PTS),INT_CN(0:INT_N_EXT_PTS)
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INT_IJK_AUX
REAL(EB),ALLOCATABLE, DIMENSION(:)   :: INT_COEF_AUX
INTEGER :: N_CVAR_START, N_CVAR_COUNT, N_FVAR_START, N_FVAR_COUNT
LOGICAL, ALLOCATABLE, DIMENSION(:) :: EP_TAG

INTEGER :: IS,I_SGN,ICD,ICD_SGN,IIF,JJF,KKF,FAXIS,IEC,IE,SKIP_FCT,IEP,JEP,KEP,INDS(1:2,IAXIS:KAXIS),AX,ICEDG,JCEDG,LOHI
REAL(EB):: DXX(2),AREA_CF,XB_IB,DEL_EP,DEL_IBEDGE,XYZ1(IAXIS:KAXIS),XYZ2(IAXIS:KAXIS)

TYPE(CC_CUTEDGE_TYPE), POINTER :: CE

REAL(EB) CPUTIME,CPUTIME_START,CPUTIME_START_LOOP
CHARACTER(100) :: MSEGS_FILE
INTEGER :: ECOUNT

! Case of periodic test 103, return. No IBM interpolation needed as there are no immersed Bodies:
IF(PERIODIC_TEST==103 .OR. PERIODIC_TEST==11 .OR. PERIODIC_TEST==7) RETURN

! Total number of cell centered variables to be exchanges into external normal probe points of CFACES.
N_INT_CVARS = INT_P_IND + N_TRACKED_SPECIES

! Total number of cell centered variables to be exchanged in cut-cell interpolation:
N_INT_CCVARS= INT_WCEN_IND

IF(GET_CUTCELLS_VERBOSE) THEN
   WRITE(LU_SETCC,*) ' '; WRITE(LU_SETCC,'(A)') ' 5. In GET_CRTCFCC_INT_STENCILS, tasks to define IBM stencils:'
   CALL CPU_TIME(CPUTIME_START)
   CPUTIME_START_LOOP = CPUTIME_START
   WRITE(LU_SETCC,'(A)',advance='no') &
   ' - Into Mesh Loop: definition of Interpolation stencils for cut-cells and faces..'
ENDIF

! Then, second mesh loop:
IF( ASSOCIATED(X1FACEP)) NULLIFY(X1FACEP)
IF( ASSOCIATED(X2FACEP)) NULLIFY(X2FACEP)
IF( ASSOCIATED(X3FACEP)) NULLIFY(X3FACEP)
IF( ASSOCIATED(X2CELLP)) NULLIFY(X2CELLP)
IF( ASSOCIATED(X3CELLP)) NULLIFY(X3CELLP)
MESHES_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! Mesh sizes:
   NXB=IBAR
   NYB=JBAR
   NZB=KBAR

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + 1     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE ! Last internal cell index.
   ISTR     = ILO_FACE - NGUARD    ! Allocation start x arrays.
   IEND     = IHI_FACE + NGUARD    ! Allocation end x arrays.

   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + 1     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE ! Last internal cell index.
   JSTR     = JLO_FACE - NGUARD    ! Allocation start y arrays.
   JEND     = JHI_FACE + NGUARD    ! Allocation end y arrays.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + 1     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE ! Last internal cell index.
   KSTR     = KLO_FACE - NGUARD    ! Allocation start z arrays.
   KEND     = KHI_FACE + NGUARD    ! Allocation end z arrays.

   ! Define grid arrays for this mesh:
   ! Populate position and cell size arrays: Uniform grid implementation.
   ! X direction:
   ALLOCATE(DXCELL(ISTR:IEND)); DXCELL(ILO_CELL-1:IHI_CELL+1) = DX(ILO_CELL-1:IHI_CELL+1)
   DO IGC=2,NGUARD
      DXCELL(ILO_CELL-IGC)=DXCELL(ILO_CELL-IGC+1)
      DXCELL(IHI_CELL+IGC)=DXCELL(IHI_CELL+IGC-1)
   ENDDO
   ALLOCATE(DXFACE(ISTR:IEND)); DXFACE(ILO_FACE:IHI_FACE)= DXN(ILO_FACE:IHI_FACE)
   DO IGC=1,NGUARD
      DXFACE(ILO_FACE-IGC)=DXFACE(ILO_FACE-IGC+1)
      DXFACE(IHI_FACE+IGC)=DXFACE(ILO_FACE+IGC-1)
   ENDDO
   ALLOCATE(XCELL(ISTR:IEND));  XCELL = 1._EB/GEOMEPS ! Initialize huge.
   XCELL(ILO_CELL-1:IHI_CELL+1) = XC(ILO_CELL-1:IHI_CELL+1)
   DO IGC=2,NGUARD
      XCELL(ILO_CELL-IGC)=XCELL(ILO_CELL-IGC+1)-DXFACE(ILO_FACE-IGC+1)
      XCELL(IHI_CELL+IGC)=XCELL(IHI_CELL+IGC-1)+DXFACE(IHI_FACE+IGC-1)
   ENDDO
   ALLOCATE(XFACE(ISTR:IEND));  XFACE = 1._EB/GEOMEPS ! Initialize huge.
   XFACE(ILO_FACE:IHI_FACE) = X(ILO_FACE:IHI_FACE)
   DO IGC=1,NGUARD
      XFACE(ILO_FACE-IGC)=XFACE(ILO_FACE-IGC+1)-DXCELL(ILO_CELL-IGC)
      XFACE(IHI_FACE+IGC)=XFACE(IHI_FACE+IGC-1)+DXCELL(IHI_CELL+IGC)
   ENDDO

   ! Y direction:
   ALLOCATE(DYCELL(JSTR:JEND)); DYCELL(JLO_CELL-1:JHI_CELL+1)= DY(JLO_CELL-1:JHI_CELL+1)
   DO IGC=2,NGUARD
      DYCELL(JLO_CELL-IGC)=DYCELL(JLO_CELL-IGC+1)
      DYCELL(JHI_CELL+IGC)=DYCELL(JHI_CELL+IGC-1)
   ENDDO
   ALLOCATE(DYFACE(JSTR:JEND)); DYFACE(JLO_FACE:JHI_FACE)= DYN(JLO_FACE:JHI_FACE)
   DO IGC=1,NGUARD
      DYFACE(JLO_FACE-IGC)=DYFACE(JLO_FACE-IGC+1)
      DYFACE(JHI_FACE+IGC)=DYFACE(JHI_FACE+IGC-1)
   ENDDO
   ALLOCATE(YCELL(JSTR:JEND));  YCELL = 1._EB/GEOMEPS ! Initialize huge.
   YCELL(JLO_CELL-1:JHI_CELL+1) = YC(JLO_CELL-1:JHI_CELL+1)
   DO IGC=2,NGUARD
      YCELL(JLO_CELL-IGC)=YCELL(JLO_CELL-IGC+1)-DYFACE(JLO_FACE-IGC+1)
      YCELL(JHI_CELL+IGC)=YCELL(JHI_CELL+IGC-1)+DYFACE(JHI_FACE+IGC-1)
   ENDDO
   ALLOCATE(YFACE(JSTR:JEND));  YFACE = 1._EB/GEOMEPS ! Initialize huge.
   YFACE(JLO_FACE:JHI_FACE) = Y(JLO_FACE:JHI_FACE)
   DO IGC=1,NGUARD
      YFACE(JLO_FACE-IGC)=YFACE(JLO_FACE-IGC+1)-DYCELL(JLO_CELL-IGC)
      YFACE(JHI_FACE+IGC)=YFACE(JHI_FACE+IGC-1)+DYCELL(JHI_CELL+IGC)
   ENDDO

   ! Z direction:
   ALLOCATE(DZCELL(KSTR:KEND)); DZCELL(KLO_CELL-1:KHI_CELL+1)= DZ(KLO_CELL-1:KHI_CELL+1)
   DO IGC=2,NGUARD
      DZCELL(KLO_CELL-IGC)=DZCELL(KLO_CELL-IGC+1)
      DZCELL(KHI_CELL+IGC)=DZCELL(KHI_CELL+IGC-1)
   ENDDO
   ALLOCATE(DZFACE(KSTR:KEND)); DZFACE(KLO_FACE:KHI_FACE)= DZN(KLO_FACE:KHI_FACE)
   DO IGC=1,NGUARD
      DZFACE(KLO_FACE-IGC)=DZFACE(KLO_FACE-IGC+1)
      DZFACE(KHI_FACE+IGC)=DZFACE(KHI_FACE+IGC-1)
   ENDDO
   ALLOCATE(ZCELL(KSTR:KEND));  ZCELL = 1._EB/GEOMEPS ! Initialize huge.
   ZCELL(KLO_CELL-1:KHI_CELL+1) = ZC(KLO_CELL-1:KHI_CELL+1)
   DO IGC=2,NGUARD
      ZCELL(KLO_CELL-IGC)=ZCELL(KLO_CELL-IGC+1)-DZFACE(KLO_FACE-IGC+1)
      ZCELL(KHI_CELL+IGC)=ZCELL(KHI_CELL+IGC-1)+DZFACE(KHI_FACE+IGC-1)
   ENDDO
   ALLOCATE(ZFACE(KSTR:KEND));  ZFACE = 1._EB/GEOMEPS ! Initialize huge.
   ZFACE(KLO_FACE:KHI_FACE) = Z(KLO_FACE:KHI_FACE)
   DO IGC=1,NGUARD
      ZFACE(KLO_FACE-IGC)=ZFACE(KLO_FACE-IGC+1)-DZCELL(KLO_CELL-IGC)
      ZFACE(KHI_FACE+IGC)=ZFACE(KHI_FACE+IGC-1)+DZCELL(KHI_CELL+IGC)
   ENDDO

   ! Initialize CC_IDRC to CC_FGSC, this is to discard from the onset faces that can not be used in interpolation
   ! stencils (i.e. not fluid points). Faces allowed for interpolation stencils must be type CC_GASPHASE:
   FCVAR(:,:,:,CC_IDRC,IAXIS:KAXIS) = FCVAR(:,:,:,CC_FGSC,IAXIS:KAXIS)

   ! 1.:
   ! Loop by CUT_CELL, define interpolation stencils in Cartesian and cut
   ! cell centroids using the corresponding cells INBOUNDARY cut-faces:
   ! to be used for interpolation of H, etc.
   TESTVAR = CC_CGSC
   CUT_CELL_LOOP2 : DO ICC=1,MESHES(NM)%N_CUTCELL_MESH

      NCELL = MESHES(NM)%CUT_CELL(ICC)%NCELL
      IJK(IAXIS:KAXIS) = MESHES(NM)%CUT_CELL(ICC)%IJK(IAXIS:KAXIS)
      I = IJK(IAXIS); J = IJK(JAXIS); K = IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      MIN_DIST_VEL = DIST_THRES*MIN(DXCELL(I),DYCELL(J),DZCELL(K))

      ! First Cartesian centroid:
      XYZ(IAXIS:KAXIS) = (/ XCELL(I), YCELL(J), ZCELL(K) /)

      NPE_LIST_START = 0
      ALLOCATE(INT_NPE(LOW_IND:HIGH_IND,0:0,1:INT_N_EXT_PTS,0:CUT_CELL(ICC)%NCELL), &
               INT_IJK(IAXIS:KAXIS,(CUT_CELL(ICC)%NCELL+1)*DELTA_INT),                      &
               INT_COEF((CUT_CELL(ICC)%NCELL+1)*DELTA_INT),INT_NOUT(IAXIS:KAXIS,0:CUT_CELL(ICC)%NCELL))

      ! Now cut-cell volumes:
      ICELL_LOOP : DO ICELL=0,NCELL

         IF(ICELL > 0) XYZ(IAXIS:KAXIS)=MESHES(NM)%CUT_CELL(ICC)%XYZCEN(IAXIS:KAXIS,ICELL)

         ! Initialize closest inboundary point data:
         DISTANCE            = 1._EB / GEOMEPS
         LASTDOTNVEC         =-1._EB / GEOMEPS
         FOUND_POINT         = .FALSE.
         XYZ_PP(IAXIS:KAXIS) = 0._EB
         FOUND_INBFC(1:3)    = 0

         JCC_LOOP : DO JCC=1,NCELL

            ! Find closest point and Inboundary cut-face:
            NFC_CC = MESHES(NM)%CUT_CELL(ICC)%CCELEM(1,JCC)
            DO CCFC=1,NFC_CC

               ICFC = MESHES(NM)%CUT_CELL(ICC)%CCELEM(CCFC+1,JCC)
               IF ( MESHES(NM)%CUT_CELL(ICC)%FACE_LIST(1,ICFC) /= CC_FTYPE_CFINB) CYCLE

               ! Inboundary face number in CUT_FACE:
               INBFC     = MESHES(NM)%CUT_CELL(ICC)%FACE_LIST(4,ICFC)
               INBFC_LOC = MESHES(NM)%CUT_CELL(ICC)%FACE_LIST(5,ICFC)

               CALL GET_CLSPT_INBCF(NM,XYZ,INBFC,INBFC_LOC,XYZ_IP,DIST,FOUNDPT,INSEG)
               IF (FOUNDPT .AND. ((DIST-DISTANCE) < GEOMEPS)) THEN
                   IF (INSEG) THEN
                       BODTRI(1:2)  = CUT_FACE(INBFC)%BODTRI(1:2,INBFC_LOC)
                       ! normal vector to boundary surface triangle:
                       IBOD    = BODTRI(1)
                       IWSEL   = BODTRI(2)
                       NVEC(IAXIS:KAXIS) = GEOMETRY(IBOD)%FACES_NORMAL(IAXIS:KAXIS,IWSEL)
                       DV(IAXIS:KAXIS) = XYZ(IAXIS:KAXIS) - XYZ_IP(IAXIS:KAXIS)
                       NORM_DV = SQRT( DV(IAXIS)**2._EB + DV(JAXIS)**2._EB + DV(KAXIS)**2._EB )
                       IF(NORM_DV > GEOMEPS) THEN ! Point in segment not same as pt to interp to.
                          DV(IAXIS:KAXIS) = (1._EB / NORM_DV) * DV(IAXIS:KAXIS)
                          DOTNVEC = NVEC(IAXIS)*DV(IAXIS) + NVEC(JAXIS)*DV(JAXIS) + NVEC(KAXIS)*DV(KAXIS)
                          IF (DOTNVEC <= LASTDOTNVEC) CYCLE
                          LASTDOTNVEC = DOTNVEC
                       ENDIF
                   ENDIF
                   DISTANCE = DIST
                   XYZ_PP(IAXIS:KAXIS)   = XYZ_IP(IAXIS:KAXIS)
                   FOUND_INBFC(1:3) = (/ CC_FTYPE_CFINB, INBFC, INBFC_LOC /) ! Inbound cut-face in CUT_FACE.
                   FOUND_POINT = .TRUE.
               ENDIF

            ENDDO

            ! If point not found, all cut-faces boundary of the icc, jcc volume
            ! are GASPHASE. There must be a SOLID point in the boundary of the
            ! underlying Cartesian cell. this is the closest point:
            IF (.NOT.FOUND_POINT) THEN
                ! Search for for CUT_CELL(icc) vertex points or other solid points:
                CALL GET_CLOSEPT_CCVT(NM,XYZ,ICC,XYZ_IP,DIST,FOUNDPT,IFCPT,IFCPT_LOC)
                IF (FOUNDPT .AND. ((DIST-DISTANCE) < GEOMEPS)) THEN
                   DISTANCE = DIST
                   XYZ_PP(IAXIS:KAXIS)   = XYZ_IP(IAXIS:KAXIS)
                   FOUND_INBFC(1:3) = (/ CC_FTYPE_SVERT, IFCPT, IFCPT_LOC /) ! SOLID vertex in CUT_FACE.
                   FOUND_POINT = .TRUE.
                ENDIF
            ENDIF

         ENDDO JCC_LOOP

         IF (.NOT.FOUND_POINT .AND. DEBUG_CC_INTERPOLATION) THEN
            IF(ICELL==0) THEN
               WRITE(LU_ERR,*) 'CF: Havent found closest point CART CELL. NM,ICC=',NM,ICC
            ELSE
               WRITE(LU_ERR,*) 'CF: Havent found closest point CUT CELL. NM,ICC,JCC=',NM,ICC,JCC
            ENDIF
         ENDIF

         ! Here test if point in boundary and interpolation point coincide:
         IF (DISTANCE <= MIN_DIST_VEL) THEN

            INT_NOUT(IAXIS:KAXIS,ICELL) = 0._EB
            INT_XN(0:INT_N_EXT_PTS) = 0._EB
            INT_CN(0:INT_N_EXT_PTS) = 0._EB; INT_CN(0) = 1._EB ! Boundary point coefficient:
            VIND = 0
            DO EP=1,INT_N_EXT_PTS  ! External point for face IFACE
                NPE_LIST_COUNT = 0
                INT_NPE(LOW_IND,VIND,EP,ICELL)  = NPE_LIST_START
                INT_NPE(HIGH_IND,VIND,EP,ICELL) = NPE_LIST_COUNT
                NPE_LIST_START = NPE_LIST_START + NPE_LIST_COUNT
            ENDDO

         ELSE ! DISTANCE <= MIN_DIST_VEL

            ! After this loop we have the closest boundary point to xyz and the
            ! cut-face it belongs. We need to use the normal out of the face (or the
            ! vertex to xyz direction to find fluid points on the stencil:
            ! The fluid points are points that lay on the plane outside in the
            ! largest Cartesian component direction of the normal.
            DIR_FCT = 1._EB
            IF (FOUND_INBFC(1) == CC_FTYPE_CFINB) THEN ! closest point in INBOUNDARY cut-face.
                BODTRI(1:2) = CUT_FACE(FOUND_INBFC(2))%BODTRI(1:2,FOUND_INBFC(3))
                ! normal vector to boundary surface triangle:
                IBOD    = BODTRI(1)
                IWSEL   = BODTRI(2)
                NVEC(IAXIS:KAXIS) = GEOMETRY(IBOD)%FACES_NORMAL(IAXIS:KAXIS,IWSEL)
                DV(IAXIS:KAXIS) = XYZ(IAXIS:KAXIS) - XYZ_PP(IAXIS:KAXIS)
                DOTNVEC = NVEC(IAXIS)*DV(IAXIS) + NVEC(JAXIS)*DV(JAXIS) + NVEC(KAXIS)*DV(KAXIS)

                IF (DOTNVEC < 0._EB) DIR_FCT = -1._EB ! if normal to triangle has opposite dir change
                                                      ! search direction.
            ENDIF

            ! Versor to GASPHASE:
            IF (DIR_FCT > 0._EB) THEN ! Versor from boundary point to centroid
                P0(IAXIS:KAXIS) = XYZ_PP(IAXIS:KAXIS)
                P1(IAXIS:KAXIS) = XYZ(IAXIS:KAXIS)
            ELSE ! Viceversa
                P0(IAXIS:KAXIS) = XYZ(IAXIS:KAXIS)
                P1(IAXIS:KAXIS) = XYZ_PP(IAXIS:KAXIS)
            ENDIF
            DV(IAXIS:KAXIS)   = DIR_FCT * ( XYZ(IAXIS:KAXIS) - XYZ_PP(IAXIS:KAXIS) )
            NORM_DV           = SQRT( DV(IAXIS)**2._EB + DV(JAXIS)**2._EB + DV(KAXIS)**2._EB )
            DV(IAXIS:KAXIS) = (1._EB / NORM_DV) * DV(IAXIS:KAXIS) ! NOUT

            CALL GET_DELN(1.001_EB,DELN,DXCELL(I),DYCELL(J),DZCELL(K),NVEC=DV,CLOSE_PT=.TRUE.)

            ! Location of interpolation point XYZ(IAXIS:KAXIS) along the DV direction, origin in
            ! boundary point XYZ_PP(IAXIS:KAXIS):
            INT_NOUT(IAXIS:KAXIS,ICELL) = DV(IAXIS:KAXIS)
            INT_XN(0)               = DIR_FCT * NORM_DV
            INT_XN(1:INT_N_EXT_PTS) = 0._EB
            ! Initialize interpolation coefficients along the normal probe direction DV
            INT_CN(0) = 0._EB; ! Boundary point interpolation coefficient
            INT_CN(1:INT_N_EXT_PTS) = 0._EB;
            VIND = 0
            DO EP=1,INT_N_EXT_PTS  ! External point for face IFACE
               INT_XN(EP) = REAL(EP,EB)*DELN
               CALL GET_INTSTENCILS_EP(.FALSE.,VIND,XYZ_PP,INT_XN(EP),DV, &
                                       NPE_LIST_START,NPE_LIST_COUNT,INT_IJK,INT_COEF)
               ! Start position for interpolation stencil related to VIND=0, of external
               ! point EP related to cut-cell ICELL:
               INT_NPE(LOW_IND,VIND,EP,ICELL)  = NPE_LIST_START
               ! Number of stencil points on stencil for said cc.
               INT_NPE(HIGH_IND,VIND,EP,ICELL) = NPE_LIST_COUNT
               NPE_LIST_START = NPE_LIST_START + NPE_LIST_COUNT
            ENDDO

         ENDIF ! DISTANCE <= MIN_DIST_VEL

         CUT_CELL(ICC)%INT_XYZBF(IAXIS:KAXIS,ICELL) = XYZ_PP(IAXIS:KAXIS) ! xyz of boundary pt.
         CUT_CELL(ICC)%INT_INBFC(1:3,ICELL)         = FOUND_INBFC(1:3)  ! which INB cut-face bndry pt belongs to.
         CUT_CELL(ICC)%INT_NOUT(IAXIS:KAXIS,ICELL)  = INT_NOUT(IAXIS:KAXIS,ICELL)
         CUT_CELL(ICC)%INT_XN(0:INT_N_EXT_PTS,ICELL)= INT_XN(0:INT_N_EXT_PTS)
         CUT_CELL(ICC)%INT_CN(0:INT_N_EXT_PTS,ICELL)= INT_CN(0:INT_N_EXT_PTS)
         ! If size of CUT_CELL(ICC)%INT_IJK,DIM=2 is less than the size of INT_IJK, reallocate:
         SZ_1 = SIZE(CUT_CELL(ICC)%INT_IJK,DIM=2)
         SZ_2 = SIZE(INT_IJK,DIM=2)
         IF(SZ_2 > SZ_1) THEN
            ALLOCATE(INT_IJK_AUX(IAXIS:KAXIS,SZ_1),INT_COEF_AUX(1:SZ_1))
            INT_IJK_AUX(IAXIS:KAXIS,1:SZ_1)  = CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,1:SZ_1)
            INT_COEF_AUX(1:SZ_1)             = CUT_CELL(ICC)%INT_COEF(1:SZ_1)
            DEALLOCATE(CUT_CELL(ICC)%INT_IJK, CUT_CELL(ICC)%INT_COEF)
            ALLOCATE(CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,SZ_2)); CUT_CELL(ICC)%INT_IJK = CC_UNDEFINED
            ALLOCATE(CUT_CELL(ICC)%INT_COEF(1:SZ_2)); CUT_CELL(ICC)%INT_COEF = 0._EB
            CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,1:SZ_1)  = INT_IJK_AUX(IAXIS:KAXIS,1:SZ_1)
            CUT_CELL(ICC)%INT_COEF(1:SZ_1)             = INT_COEF_AUX(1:SZ_1)
            DEALLOCATE(INT_IJK_AUX,INT_COEF_AUX)
            DEALLOCATE(CUT_CELL(ICC)%INT_CCVARS,CUT_CELL(ICC)%INT_NOMIND)
            ALLOCATE(CUT_CELL(ICC)%INT_CCVARS(1:N_INT_CCVARS,SZ_2)); CUT_CELL(ICC)%INT_CCVARS=0._EB
            ALLOCATE(CUT_CELL(ICC)%INT_NOMIND(LOW_IND:HIGH_IND,SZ_2)); CUT_CELL(ICC)%INT_NOMIND = CC_UNDEFINED
         ENDIF
         VIND = 0
         DO EP=1,INT_N_EXT_PTS  ! External point for CELL ICELL
            INT_NPE_LO = INT_NPE(LOW_IND,VIND,EP,ICELL)
            INT_NPE_HI = INT_NPE(HIGH_IND,VIND,EP,ICELL)
            CUT_CELL(ICC)%INT_NPE(LOW_IND,VIND,EP,ICELL) = INT_NPE_LO
            CUT_CELL(ICC)%INT_NPE(HIGH_IND,VIND,EP,ICELL)= INT_NPE_HI
            DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
               CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,INPE)  = INT_IJK(IAXIS:KAXIS,INPE)
               CUT_CELL(ICC)%INT_COEF(INPE)             = INT_COEF(INPE)
            ENDDO
         ENDDO

      ENDDO ICELL_LOOP
      DEALLOCATE(INT_NPE,INT_IJK,INT_COEF,INT_NOUT)

   ENDDO CUT_CELL_LOOP2

   ! Compute stencils for RCEDGES, regular edges connecting cut and regular faces, and IBEDGES, solid edges next to cut-faces:
   RCEDGE_LOOP_1 : DO IEDGE=1,MESHES(NM)%CC_NRCEDGE
      ALLOCATE(CC_RCEDGE(IEDGE)%XB_IB(-2:2),CC_RCEDGE(IEDGE)%SURF_INDEX(-2:2),&
      CC_RCEDGE(IEDGE)%DUIDXJ(-2:2),CC_RCEDGE(IEDGE)%MU_DUIDXJ(-2:2))
      ALLOCATE(CC_RCEDGE(IEDGE)%INT_NPE(LOW_IND:HIGH_IND,0:KAXIS,1:INT_N_EXT_PTS,-2:2))
      CC_RCEDGE(IEDGE)%XB_IB(-2:2)      = 0._EB
      CC_RCEDGE(IEDGE)%SURF_INDEX(-2:2) = -1
      ! CC_RCEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(-2:2) = .FALSE. ! Process Orientation in double loop.
      ! CC_RCEDGE(IEDGE)%EDGE_IN_MESH(-2:2)             = .TRUE.  ! Always true for RCEDGES, no need to mesh_cc_exchange variables.
      CC_RCEDGE(IEDGE)%INT_NPE                          = 0       ! Required to avoid segfault in comm.

      IE  = MESHES(NM)%CC_RCEDGE(IEDGE)%IE
      II  = EDGE(IE)%I
      JJ  = EDGE(IE)%J
      KK  = EDGE(IE)%K
      IEC = EDGE(IE)%AXIS

      ! First: Loop over all possible face orientations of edge to define XB_IB, SURF_INDEX, PROCESS_EDGE_ORIENTATION:
      ORIENTATION_LOOP_RC_1: DO IS=1,3
         IF (IS==IEC) CYCLE ORIENTATION_LOOP_RC_1
         SIGN_LOOP_RC_1: DO I_SGN=-1,1,2

            ! Determine Index_Coordinate_Direction
            ! IEC=1, ICD=1 refers to DWDY; ICD=2 refers to DVDZ
            ! IEC=2, ICD=1 refers to DUDZ; ICD=2 refers to DWDX
            ! IEC=3, ICD=1 refers to DVDX; ICD=2 refers to DUDY

            IF (IS>IEC) ICD = IS-IEC
            IF (IS<IEC) ICD = IS-IEC+3
            ICD_SGN = I_SGN * ICD

            ! With ICD_SGN check if face:
            ! IBEDGE IEC=IAXIS => ICD_SGN=-2 => FACE  low Z normal to JAXIS.
            !                     ICD_SGN=-1 => FACE  low Y normal to KAXIS.
            !                     ICD_SGN= 1 => FACE high Y normal to KAXIS.
            !                     ICD_SGN= 2 => FACE high Z normal to JAXIS.
            ! IBEDGE IEC=JAXIS => ICD_SGN=-2 => FACE  low X normal to KAXIS.
            !                     ICD_SGN=-1 => FACE  low Z normal to IAXIS.
            !                     ICD_SGN= 1 => FACE high Z normal to IAXIS.
            !                     ICD_SGN= 2 => FACE high X normal to KAXIS.
            ! IBEDGE IEC=KAXIS => ICD_SGN=-2 => FACE  low Y normal to IAXIS.
            !                     ICD_SGN=-1 => FACE  low X normal to JAXIS.
            !                     ICD_SGN= 1 => FACE high X normal to JAXIS.
            !                     ICD_SGN= 2 => FACE high Y normal to IAXIS.
            ! is GASPHASE cut-face.
            SELECT CASE(IEC)
               CASE(IAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE( 1); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=KAXIS
                     CASE( 2); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=JAXIS
                  END SELECT
                  ! Compute XB_IB: For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DY(JJF); DXX(2)  = DZ(KKF); DEL_IBEDGE = DX(IIF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(IAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(IAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                  ENDIF
                  IF (FAXIS==JAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(2) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           ! Low side cut-cell: Load first cut-face SURF_INDEX:
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF,JJF+1,KKF,CC_IDCF)>0) THEN
                           ! High side cut-cell: Load first cut-face SURF_INDEX:
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF+1,KKF,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ELSE ! IF(FAXIS==KAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(1) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF,JJF,KKF+1,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF+1,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ENDIF

               CASE(JAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE( 1); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=IAXIS
                     CASE( 2); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                  END SELECT
                  ! Compute XB_IB: For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DZ(KKF); DXX(2)  = DX(IIF); DEL_IBEDGE = DY(JJF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(JAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(JAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                  ENDIF
                  IF (FAXIS==KAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(2) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF,JJF,KKF+1,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF+1,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ELSE ! IF(FAXIS==IAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(1) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF+1,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF+1,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ENDIF

               CASE(KAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 1); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 2); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=IAXIS
                  END SELECT
                  ! Compute XB_IB: For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DX(IIF); DXX(2)  = DY(JJF); DEL_IBEDGE = DZ(KKF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(KAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(KAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                  ENDIF
                  IF (FAXIS==IAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(2) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF+1,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF+1,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ELSE ! IF(FAXIS==JAXIS) THEN
                     CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN) = DXX(1) ! Twice Distance to velocity collocation point.
                     IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)==CC_CUTCFE) THEN
                        CC_RCEDGE(IEDGE)%XB_IB(ICD_SGN)=(AREA_CF/DEL_IBEDGE)
                        ! SURF_INDEX:
                        IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                        ELSEIF(CCVAR(IIF,JJF+1,KKF,CC_IDCF)>0) THEN
                           CC_RCEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF+1,KKF,CC_IDCF))%SURF_INDEX(1)
                        ENDIF
                     ENDIF
                  ENDIF

             END SELECT

          ENDDO SIGN_LOOP_RC_1
       ENDDO ORIENTATION_LOOP_RC_1

   ENDDO RCEDGE_LOOP_1

   ! Dummy allocation for now:
   IBEDGE_LOOP1 : DO IEDGE=1,MESHES(NM)%CC_NIBEDGE

      ALLOCATE(CC_IBEDGE(IEDGE)%XB_IB(-2:2),CC_IBEDGE(IEDGE)%SURF_INDEX(-2:2),&
               CC_IBEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(-2:2),CC_IBEDGE(IEDGE)%EDGE_IN_MESH(-2:2),&
               CC_IBEDGE(IEDGE)%SIDE_IN_GEOM(-2:2))
      ALLOCATE(CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND:HIGH_IND,0:KAXIS,1:INT_N_EXT_PTS,-2:2))
      CC_IBEDGE(IEDGE)%XB_IB(-2:2)      = 0._EB
      CC_IBEDGE(IEDGE)%SURF_INDEX(-2:2) = 0
      CC_IBEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(-2:2) = .FALSE. ! Process Orientation in double loop.
      CC_IBEDGE(IEDGE)%EDGE_IN_MESH(-2:2)             = .FALSE. ! If true, no need to mesh_cc_exchange variables.
      CC_IBEDGE(IEDGE)%SIDE_IN_GEOM(-2:2)             = .TRUE.  ! If true, this side of Cartesian edge looks inside the geometry.
                                                                ! else there is a boundary edge on its position for this side.
      CC_IBEDGE(IEDGE)%INT_NPE                        = 0 ! Required to avoid segfault in comm.

      IE  = MESHES(NM)%CC_IBEDGE(IEDGE)%IE
      II  = EDGE(IE)%I
      JJ  = EDGE(IE)%J
      KK  = EDGE(IE)%K
      IEC = EDGE(IE)%AXIS

      ! Edge nodes location:
      XYZ1 = (/ X(II), Y(JJ), Z(KK) /); XYZ2 = XYZ1
      SELECT CASE(IEC)
         CASE(IAXIS); XYZ1(IAXIS) = X(II-1)
         CASE(JAXIS); XYZ1(JAXIS) = Y(JJ-1)
         CASE(KAXIS); XYZ1(KAXIS) = Z(KK-1)
      END SELECT

      ! First: Loop over all possible face orientations of edge to define XB_IB, SURF_INDEX, PROCESS_EDGE_ORIENTATION:
      ORIENTATION_LOOP_1: DO IS=1,3
         IF (IS==IEC) CYCLE ORIENTATION_LOOP_1
         SIGN_LOOP_1: DO I_SGN=-1,1,2

            ! Determine Index_Coordinate_Direction
            ! IEC=1, ICD=1 refers to DWDY; ICD=2 refers to DVDZ
            ! IEC=2, ICD=1 refers to DUDZ; ICD=2 refers to DWDX
            ! IEC=3, ICD=1 refers to DVDX; ICD=2 refers to DUDY

            IF (IS>IEC) ICD = IS-IEC
            IF (IS<IEC) ICD = IS-IEC+3
            ICD_SGN = I_SGN * ICD

            ! With ICD_SGN check if face:
            ! IBEDGE IEC=IAXIS => ICD_SGN=-2 => FACE  low Z normal to JAXIS.
            !                     ICD_SGN=-1 => FACE  low Y normal to KAXIS.
            !                     ICD_SGN= 1 => FACE high Y normal to KAXIS.
            !                     ICD_SGN= 2 => FACE high Z normal to JAXIS.
            ! IBEDGE IEC=JAXIS => ICD_SGN=-2 => FACE  low X normal to KAXIS.
            !                     ICD_SGN=-1 => FACE  low Z normal to IAXIS.
            !                     ICD_SGN= 1 => FACE high Z normal to IAXIS.
            !                     ICD_SGN= 2 => FACE high X normal to KAXIS.
            ! IBEDGE IEC=KAXIS => ICD_SGN=-2 => FACE  low Y normal to IAXIS.
            !                     ICD_SGN=-1 => FACE  low X normal to JAXIS.
            !                     ICD_SGN= 1 => FACE high X normal to JAXIS.
            !                     ICD_SGN= 2 => FACE high Y normal to IAXIS.
            ! is GASPHASE cut-face.
            SELECT CASE(IEC)
               CASE(IAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE( 1); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=KAXIS
                     CASE( 2); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=JAXIS
                  END SELECT
                  ! Drop if face is not type CUTCFE:
                  IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)/=CC_CUTCFE) CYCLE SIGN_LOOP_1
                  ! Compute XB_IB, SURF_INDEX:
                  ! For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DY(JJF); DXX(2)  = DZ(KKF); DEL_IBEDGE = DX(IIF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(IAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(IAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                     ICEDG = FCVAR(IIF,JJF,KKF,CC_IDCE,FAXIS)
                     IF(ICEDG>0) THEN
                        CE => CUT_EDGE(ICEDG)
                        JEDG_LOOP_1 : DO JCEDG=1,CE%NEDGE
                           DO LOHI=1,2; DO AX=IAXIS,KAXIS
                           IF( ABS(XYZ1(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS .AND. &
                               ABS(XYZ2(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS ) CYCLE JEDG_LOOP_1
                           ENDDO; ENDDO
                           CC_IBEDGE(IEDGE)%SIDE_IN_GEOM(ICD_SGN) = .FALSE. ! This side of Cartesian edge looks into the gasphase.
                        ENDDO JEDG_LOOP_1
                     ENDIF
                  ENDIF
                  IF (FAXIS==JAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(2)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        ! Low side cut-cell: Load first cut-face SURF_INDEX:
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF,JJF+1,KKF,CC_IDCF)>0) THEN
                        ! High side cut-cell: Load first cut-face SURF_INDEX:
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF+1,KKF,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ELSE ! IF(FAXIS==KAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(1)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF,JJF,KKF+1,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF+1,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ENDIF

                  ! Now search where the EP external stress edge will be defined:
                  XB_IB = CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN)
                  SKIP_FCT = 1
                  IF (FAXIS==JAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     JEP=JJ; KEP=KK+SKIP_FCT*I_SGN
                  ELSE ! IF(FAXIS==KAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     JEP=JJ+SKIP_FCT*I_SGN; KEP=KK
                  ENDIF
                  IF( JEP<=JBAR .AND. JEP>=0 .AND. KEP<=KBAR .AND. KEP>=0 ) CC_IBEDGE(IEDGE)%EDGE_IN_MESH(ICD_SGN) = .TRUE.

               CASE(JAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE( 1); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=IAXIS
                     CASE( 2); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                  END SELECT
                  IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)/=CC_CUTCFE) CYCLE SIGN_LOOP_1
                  ! Compute XB_IB, SURF_INDEX:
                  ! For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DZ(KKF); DXX(2)  = DX(IIF); DEL_IBEDGE = DY(JJF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(JAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(JAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                     ICEDG = FCVAR(IIF,JJF,KKF,CC_IDCE,FAXIS)
                     IF(ICEDG>0) THEN
                        CE => CUT_EDGE(ICEDG)
                        JEDG_LOOP_2 : DO JCEDG=1,CE%NEDGE
                           DO LOHI=1,2; DO AX=IAXIS,KAXIS
                           IF( ABS(XYZ1(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS .AND. &
                               ABS(XYZ2(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS ) CYCLE JEDG_LOOP_2
                           ENDDO; ENDDO
                           CC_IBEDGE(IEDGE)%SIDE_IN_GEOM(ICD_SGN) = .FALSE. ! This side of Cartesian edge looks into the gasphase.
                        ENDDO JEDG_LOOP_2
                     ENDIF
                  ENDIF
                  IF (FAXIS==KAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(2)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF,JJF,KKF+1,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF+1,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ELSE ! IF(FAXIS==IAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(1)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF+1,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF+1,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ENDIF

                  ! Now search where the EP external stress edge will be defined:
                  XB_IB = CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN)
                  SKIP_FCT = 1
                  IF (FAXIS==KAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; KEP=KK
                  ELSE ! IF(FAXIS==IAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; KEP=KK+SKIP_FCT*I_SGN
                  ENDIF
                  IF( IEP<=IBAR .AND. IEP>=0 .AND. KEP<=KBAR .AND. KEP>=0 ) CC_IBEDGE(IEDGE)%EDGE_IN_MESH(ICD_SGN) = .TRUE.

               CASE(KAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 1); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 2); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=IAXIS
                  END SELECT
                  IF(FCVAR(IIF,JJF,KKF,CC_FGSC,FAXIS)/=CC_CUTCFE) CYCLE SIGN_LOOP_1
                  ! Compute XB_IB, SURF_INDEX:
                  ! For XB_IB we use the sum of gas cut-faces in the face and compare it with the face AREA:
                  AREA_CF = 0._EB; DXX(1)  = DX(IIF); DXX(2)  = DY(JJF); DEL_IBEDGE = DZ(KKF)
                  ICF = FCVAR(IIF,JJF,KKF,CC_IDCF,FAXIS)
                  IF(ICF>0) THEN
                     AREA_CF = SUM(CUT_FACE(ICF)%AREA(1:CUT_FACE(ICF)%NFACE))
                     DEL_IBEDGE = ABS(MAXVAL(CUT_FACE(ICF)%XYZVERT(KAXIS,1:CUT_FACE(ICF)%NVERT)) - &
                                      MINVAL(CUT_FACE(ICF)%XYZVERT(KAXIS,1:CUT_FACE(ICF)%NVERT))) + TWENTY_EPSILON_EB
                     ICEDG = FCVAR(IIF,JJF,KKF,CC_IDCE,FAXIS)
                     IF(ICEDG>0) THEN
                        CE => CUT_EDGE(ICEDG)
                        JEDG_LOOP_3 : DO JCEDG=1,CE%NEDGE
                           DO LOHI=1,2; DO AX=IAXIS,KAXIS
                           IF( ABS(XYZ1(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS .AND. &
                               ABS(XYZ2(AX)-CE%XYZVERT(AX,CE%CEELEM(LOHI,JCEDG)))>GEOMEPS ) CYCLE JEDG_LOOP_3
                           ENDDO; ENDDO
                           CC_IBEDGE(IEDGE)%SIDE_IN_GEOM(ICD_SGN) = .FALSE. ! This side of Cartesian edge looks into the gasphase.
                        ENDDO JEDG_LOOP_3
                     ENDIF
                  ENDIF
                  IF (FAXIS==IAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(2)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF+1,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF+1,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ELSE ! IF(FAXIS==JAXIS) THEN
                     ! XB_IB:
                     CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN) = -(DXX(1)-AREA_CF/DEL_IBEDGE) !-ve dist Bound to IBEDGE opposed to normal.
                     ! SURF_INDEX:
                     IF (CCVAR(IIF,JJF,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF,KKF,CC_IDCF))%SURF_INDEX(1)
                     ELSEIF(CCVAR(IIF,JJF+1,KKF,CC_IDCF)>0) THEN
                        CC_IBEDGE(IEDGE)%SURF_INDEX(ICD_SGN) = CUT_FACE(CCVAR(IIF,JJF+1,KKF,CC_IDCF))%SURF_INDEX(1)
                     ENDIF
                  ENDIF

                  ! Now search where the EP external stress edge will be defined:
                  XB_IB = CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN)
                  SKIP_FCT = 1
                  IF (FAXIS==IAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ+SKIP_FCT*I_SGN
                  ELSE ! IF(FAXIS==JAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; JEP=JJ
                  ENDIF
                  IF( IEP<=IBAR .AND. IEP>=0 .AND. JEP<=JBAR .AND. JEP>=0 ) CC_IBEDGE(IEDGE)%EDGE_IN_MESH(ICD_SGN) = .TRUE.

            END SELECT

            CC_IBEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(ICD_SGN) = .TRUE.

         ENDDO SIGN_LOOP_1
      ENDDO ORIENTATION_LOOP_1

      ! If the edge orientation is not EDGE_IN_MESH, velocity, MU data for EP is communicated:
      NPE_LIST_START = 0
      ALLOCATE(INT_NPE(LOW_IND:HIGH_IND,0:KAXIS,1:INT_N_EXT_PTS,-2:2),INT_IJK(IAXIS:KAXIS,32)); INT_NPE = 0; INT_IJK = 0;
      ALLOCATE(INT_DCOEF(32,1)); INT_DCOEF = 0._EB
      ! First cell centered Variable MU:
      EP   = 1; N_CVAR_START = NPE_LIST_START
      ORIENTATION_LOOP_2: DO IS=1,3
         IF (IS==IEC) CYCLE ORIENTATION_LOOP_2
         SIGN_LOOP_2: DO I_SGN=-1,1,2
            IF (IS>IEC) ICD = IS-IEC
            IF (IS<IEC) ICD = IS-IEC+3
            ICD_SGN = I_SGN * ICD

            IF(.NOT.CC_IBEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(ICD_SGN)) CYCLE SIGN_LOOP_2
            IF(CC_IBEDGE(IEDGE)%EDGE_IN_MESH(ICD_SGN)) CYCLE SIGN_LOOP_2

            XB_IB = CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN)

            SELECT CASE(IEC)
               CASE(IAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE( 1); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=KAXIS
                     CASE( 2); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=JAXIS
                  END SELECT
                  DXX(1)  = DY(JJF); DXX(2)  = DZ(KKF)
                  SKIP_FCT = 1
                  IF (FAXIS==JAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ; KEP=KK+SKIP_FCT*I_SGN
                  ELSE ! IF(FAXIS==KAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ+SKIP_FCT*I_SGN; KEP=KK
                  ENDIF
                  ! Add I,J,K locations of cells:
                  INDS(1:2,IAXIS) = (/0, 0/)
                  INDS(1:2,JAXIS) = (/0, 1/)
                  INDS(1:2,KAXIS) = (/0, 1/)

               CASE(JAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE( 1); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=IAXIS
                     CASE( 2); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                  END SELECT
                  DXX(1)  = DZ(KKF); DXX(2)  = DX(IIF)
                  SKIP_FCT = 1
                  IF (FAXIS==KAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; JEP=JJ; KEP=KK
                  ELSE ! IF(FAXIS==IAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ; KEP=KK+SKIP_FCT*I_SGN
                  ENDIF
                  ! Add I,J,K locations of cells:
                  INDS(1:2,IAXIS) = (/0, 1/)
                  INDS(1:2,JAXIS) = (/0, 0/)
                  INDS(1:2,KAXIS) = (/0, 1/)

               CASE(KAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 1); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 2); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=IAXIS
                  END SELECT
                  DXX(1)  = DX(IIF); DXX(2)  = DY(JJF)
                  SKIP_FCT = 1
                  IF (FAXIS==IAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ+SKIP_FCT*I_SGN; KEP=KK
                  ELSE ! IF(FAXIS==JAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; JEP=JJ; KEP=KK
                  ENDIF
                  ! Add I,J,K locations of cells:
                  INDS(1:2,IAXIS) = (/0, 1/)
                  INDS(1:2,JAXIS) = (/0, 1/)
                  INDS(1:2,KAXIS) = (/0, 0/)

            END SELECT

            ! ADD all
            VIND = 0; NPE_LIST_COUNT  = 0
            DO K=INDS(1,KAXIS),INDS(2,KAXIS)
               DO J=INDS(1,JAXIS),INDS(2,JAXIS)
                  DO I=INDS(1,IAXIS),INDS(2,IAXIS)
                     ! IF(CELL(CELL_INDEX(IEP+I,JEP+J,KEP+K))%SOLID) CYCLE ! Cycle solid cells. Can't use it here as is (overrun).
                     IF(CCVAR(IEP+I,JEP+J,KEP+K,CC_CGSC)==CC_SOLID) CYCLE
                     NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                     INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP+I,JEP+J,KEP+K/)
                  ENDDO
               ENDDO
            ENDDO
            ! Start position and number of points for cell centered vars related to EP edge of ICD_SGN orientation:
            INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
            INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
            NPE_LIST_START = NPE_LIST_START + NPE_LIST_COUNT

         ENDDO SIGN_LOOP_2
      ENDDO ORIENTATION_LOOP_2
      ! Number of cell centered stencil points:
      N_CVAR_COUNT = NPE_LIST_START

      ! Now add Face Variables for the two directions normal to IEC:
      N_FVAR_START = N_CVAR_START + N_CVAR_COUNT
      ORIENTATION_LOOP_3: DO IS=1,3
         IF (IS==IEC) CYCLE ORIENTATION_LOOP_3
         SIGN_LOOP_3: DO I_SGN=-1,1,2
            IF (IS>IEC) ICD = IS-IEC
            IF (IS<IEC) ICD = IS-IEC+3
            ICD_SGN = I_SGN * ICD

            IF(.NOT.CC_IBEDGE(IEDGE)%PROCESS_EDGE_ORIENTATION(ICD_SGN)) CYCLE SIGN_LOOP_3
            IF(CC_IBEDGE(IEDGE)%EDGE_IN_MESH(ICD_SGN)) CYCLE SIGN_LOOP_3

            XB_IB = CC_IBEDGE(IEDGE)%XB_IB(ICD_SGN)

            SELECT CASE(IEC)
               CASE(IAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE( 1); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=KAXIS
                     CASE( 2); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=JAXIS
                  END SELECT
                  DXX(1)  = DY(JJF); DXX(2)  = DZ(KKF)
                  SKIP_FCT = 1
                  IF (FAXIS==JAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ; KEP=KK+SKIP_FCT*I_SGN
                  ELSE ! IF(FAXIS==KAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ+SKIP_FCT*I_SGN; KEP=KK
                  ENDIF

                  ! V velocities in EP for KEP,KEP+1:
                  VIND = JAXIS; NPE_LIST_COUNT = 0
                  DO K=0,1
                    NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                    INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP  ,JEP  ,KEP+K/)
                    INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*K-1,EB)/DXX(2)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

                  ! W Velocities in EP for JEP,JEP+1:
                  VIND = KAXIS; NPE_LIST_COUNT = 0
                  DO J=0,1
                     NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                     INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP  ,JEP+J,KEP  /)
                     INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*J-1,EB)/DXX(1)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

               CASE(JAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE( 1); IIF=II  ; JJF=JJ  ; KKF=KK+1; FAXIS=IAXIS
                     CASE( 2); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=KAXIS
                  END SELECT
                  DXX(1)  = DZ(KKF); DXX(2)  = DX(IIF)
                  SKIP_FCT = 1
                  IF (FAXIS==KAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; JEP=JJ; KEP=KK
                  ELSE ! IF(FAXIS==IAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ; KEP=KK+SKIP_FCT*I_SGN
                  ENDIF

                  ! W Velocities in EP for IEP,IEP+1:
                  VIND = KAXIS; NPE_LIST_COUNT = 0
                  DO I=0,1
                     NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                     INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP+I,JEP  ,KEP  /)
                     INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*I-1,EB)/DXX(2)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

                  ! U Velocities in EP for KEP,KEP+1:
                  VIND = IAXIS; NPE_LIST_COUNT = 0
                  DO K=0,1
                     NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                     INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP  ,JEP  ,KEP+K/)
                     INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*K-1,EB)/DXX(1)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

               CASE(KAXIS)
                  ! Define Face indexes and normal axis FAXIS.
                  SELECT CASE(ICD_SGN)
                     CASE(-2); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=IAXIS
                     CASE(-1); IIF=II  ; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 1); IIF=II+1; JJF=JJ  ; KKF=KK  ; FAXIS=JAXIS
                     CASE( 2); IIF=II  ; JJF=JJ+1; KKF=KK  ; FAXIS=IAXIS
                  END SELECT
                  DXX(1)  = DX(IIF); DXX(2)  = DY(JJF)
                  SKIP_FCT = 1
                  IF (FAXIS==IAXIS) THEN
                     DEL_EP = DXX(2) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(2) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II; JEP=JJ+SKIP_FCT*I_SGN; KEP=KK
                  ELSE ! IF(FAXIS==JAXIS) THEN
                     DEL_EP = DXX(1) - ABS(XB_IB)
                     IF( DEL_EP < THRES_FCT_EP*DXX(1) ) SKIP_FCT = 2 ! Pick next EP point +2*I_SGN
                     IEP=II+SKIP_FCT*I_SGN; JEP=JJ; KEP=KK
                  ENDIF

                  ! U Velocities in EP for JEP,JEP+1:
                  VIND = IAXIS; NPE_LIST_COUNT = 0
                  DO J=0,1
                     NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                     INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP  ,JEP+J,KEP  /)
                     INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*J-1,EB)/DXX(2)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

                  ! V velocities in EP for IEP,IEP+1:
                  VIND = JAXIS; NPE_LIST_COUNT = 0
                  DO I=0,1
                    NPE_LIST_COUNT = NPE_LIST_COUNT + 1
                    INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/IEP+I,JEP  ,KEP  /)
                    INT_DCOEF(NPE_LIST_START+NPE_LIST_COUNT,1) = REAL(2*I-1,EB)/DXX(1)
                  ENDDO
                  INT_NPE(LOW_IND,VIND,EP,ICD_SGN)  = NPE_LIST_START
                  INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_LIST_COUNT
                  NPE_LIST_START                    = NPE_LIST_START + NPE_LIST_COUNT

            END SELECT

         ENDDO SIGN_LOOP_3
      ENDDO ORIENTATION_LOOP_3
      N_FVAR_COUNT = NPE_LIST_START - N_FVAR_START

      IF (NPE_LIST_START > 0) THEN
         ! Allocate INT_IJK, INT_CVARS, INT_FVARS:
         ALLOCATE(CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,NPE_LIST_START))
         ALLOCATE(CC_IBEDGE(IEDGE)%INT_CVARS(1:N_INT_EP_CVARS,N_CVAR_START+1:N_CVAR_START+N_CVAR_COUNT))
         ALLOCATE(CC_IBEDGE(IEDGE)%INT_FVARS(1:N_INT_EP_FVARS,N_FVAR_START+1:N_FVAR_START+N_FVAR_COUNT))
         CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND:HIGH_IND,0:KAXIS,1:INT_N_EXT_PTS,-2:2) = &
                           INT_NPE(LOW_IND:HIGH_IND,0:KAXIS,1:INT_N_EXT_PTS,-2:2)
         CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,1:NPE_LIST_START) = INT_IJK(IAXIS:KAXIS,1:NPE_LIST_START)
         CC_IBEDGE(IEDGE)%INT_CVARS = 0._EB; CC_IBEDGE(IEDGE)%INT_FVARS = 0._EB
         ALLOCATE(CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,NPE_LIST_START)); CC_IBEDGE(IEDGE)%INT_NOMIND = CC_UNDEFINED
         ALLOCATE(CC_IBEDGE(IEDGE)%INT_DCOEF(NPE_LIST_START,1));
         CC_IBEDGE(IEDGE)%INT_DCOEF(1:NPE_LIST_START,1) = INT_DCOEF(1:NPE_LIST_START,1)
      ENDIF

      DEALLOCATE(INT_NPE,INT_IJK,INT_DCOEF)

   ENDDO IBEDGE_LOOP1

   ! Up to this point we have cut-cells, regular, immersed edges.
   ! 1. CUT_CELL
   ! 2. CC_RCEDGE
   ! 3. CC_IBEDGE

   DO NOM=1,NMESHES
      ! Also considers the case NOM==NM as a regular case.
      ! Face Variables:
      ALLOCATE(MESHES(NM)%OMESH(NOM)%IIO_FC_R(DELTA_FC))
      ALLOCATE(MESHES(NM)%OMESH(NOM)%JJO_FC_R(DELTA_FC))
      ALLOCATE(MESHES(NM)%OMESH(NOM)%KKO_FC_R(DELTA_FC))
      ALLOCATE(MESHES(NM)%OMESH(NOM)%AXS_FC_R(DELTA_FC))
      ! Cell Variables:
      ALLOCATE(MESHES(NM)%OMESH(NOM)%IIO_CC_R(DELTA_FC))
      ALLOCATE(MESHES(NM)%OMESH(NOM)%JJO_CC_R(DELTA_FC))
      ALLOCATE(MESHES(NM)%OMESH(NOM)%KKO_CC_R(DELTA_FC))
   ENDDO

   ! Figure out which Regular face locations for this mesh are required for interpolation:
   ALLOCATE(IJKFACE2(LOW_IND:HIGH_IND,ISTR:IEND,JSTR:JEND,KSTR:KEND,IAXIS:KAXIS)); IJKFACE2 = CC_UNDEFINED

   ! Figure out which other meshes this mesh will receive face centered variables from:
   ! 1. RCEDGES:
   DO IEDGE=1,MESHES(NM)%CC_NRCEDGE
      DO EP=1,INT_N_EXT_PTS  ! External point for IEDGE
         DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
            INT_NPE_LO = CC_RCEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,0)
            INT_NPE_HI = CC_RCEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0); IF (INT_NPE_HI<1) CYCLE
            X1AXIS = VIND
            ALLOCATE(EP_TAG(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)); EP_TAG(:)=.FALSE.
            DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
               I = CC_RCEDGE(IEDGE)%INT_IJK(IAXIS,INPE)
               J = CC_RCEDGE(IEDGE)%INT_IJK(JAXIS,INPE)
               K = CC_RCEDGE(IEDGE)%INT_IJK(KAXIS,INPE)
               SELECT CASE(X1AXIS)
               CASE(IAXIS)
                  IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                     FLGX = (I >= ILO_FACE) .AND. (I <= IHI_FACE)
                     FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                     FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                     INNM = FLGX .AND. FLGY .AND. FLGZ
                     IF (INNM) THEN
                        NOM=NM; IIO=I; JJO=J; KKO=K
                     ELSE
                        CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XFACE(I),YCELL(J),ZCELL(K),NOM,IIO,JJO,KKO)
                     ENDIF
                     IF(NOM==0) EP_TAG(INPE) = .TRUE.
                     CALL ASSIGN_TO_FC_R
                  ENDIF
               CASE(JAXIS)
                  IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                     FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                     FLGY = (J >= JLO_FACE) .AND. (J <= JHI_FACE)
                     FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                     INNM = FLGX .AND. FLGY .AND. FLGZ
                     IF (INNM) THEN
                        NOM=NM; IIO=I; JJO=J; KKO=K
                     ELSE
                        CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XCELL(I),YFACE(J),ZCELL(K),NOM,IIO,JJO,KKO)
                     ENDIF
                     IF(NOM==0) EP_TAG(INPE) = .TRUE.
                     CALL ASSIGN_TO_FC_R
                  ENDIF
               CASE(KAXIS)
                  IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                     FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                     FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                     FLGZ = (K >= KLO_FACE) .AND. (K <= KHI_FACE)
                     INNM = FLGX .AND. FLGY .AND. FLGZ
                     IF (INNM) THEN
                        NOM=NM; IIO=I; JJO=J; KKO=K
                     ELSE
                        CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XCELL(I),YCELL(J),ZFACE(K),NOM,IIO,JJO,KKO)
                     ENDIF
                     IF(NOM==0) EP_TAG(INPE) = .TRUE.
                     CALL ASSIGN_TO_FC_R
                  ENDIF
               END SELECT
               CC_RCEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE) = IJKFACE2(LOW_IND:HIGH_IND,I,J,K,X1AXIS)
            ENDDO
            ! Now restrict count on cut-face :
            IF(ANY(EP_TAG .EQV. .TRUE.)) CALL RESTRICT_EP(CC_ETYPE_RCGAS)
            DEALLOCATE(EP_TAG)
            ! Compute derivative coefficients.
            CALL COMPUTE_DCOEF(CC_ETYPE_RCGAS)
         ENDDO
      ENDDO
   ENDDO
   ! 2. IBEDGES:
   DO IEDGE=1,MESHES(NM)%CC_NIBEDGE
      DO ICD_SGN=-2,2
         IF(ICD_SGN==0) CYCLE
         DO EP=1,INT_N_EXT_PTS  ! External point for IEDGE
            DO VIND=IAXIS,KAXIS ! Velocity component U, V or W for external point EP
               INT_NPE_HI = CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN); IF (INT_NPE_HI<1) CYCLE
               INT_NPE_LO = CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
               X1AXIS = VIND
               ALLOCATE(EP_TAG(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)); EP_TAG(:)=.FALSE.
               DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
                  I = CC_IBEDGE(IEDGE)%INT_IJK(IAXIS,INPE)
                  J = CC_IBEDGE(IEDGE)%INT_IJK(JAXIS,INPE)
                  K = CC_IBEDGE(IEDGE)%INT_IJK(KAXIS,INPE)
                  SELECT CASE(X1AXIS)
                  CASE(IAXIS)
                     IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                        FLGX = (I >= ILO_FACE) .AND. (I <= IHI_FACE)
                        FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                        FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                        INNM = FLGX .AND. FLGY .AND. FLGZ
                        IF (INNM) THEN
                           NOM=NM; IIO=I; JJO=J; KKO=K
                        ELSE
                           CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XFACE(I),YCELL(J),ZCELL(K),NOM,IIO,JJO,KKO)
                        ENDIF
                        IF(NOM==0) EP_TAG(INPE) = .TRUE.
                        CALL ASSIGN_TO_FC_R
                     ENDIF
                  CASE(JAXIS)
                     IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                        FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                        FLGY = (J >= JLO_FACE) .AND. (J <= JHI_FACE)
                        FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                        INNM = FLGX .AND. FLGY .AND. FLGZ
                        IF (INNM) THEN
                           NOM=NM; IIO=I; JJO=J; KKO=K
                        ELSE
                           CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XCELL(I),YFACE(J),ZCELL(K),NOM,IIO,JJO,KKO)
                        ENDIF
                        IF(NOM==0) EP_TAG(INPE) = .TRUE.
                        CALL ASSIGN_TO_FC_R
                     ENDIF
                  CASE(KAXIS)
                     IF (IJKFACE2(LOW_IND,I,J,K,X1AXIS) < 1 ) THEN
                        FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                        FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                        FLGZ = (K >= KLO_FACE) .AND. (K <= KHI_FACE)
                        INNM = FLGX .AND. FLGY .AND. FLGZ
                        IF (INNM) THEN
                           NOM=NM; IIO=I; JJO=J; KKO=K
                        ELSE
                           CALL SEARCH_OTHER_MESHES_FACE(NM,X1AXIS,XCELL(I),YCELL(J),ZFACE(K),NOM,IIO,JJO,KKO)
                        ENDIF
                        IF(NOM==0) EP_TAG(INPE) = .TRUE.
                        CALL ASSIGN_TO_FC_R
                     ENDIF
                  END SELECT
                  CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE) = IJKFACE2(LOW_IND:HIGH_IND,I,J,K,X1AXIS)
               ENDDO
               ! Now restrict count on cut-face :
               IF(ANY(EP_TAG .EQV. .TRUE.)) CALL RESTRICT_EP(CC_ETYPE_EP)
               DEALLOCATE(EP_TAG)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(IJKFACE2)

   ! Now Cell Variables:
   ALLOCATE(IJKCELL(LOW_IND:HIGH_IND,ISTR:IEND,JSTR:JEND,KSTR:KEND)); IJKCELL = CC_UNDEFINED

   ! 1. Cut-cells:
   VIND = 0
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      I = CUT_CELL(ICC)%IJK(IAXIS)
      J = CUT_CELL(ICC)%IJK(JAXIS)
      K = CUT_CELL(ICC)%IJK(KAXIS)
      ! Don't count cut-cells inside an OBST:
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO ICELL=0,CUT_CELL(ICC)%NCELL
         DO EP=1,INT_N_EXT_PTS  ! External point for face IFACE
            INT_NPE_LO = CUT_CELL(ICC)%INT_NPE(LOW_IND,VIND,EP,ICELL)
            INT_NPE_HI = CUT_CELL(ICC)%INT_NPE(HIGH_IND,VIND,EP,ICELL)
            X1AXIS = VIND
            ALLOCATE(EP_TAG(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)); EP_TAG(:)=.FALSE.
            DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
               I = CUT_CELL(ICC)%INT_IJK(IAXIS,INPE)
               J = CUT_CELL(ICC)%INT_IJK(JAXIS,INPE)
               K = CUT_CELL(ICC)%INT_IJK(KAXIS,INPE)
               ! If cell not counted yet:
               IF (IJKCELL(LOW_IND,I,J,K) < 1 ) THEN
                  FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                  FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                  FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                  INNM = FLGX .AND. FLGY .AND. FLGZ
                  IF (INNM) THEN
                     NOM=NM; IIO=I; JJO=J; KKO=K
                  ELSE
                     CALL SEARCH_OTHER_MESHES(XCELL(I),YCELL(J),ZCELL(K),NOM,IIO,JJO,KKO)
                  ENDIF
                  IF(NOM==0) EP_TAG(INPE) = .TRUE.
                  CALL ASSIGN_TO_CC_R
               ENDIF
               CUT_CELL(ICC)%INT_NOMIND(LOW_IND:HIGH_IND,INPE) = IJKCELL(LOW_IND:HIGH_IND,I,J,K)
            ENDDO
            ! Now restrict count on cut-face :
            IF(ANY(EP_TAG .EQV. .TRUE.)) CALL RESTRICT_EP(CC_FTYPE_CCGAS)
            DEALLOCATE(EP_TAG)
         ENDDO
      ENDDO
   ENDDO

   ! 2. Cell-centered variables for IBEDGES:
   VIND = 0
   DO IEDGE=1,MESHES(NM)%CC_NIBEDGE
      DO ICD_SGN=-2,2
         IF(ICD_SGN==0) CYCLE
         DO EP=1,INT_N_EXT_PTS  ! External point for IEDGE
            INT_NPE_HI = CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN); IF (INT_NPE_HI<1) CYCLE
            INT_NPE_LO = CC_IBEDGE(IEDGE)%INT_NPE(LOW_IND,VIND,EP,ICD_SGN)
            X1AXIS = VIND
            ALLOCATE(EP_TAG(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)); EP_TAG(:)=.FALSE.
            DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
               I = CC_IBEDGE(IEDGE)%INT_IJK(IAXIS,INPE)
               J = CC_IBEDGE(IEDGE)%INT_IJK(JAXIS,INPE)
               K = CC_IBEDGE(IEDGE)%INT_IJK(KAXIS,INPE)
               ! If cell not counted yet:
               IF (IJKCELL(LOW_IND,I,J,K) < 1 ) THEN
                  FLGX = (I >= ILO_CELL) .AND. (I <= IHI_CELL)
                  FLGY = (J >= JLO_CELL) .AND. (J <= JHI_CELL)
                  FLGZ = (K >= KLO_CELL) .AND. (K <= KHI_CELL)
                  INNM = FLGX .AND. FLGY .AND. FLGZ
                  IF (INNM) THEN
                     NOM=NM; IIO=I; JJO=J; KKO=K
                  ELSE
                     CALL SEARCH_OTHER_MESHES(XCELL(I),YCELL(J),ZCELL(K),NOM,IIO,JJO,KKO)
                  ENDIF
                  IF(NOM==0) EP_TAG(INPE) = .TRUE.
                  CALL ASSIGN_TO_CC_R
               ENDIF
               CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE) = IJKCELL(LOW_IND:HIGH_IND,I,J,K)
            ENDDO
            ! Now restrict count on cut-face :
            IF(ANY(EP_TAG .EQV. .TRUE.)) CALL RESTRICT_EP(CC_ETYPE_EP)
            DEALLOCATE(EP_TAG)
         ENDDO
      ENDDO
   ENDDO

   ! Add ghost-cells which are of type CC_CUTCFE or next to one cell type CC_CUTCFE:
   ! First record size of interpolation cells to be reveiced from OMESHES:
   DO NOM=1,NMESHES
      OMESH(NOM)%NCC_INT_R=OMESH(NOM)%NFCC_R(2)
   ENDDO
   ! Now loop INTERPOLATED WALL_CELLs:
   EXT_WALL_LOOP : DO IW=1,N_EXTERNAL_WALL_CELLS

      WC=>WALL(IW)
      EWC=>EXTERNAL_WALL(IW)
      IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXT_WALL_LOOP

      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      II = BC%II
      JJ = BC%JJ
      KK = BC%KK
      NOM = EWC%NOM
      IF (NOM <= 0) CYCLE EXT_WALL_LOOP

      IF(ANY(CCVAR(II-1:II+1,JJ-1:JJ+1,KK-1:KK+1,CC_CGSC)==CC_CUTCFE)) THEN
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                OMESH(NOM)%NFCC_R(2)= OMESH(NOM)%NFCC_R(2) + 1
                SIZE_REC=SIZE(OMESH(NOM)%IIO_CC_R,DIM=1)
                IF(OMESH(NOM)%NFCC_R(2) > SIZE_REC) THEN
                    ALLOCATE(IIO_CC_R_AUX(SIZE_REC),JJO_CC_R_AUX(SIZE_REC),KKO_CC_R_AUX(SIZE_REC));
                    IIO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%IIO_CC_R(1:SIZE_REC)
                    JJO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%JJO_CC_R(1:SIZE_REC)
                    KKO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%KKO_CC_R(1:SIZE_REC)
                    DEALLOCATE(OMESH(NOM)%IIO_CC_R); ALLOCATE(OMESH(NOM)%IIO_CC_R(SIZE_REC+DELTA_FC))
                    OMESH(NOM)%IIO_CC_R(1:SIZE_REC)=IIO_CC_R_AUX(1:SIZE_REC)
                    DEALLOCATE(OMESH(NOM)%JJO_CC_R); ALLOCATE(OMESH(NOM)%JJO_CC_R(SIZE_REC+DELTA_FC))
                    OMESH(NOM)%JJO_CC_R(1:SIZE_REC)=JJO_CC_R_AUX(1:SIZE_REC)
                    DEALLOCATE(OMESH(NOM)%KKO_CC_R); ALLOCATE(OMESH(NOM)%KKO_CC_R(SIZE_REC+DELTA_FC))
                    OMESH(NOM)%KKO_CC_R(1:SIZE_REC)=KKO_CC_R_AUX(1:SIZE_REC)
                    DEALLOCATE(IIO_CC_R_AUX,JJO_CC_R_AUX,KKO_CC_R_AUX)
                ENDIF
                OMESH(NOM)%IIO_CC_R(OMESH(NOM)%NFCC_R(2)) = IIO
                OMESH(NOM)%JJO_CC_R(OMESH(NOM)%NFCC_R(2)) = JJO
                OMESH(NOM)%KKO_CC_R(OMESH(NOM)%NFCC_R(2)) = KKO
               ENDDO
            ENDDO
         ENDDO
       ENDIF
   ENDDO EXT_WALL_LOOP
   DEALLOCATE(IJKCELL)

   ! WRITE(LU_ERR,*) ' MY_RANK,   NM,   NOM,   OMESH(NOM)%NFC_R,  OMESH(NOM)%NCC_R'
   ! DO NOM=1,NMESHES
   !    WRITE(LU_ERR,*) MY_RANK,NM,NOM,OMESH(NOM)%NFC_R,OMESH(NOM)%NCC_R
   ! ENDDO
   ! WRITE(LU_ERR,*) ' '

   ! Quality control:
   ! print*, 'MESHES(NM)%CC_NRCELL_H=',MESHES(NM)%CC_NRCELL_H
   ! IRC=176 ! Last entry for mesh 24x24x24 on sphre_air_demo_1.fds
   ! print*,' '
   ! print*,'RCELL=',IRC
   ! print*,'IJK=',MESHES(NM)%CC_RCELL_H(IRC)%IJK(IAXIS:KAXIS)
   ! print*,'NCCELL=',MESHES(NM)%CC_RCELL_H(IRC)%NCCELL
   ! print*,'CELL_LIST=',MESHES(NM)%CC_RCELL_H(IRC)%CELL_LIST(1:MESHES(NM)%CC_RCELL_H(IRC)%NCCELL)
   ! print*,'INBFC_CARTCEN(1:3)=',MESHES(NM)%CC_RCELL_H(IRC)%INBFC_CARTCEN(1:3)
   ! print*,'INTCOEF_CARTCEN(1:5)=',MESHES(NM)%CC_RCELL_H(IRC)%INTCOEF_CARTCEN(1:5)


   ! Deallocate arrays:
   ! Face centered positions and cell sizes:
   IF (ALLOCATED(XFACE)) DEALLOCATE(XFACE)
   IF (ALLOCATED(YFACE)) DEALLOCATE(YFACE)
   IF (ALLOCATED(ZFACE)) DEALLOCATE(ZFACE)
   IF (ALLOCATED(DXFACE)) DEALLOCATE(DXFACE)
   IF (ALLOCATED(DYFACE)) DEALLOCATE(DYFACE)
   IF (ALLOCATED(DZFACE)) DEALLOCATE(DZFACE)

   ! Cell centered positions and cell sizes:
   IF (ALLOCATED(XCELL)) DEALLOCATE(XCELL)
   IF (ALLOCATED(YCELL)) DEALLOCATE(YCELL)
   IF (ALLOCATED(ZCELL)) DEALLOCATE(ZCELL)
   IF (ALLOCATED(DXCELL)) DEALLOCATE(DXCELL)
   IF (ALLOCATED(DYCELL)) DEALLOCATE(DYCELL)
   IF (ALLOCATED(DZCELL)) DEALLOCATE(DZCELL)

ENDDO MESHES_LOOP

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START_LOOP,' sec.'
   WRITE(LU_SETCC,'(A)') &
   ' - Into FILL_IJKO_INTERP_STENCILS MPI communication..'
   CALL CPU_TIME(CPUTIME_START_LOOP)
ENDIF

! Finally Exchange info on messages to send among MPI processes:
! Populates OMESH(NOM)% : NFCC_S, IIO_FCC_S, JJO_FCC_S, KKO_FCC_S, AXS_FCC_S
CALL FILL_IJKO_INTERP_STENCILS

! Fill unpacking arrays
CALL CC_EXCHANGE_UNPACKING_ARRAYS


IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') '   Done FILL_IJKO_INTERP_STENCILS. Time taken : ',CPUTIME-CPUTIME_START_LOOP,' sec.'
ENDIF

IF (DEBUG_CC_INTERPOLATION) THEN
   ! Write IBSEGS normals:
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      WRITE(MSEGS_FILE,'(A,A,I4.4,A)') TRIM(CHID),'_ibsegns_mesh_',NM,'.dat'
      LU_DB_CCIB = GET_FILE_NUMBER()
      OPEN(LU_DB_CCIB,FILE=TRIM(MSEGS_FILE),STATUS='UNKNOWN')
      DO ECOUNT=1,MESHES(NM)%CC_NIBEDGE
         WRITE(LU_DB_CCIB,'(5F13.8)') 0.,0.,0.,0.,0.
      ENDDO
      CLOSE(LU_DB_CCIB)
   ENDDO
ENDIF

RETURN

CONTAINS

! ----------------------------- COMPUTE_DCOEF ----------------------------------

SUBROUTINE COMPUTE_DCOEF(DATA_IN)

INTEGER, INTENT(IN) :: DATA_IN

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)   :: MASK_IJK
REAL(EB),ALLOCATABLE, DIMENSION(:,:,:)   :: N2
INTEGER  :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER  :: II,JJ,KK,DUMAXIS,COUNT,NEDGI
REAL(EB),ALLOCATABLE, DIMENSION(:,:,:,:) :: RAW_DCOEF
REAL(EB) :: XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS:KAXIS),XYZE(IAXIS:KAXIS),X_P(IAXIS:KAXIS)
LOGICAL :: EVAL=.FALSE.

IF(DATA_IN==CC_ETYPE_RCGAS) THEN
   IF(ALLOCATED(CC_RCEDGE(IEDGE)%INT_DCOEF)) EVAL = .TRUE.
ELSE
   IF(ALLOCATED(CUT_FACE(ICF)%INT_DCOEF)) EVAL = .TRUE.
ENDIF

IF(EVAL) THEN
   IF(DATA_IN==CC_ETYPE_RCGAS) THEN
      NPE_COUNT = CC_RCEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0)
      IF(NPE_COUNT==0) RETURN
      ! Zero DCOEF if any of the box edges vertices is missing:
      ILO = MINVAL(CC_RCEDGE(IEDGE)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      IHI = MAXVAL(CC_RCEDGE(IEDGE)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JLO = MINVAL(CC_RCEDGE(IEDGE)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JHI = MAXVAL(CC_RCEDGE(IEDGE)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KLO = MINVAL(CC_RCEDGE(IEDGE)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KHI = MAXVAL(CC_RCEDGE(IEDGE)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
   ELSEIF(DATA_IN==CC_ETYPE_SCINB) THEN
      NPE_COUNT = CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0)
      IF(NPE_COUNT==0) RETURN
      ! Zero DCOEF if any of the box edges vertices is missing:
      ILO = MINVAL(CC_IBEDGE(IEDGE)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      IHI = MAXVAL(CC_IBEDGE(IEDGE)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JLO = MINVAL(CC_IBEDGE(IEDGE)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JHI = MAXVAL(CC_IBEDGE(IEDGE)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KLO = MINVAL(CC_IBEDGE(IEDGE)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KHI = MAXVAL(CC_IBEDGE(IEDGE)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
   ELSE
      NPE_COUNT = CUT_FACE(ICF)%INT_NPE(HIGH_IND,VIND,EP,IFACE)
      IF(NPE_COUNT==0) RETURN
      ! Zero DCOEF if any of the box edges vertices is missing:
      ILO = MINVAL(CUT_FACE(ICF)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      IHI = MAXVAL(CUT_FACE(ICF)%INT_IJK(IAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JLO = MINVAL(CUT_FACE(ICF)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      JHI = MAXVAL(CUT_FACE(ICF)%INT_IJK(JAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KLO = MINVAL(CUT_FACE(ICF)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
      KHI = MAXVAL(CUT_FACE(ICF)%INT_IJK(KAXIS,INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT))
   ENDIF

   ALLOCATE(INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI))
   INT_DCOEF = 0._EB
   ALLOCATE(MASK_IJK(ILO:IHI,JLO:JHI,KLO:KHI)); MASK_IJK = 0
   ALLOCATE(RAW_DCOEF(IAXIS:KAXIS,ILO:IHI,JLO:JHI,KLO:KHI)); RAW_DCOEF=0._EB;
   ALLOCATE(N2(ILO:IHI,JLO:JHI,KLO:KHI))

   IF(VIND==IAXIS) THEN
      XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS) = (/ XFACE(ILO), XFACE(IHI) /)
   ELSE
      XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS) = (/ XCELL(ILO), XCELL(IHI) /)
   ENDIF
   IF(VIND==JAXIS) THEN
      XYZ_LOHI(LOW_IND:HIGH_IND,JAXIS) = (/ YFACE(JLO), YFACE(JHI) /)
   ELSE
      XYZ_LOHI(LOW_IND:HIGH_IND,JAXIS) = (/ YCELL(JLO), YCELL(JHI) /)
   ENDIF
   IF(VIND==KAXIS) THEN
      XYZ_LOHI(LOW_IND:HIGH_IND,KAXIS) = (/ ZFACE(KLO), ZFACE(KHI) /)
   ELSE
      XYZ_LOHI(LOW_IND:HIGH_IND,KAXIS) = (/ ZCELL(KLO), ZCELL(KHI) /)
   ENDIF

   ! Define external point:
   IF(DATA_IN==CC_ETYPE_RCGAS) THEN
      XYZE(IAXIS:KAXIS) = CC_RCEDGE(IEDGE)%INT_XYZBF(IAXIS:KAXIS,0) + &
      CC_RCEDGE(IEDGE)%INT_XN(EP,0)*CC_RCEDGE(IEDGE)%INT_NOUT(IAXIS:KAXIS,0)
   ELSEIF(DATA_IN==CC_ETYPE_SCINB) THEN
      XYZE(IAXIS:KAXIS) = CC_IBEDGE(IEDGE)%INT_XYZBF(IAXIS:KAXIS,0) + &
      CC_IBEDGE(IEDGE)%INT_XN(EP,0)*CC_IBEDGE(IEDGE)%INT_NOUT(IAXIS:KAXIS,0)
   ELSE
      XYZE(IAXIS:KAXIS) = CUT_FACE(ICF)%INT_XYZBF(IAXIS:KAXIS,IFACE) + &
      CUT_FACE(ICF)%INT_XN(EP,IFACE)*CUT_FACE(ICF)%INT_NOUT(IAXIS:KAXIS,IFACE)
   ENDIF

   ! Masked Trilinear interpolation coefficients:
   X_P(IAXIS:KAXIS) = 0._EB
   DO DUMAXIS=IAXIS,KAXIS
      IF(ABS(XYZ_LOHI(HIGH_IND,DUMAXIS)-XYZ_LOHI(LOW_IND,DUMAXIS)) > TWENTY_EPSILON_EB) &
      X_P(DUMAXIS) = (XYZE(DUMAXIS)-XYZ_LOHI(LOW_IND,DUMAXIS)) / &
                     (XYZ_LOHI(HIGH_IND,DUMAXIS)-XYZ_LOHI(LOW_IND,DUMAXIS))
   ENDDO

   IF(DATA_IN==CC_ETYPE_RCGAS) THEN
      DO INPE=INT_NPE_LO+1,INT_NPE_LO+NPE_COUNT
      MASK_IJK(CC_RCEDGE(IEDGE)%INT_IJK(IAXIS,INPE),CC_RCEDGE(IEDGE)%INT_IJK(JAXIS,INPE),CC_RCEDGE(IEDGE)%INT_IJK(KAXIS,INPE))=1
      ENDDO
   ELSEIF(DATA_IN==CC_ETYPE_SCINB) THEN
      DO INPE=INT_NPE_LO+1,INT_NPE_LO+NPE_COUNT
      MASK_IJK(CC_IBEDGE(IEDGE)%INT_IJK(IAXIS,INPE),CC_IBEDGE(IEDGE)%INT_IJK(JAXIS,INPE),CC_IBEDGE(IEDGE)%INT_IJK(KAXIS,INPE))=1
      ENDDO
   ELSE
      DO INPE=INT_NPE_LO+1,INT_NPE_LO+NPE_COUNT
      MASK_IJK(CUT_FACE(ICF)%INT_IJK(IAXIS,INPE),CUT_FACE(ICF)%INT_IJK(JAXIS,INPE),CUT_FACE(ICF)%INT_IJK(KAXIS,INPE))=1
      ENDDO
   ENDIF

   ! d/dx : First look at which Points are present as both ends of X edges:
   NEDGI = 0
   DO KK = KLO,KHI
      DO JJ = JLO,JHI
         IF(MASK_IJK(ILO,JJ,KK) == 1 .AND. MASK_IJK(IHI,JJ,KK) == 1) THEN ! Both points on edge are populated:
            NEDGI = NEDGI + 1
            RAW_DCOEF(IAXIS,IHI,JJ,KK) = 1._EB/DXCELL(ILO)
            RAW_DCOEF(IAXIS,ILO,JJ,KK) =-1._EB/DXCELL(ILO)
         ENDIF
      ENDDO
   ENDDO
   ! Regarding the number of Edges interpolate:
   IF (NEDGI > 0 .AND. NEDGI < 4) THEN ! Simple average:
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(IAXIS,II,JJ,KK) = RAW_DCOEF(IAXIS,II,JJ,KK)/REAL(NEDGI,EB)
            ENDDO
         ENDDO
      ENDDO
   ELSE ! Bilinear in y, z directions:
      N2(ILO:IHI,JLO,KLO) = (1._EB-X_P(JAXIS))*(1._EB-X_P(KAXIS))
      N2(ILO:IHI,JHI,KLO) = (      X_P(JAXIS))*(1._EB-X_P(KAXIS))
      N2(ILO:IHI,JLO,KHI) = (1._EB-X_P(JAXIS))*(      X_P(KAXIS))
      N2(ILO:IHI,JHI,KHI) = (      X_P(JAXIS))*(      X_P(KAXIS))
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(IAXIS,II,JJ,KK) = N2(II,JJ,KK)*RAW_DCOEF(IAXIS,II,JJ,KK)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   ! d/dy : First look at which Points are present as both ends of Y edges:
   NEDGI = 0
   DO KK = KLO,KHI
      DO II = ILO,IHI
         IF(MASK_IJK(II,JLO,KK) == 1 .AND. MASK_IJK(II,JHI,KK) == 1) THEN ! Both points on edge are populated:
            NEDGI = NEDGI + 1
            RAW_DCOEF(JAXIS,II,JHI,KK) = 1._EB/DYCELL(JLO)
            RAW_DCOEF(JAXIS,II,JLO,KK) =-1._EB/DYCELL(JLO)
         ENDIF
      ENDDO
   ENDDO
   ! Regarding the number of Edges interpolate:
   IF (NEDGI > 0 .AND. NEDGI < 4) THEN ! Simple average:
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(JAXIS,II,JJ,KK) = RAW_DCOEF(JAXIS,II,JJ,KK)/REAL(NEDGI,EB)
            ENDDO
         ENDDO
      ENDDO
   ELSE ! Bilinear in x, z directions:
      N2(ILO,JLO:JHI,KLO) = (1._EB-X_P(IAXIS))*(1._EB-X_P(KAXIS))
      N2(IHI,JLO:JHI,KLO) = (      X_P(IAXIS))*(1._EB-X_P(KAXIS))
      N2(ILO,JLO:JHI,KHI) = (1._EB-X_P(IAXIS))*(      X_P(KAXIS))
      N2(IHI,JLO:JHI,KHI) = (      X_P(IAXIS))*(      X_P(KAXIS))
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(JAXIS,II,JJ,KK) = N2(II,JJ,KK)*RAW_DCOEF(JAXIS,II,JJ,KK)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   ! d/dz : First look at which Points are present as both ends of Z edges:
   NEDGI = 0
   DO JJ = JLO,JHI
      DO II = ILO,IHI
         IF(MASK_IJK(II,JJ,KLO) == 1 .AND. MASK_IJK(II,JJ,KHI) == 1) THEN ! Both points on edge are populated:
            NEDGI = NEDGI + 1
            RAW_DCOEF(KAXIS,II,JJ,KHI) = 1._EB/DZCELL(KLO)
            RAW_DCOEF(KAXIS,II,JJ,KLO) =-1._EB/DZCELL(KLO)
         ENDIF
      ENDDO
   ENDDO
   ! Regarding the number of Edges interpolate:
   IF (NEDGI > 0 .AND. NEDGI < 4) THEN ! Simple average:
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(KAXIS,II,JJ,KK) = RAW_DCOEF(KAXIS,II,JJ,KK)/REAL(NEDGI,EB)
            ENDDO
         ENDDO
      ENDDO
   ELSE ! Bilinear in x, y directions:
      N2(ILO,JLO,KLO:KHI) = (1._EB-X_P(IAXIS))*(1._EB-X_P(JAXIS))
      N2(IHI,JLO,KLO:KHI) = (      X_P(IAXIS))*(1._EB-X_P(JAXIS))
      N2(ILO,JHI,KLO:KHI) = (1._EB-X_P(IAXIS))*(      X_P(JAXIS))
      N2(IHI,JHI,KLO:KHI) = (      X_P(IAXIS))*(      X_P(JAXIS))
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_DCOEF(KAXIS,II,JJ,KK) = N2(II,JJ,KK)*RAW_DCOEF(KAXIS,II,JJ,KK)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
   ! Finally populate INT_DCOEF:
   COUNT = 0
   DO KK = KLO,KHI
      DO JJ = JLO,JHI
         DO II = ILO,IHI
            IF(MASK_IJK(II,JJ,KK) == 1) THEN
               COUNT = COUNT + 1
               INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+COUNT) = RAW_DCOEF(IAXIS:KAXIS,II,JJ,KK)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   IF(DATA_IN==CC_ETYPE_RCGAS) THEN
      CC_RCEDGE(IEDGE)%INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
      INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   ELSEIF(DATA_IN==CC_ETYPE_SCINB) THEN
      CC_IBEDGE(IEDGE)%INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
      INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   ELSE
      CUT_FACE(ICF)%INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
      INT_DCOEF(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   ENDIF
   DEALLOCATE(INT_DCOEF,MASK_IJK,N2,RAW_DCOEF)
ENDIF

RETURN
END SUBROUTINE COMPUTE_DCOEF

! ------------------------------ RESTRICT_EP ----------------------------------

SUBROUTINE RESTRICT_EP(CFRC_FLG)

INTEGER, INTENT(IN) :: CFRC_FLG

REAL(EB):: PROD_COEF

! Apply restriction to stencil points with NOM>0:
ALLOCATE(INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI),   &
         INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI),              &
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI))
INT_IJK = CC_UNDEFINED; INT_COEF = 0._EB; INT_NOMIND = CC_UNDEFINED
NPE_COUNT = 0
PROD_COEF= 0._EB
SELECT CASE(CFRC_FLG)
CASE(CC_FTYPE_CFGAS,CC_FTYPE_CFINB)
   DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
      IF(.NOT.EP_TAG(INPE)) THEN ! Point has a NOM /= 0
         NPE_COUNT = NPE_COUNT + 1
         INT_IJK(IAXIS:KAXIS,INT_NPE_LO+NPE_COUNT) = CUT_FACE(ICF)%INT_IJK(IAXIS:KAXIS,INPE)
         INT_COEF(INT_NPE_LO+NPE_COUNT) = CUT_FACE(ICF)%INT_COEF(INPE)
         PROD_COEF = PROD_COEF + INT_COEF(INT_NPE_LO+NPE_COUNT)
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+NPE_COUNT) = &
         CUT_FACE(ICF)%INT_NOMIND(LOW_IND:HIGH_IND,INPE)
      ENDIF
   ENDDO
   IF (ABS(PROD_COEF) < TWENTY_EPSILON_EB) THEN ! Any viable points throught EP_TAG have been discarded by IJKFACE2.
      INT_IJK=CC_UNDEFINED; INT_COEF=0._EB; INT_NOMIND=CC_UNDEFINED; NPE_COUNT=0
   ENDIF
   CUT_FACE(ICF)%INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CUT_FACE(ICF)%INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CUT_FACE(ICF)%INT_NPE(HIGH_IND,VIND,EP,IFACE) = NPE_COUNT

   IF (NPE_COUNT > 0) &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT) = INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT)/PROD_COEF

   CUT_FACE(ICF)%INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

CASE(CC_ETYPE_RCGAS)
   DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
      IF(.NOT.EP_TAG(INPE)) THEN ! Point has a NOM /= 0
         NPE_COUNT = NPE_COUNT + 1
         INT_IJK(IAXIS:KAXIS,INT_NPE_LO+NPE_COUNT) = CC_RCEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INPE)
         INT_COEF(INT_NPE_LO+NPE_COUNT) = CC_RCEDGE(IEDGE)%INT_COEF(INPE)
         PROD_COEF = PROD_COEF + INT_COEF(INT_NPE_LO+NPE_COUNT)
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+NPE_COUNT) = &
         CC_RCEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE)
      ENDIF
   ENDDO
   IF (ABS(PROD_COEF) < TWENTY_EPSILON_EB) THEN ! Any viable points throught EP_TAG have been discarded by IJKFACE2.
      INT_IJK=CC_UNDEFINED; INT_COEF=0._EB; INT_NOMIND=CC_UNDEFINED; NPE_COUNT=0
   ENDIF
   CC_RCEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_RCEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_RCEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0) = NPE_COUNT

   IF (NPE_COUNT > 0) &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT) = INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT)/PROD_COEF

   CC_RCEDGE(IEDGE)%INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

CASE(CC_ETYPE_SCINB)
   DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
      IF(.NOT.EP_TAG(INPE)) THEN ! Point has a NOM /= 0
         NPE_COUNT = NPE_COUNT + 1
         INT_IJK(IAXIS:KAXIS,INT_NPE_LO+NPE_COUNT) = CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INPE)
         INT_COEF(INT_NPE_LO+NPE_COUNT) = CC_IBEDGE(IEDGE)%INT_COEF(INPE)
         PROD_COEF = PROD_COEF + INT_COEF(INT_NPE_LO+NPE_COUNT)
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+NPE_COUNT) = &
         CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE)
      ENDIF
   ENDDO
   IF (ABS(PROD_COEF) < TWENTY_EPSILON_EB) THEN ! Any viable points throught EP_TAG have been discarded by IJKFACE2.
      INT_IJK=CC_UNDEFINED; INT_COEF=0._EB; INT_NOMIND=CC_UNDEFINED; NPE_COUNT=0
   ENDIF
   CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,0) = NPE_COUNT

   IF (NPE_COUNT > 0) &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT) = INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT)/PROD_COEF

   CC_IBEDGE(IEDGE)%INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

CASE(CC_ETYPE_EP)
   DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
      IF(.NOT.EP_TAG(INPE)) THEN ! Point has a NOM /= 0
         NPE_COUNT = NPE_COUNT + 1
         INT_IJK(IAXIS:KAXIS,INT_NPE_LO+NPE_COUNT) = CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INPE)
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+NPE_COUNT) = &
         CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INPE)
         INT_COEF(INT_NPE_LO+NPE_COUNT) = CC_IBEDGE(IEDGE)%INT_DCOEF(INPE,1)
      ENDIF
   ENDDO
   CC_IBEDGE(IEDGE)%INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_IBEDGE(IEDGE)%INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_IBEDGE(IEDGE)%INT_DCOEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI,1) = INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CC_IBEDGE(IEDGE)%INT_NPE(HIGH_IND,VIND,EP,ICD_SGN) = NPE_COUNT

CASE(CC_FTYPE_RCGAS) ! Skip.
CASE(CC_FTYPE_CCGAS)
   DO INPE=INT_NPE_LO+1,INT_NPE_LO+INT_NPE_HI
      IF(.NOT.EP_TAG(INPE)) THEN ! Point has a NOM /= 0
         NPE_COUNT = NPE_COUNT + 1
         INT_IJK(IAXIS:KAXIS,INT_NPE_LO+NPE_COUNT) = CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,INPE)
         INT_COEF(INT_NPE_LO+NPE_COUNT) = CUT_CELL(ICC)%INT_COEF(INPE)
         PROD_COEF = PROD_COEF + INT_COEF(INT_NPE_LO+NPE_COUNT)
         INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+NPE_COUNT) = &
         CUT_CELL(ICC)%INT_NOMIND(LOW_IND:HIGH_IND,INPE)
      ENDIF
   ENDDO
   IF (ABS(PROD_COEF) < TWENTY_EPSILON_EB) THEN ! Any viable points throught EP_TAG have been discarded by IJKCELL.
      INT_IJK=CC_UNDEFINED; INT_COEF=0._EB; INT_NOMIND=CC_UNDEFINED; NPE_COUNT=0
   ENDIF
   CUT_CELL(ICC)%INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_IJK(IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CUT_CELL(ICC)%INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_NOMIND(LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   CUT_CELL(ICC)%INT_NPE(HIGH_IND,VIND,EP,ICELL) = NPE_COUNT

   IF (NPE_COUNT > 0) &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT) = INT_COEF(INT_NPE_LO+1:INT_NPE_LO+NPE_COUNT)/PROD_COEF

   CUT_CELL(ICC)%INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI) = &
   INT_COEF(INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

END SELECT

DEALLOCATE(INT_IJK,INT_NOMIND,INT_COEF)

RETURN
END SUBROUTINE RESTRICT_EP

! ------------------------------- GET_DELN ------------------------------------

SUBROUTINE GET_DELN(FCTN_IN,DELN,DXLOC,DYLOC,DZLOC,NVEC,CLOSE_PT)

REAL(EB), INTENT(IN) :: FCTN_IN,DXLOC,DYLOC,DZLOC
REAL(EB), OPTIONAL, INTENT(IN) :: NVEC(MAX_DIM)
LOGICAL, OPTIONAL, INTENT(IN) :: CLOSE_PT
REAL(EB), INTENT(OUT) :: DELN

! Local Variables:
REAL(EB) :: FCTN
FCTN = FCTN_IN
IF (PRESENT(NVEC)) THEN
   IF( .NOT.PRESENT(CLOSE_PT) .AND. (ABS(NVEC(IAXIS))>GEOMEPS) .AND. &
   (ABS(NVEC(JAXIS))>GEOMEPS) .AND. (ABS(NVEC(KAXIS))>GEOMEPS)) FCTN=SQRT(3._EB)
   ! IF(PRESENT(CLOSE_PT)) THEN
   !    IF(CLOSE_PT) FCTN = 1._EB
   ! ENDIF
   DELN = FCTN*(DXLOC*ABS(NVEC(IAXIS))+DYLOC*ABS(NVEC(JAXIS))+DZLOC*ABS(NVEC(KAXIS)))
ELSE
   DELN = SQRT(DXLOC**2._EB+DYLOC**2._EB+DZLOC**2._EB)
ENDIF
RETURN
END SUBROUTINE GET_DELN

! ---------------------------- ASSIGN_TO_CC_R ---------------------------------

SUBROUTINE ASSIGN_TO_CC_R

 IF (NOM > 0) THEN ! Add to IIO_FC_R,JJO_FC_R,KKO_FC_R,AXIS_FC_R list,
                   ! and add 1 to NFC_R for OMESH(NOM).
    ! Use Automatic reallocation:
    OMESH(NOM)%NFCC_R(2)= OMESH(NOM)%NFCC_R(2) + 1
    SIZE_REC=SIZE(OMESH(NOM)%IIO_CC_R,DIM=1)
    IF(OMESH(NOM)%NFCC_R(2) > SIZE_REC) THEN
        ALLOCATE(IIO_CC_R_AUX(SIZE_REC),JJO_CC_R_AUX(SIZE_REC),KKO_CC_R_AUX(SIZE_REC));
        IIO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%IIO_CC_R(1:SIZE_REC)
        JJO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%JJO_CC_R(1:SIZE_REC)
        KKO_CC_R_AUX(1:SIZE_REC)=OMESH(NOM)%KKO_CC_R(1:SIZE_REC)
        DEALLOCATE(OMESH(NOM)%IIO_CC_R); ALLOCATE(OMESH(NOM)%IIO_CC_R(SIZE_REC+DELTA_FC))
        OMESH(NOM)%IIO_CC_R(1:SIZE_REC)=IIO_CC_R_AUX(1:SIZE_REC)
        DEALLOCATE(OMESH(NOM)%JJO_CC_R); ALLOCATE(OMESH(NOM)%JJO_CC_R(SIZE_REC+DELTA_FC))
        OMESH(NOM)%JJO_CC_R(1:SIZE_REC)=JJO_CC_R_AUX(1:SIZE_REC)
        DEALLOCATE(OMESH(NOM)%KKO_CC_R); ALLOCATE(OMESH(NOM)%KKO_CC_R(SIZE_REC+DELTA_FC))
        OMESH(NOM)%KKO_CC_R(1:SIZE_REC)=KKO_CC_R_AUX(1:SIZE_REC)
        DEALLOCATE(IIO_CC_R_AUX,JJO_CC_R_AUX,KKO_CC_R_AUX)
    ENDIF
    OMESH(NOM)%IIO_CC_R(OMESH(NOM)%NFCC_R(2)) = IIO
    OMESH(NOM)%JJO_CC_R(OMESH(NOM)%NFCC_R(2)) = JJO
    OMESH(NOM)%KKO_CC_R(OMESH(NOM)%NFCC_R(2)) = KKO
    IJKCELL(LOW_IND:HIGH_IND,I,J,K) = (/ NOM, OMESH(NOM)%NFCC_R(2) /)
 ENDIF

RETURN
END SUBROUTINE ASSIGN_TO_CC_R

! ---------------------------- ASSIGN_TO_FC_R ---------------------------------

SUBROUTINE ASSIGN_TO_FC_R

IF (NOM > 0) THEN
   ! Use Automatic reallocation:
   OMESH(NOM)%NFCC_R(1)= OMESH(NOM)%NFCC_R(1) + 1
   SIZE_REC=SIZE(OMESH(NOM)%IIO_FC_R,DIM=1)
   IF(OMESH(NOM)%NFCC_R(1) > SIZE_REC) THEN
       ALLOCATE(IIO_FC_R_AUX(SIZE_REC),JJO_FC_R_AUX(SIZE_REC),KKO_FC_R_AUX(SIZE_REC));
       ALLOCATE(AXS_FC_R_AUX(SIZE_REC))
       IIO_FC_R_AUX(1:SIZE_REC)=OMESH(NOM)%IIO_FC_R(1:SIZE_REC)
       JJO_FC_R_AUX(1:SIZE_REC)=OMESH(NOM)%JJO_FC_R(1:SIZE_REC)
       KKO_FC_R_AUX(1:SIZE_REC)=OMESH(NOM)%KKO_FC_R(1:SIZE_REC)
       AXS_FC_R_AUX(1:SIZE_REC)=OMESH(NOM)%AXS_FC_R(1:SIZE_REC)
       DEALLOCATE(OMESH(NOM)%IIO_FC_R); ALLOCATE(OMESH(NOM)%IIO_FC_R(SIZE_REC+DELTA_FC))
       OMESH(NOM)%IIO_FC_R(1:SIZE_REC)=IIO_FC_R_AUX(1:SIZE_REC)
       DEALLOCATE(OMESH(NOM)%JJO_FC_R); ALLOCATE(OMESH(NOM)%JJO_FC_R(SIZE_REC+DELTA_FC))
       OMESH(NOM)%JJO_FC_R(1:SIZE_REC)=JJO_FC_R_AUX(1:SIZE_REC)
       DEALLOCATE(OMESH(NOM)%KKO_FC_R); ALLOCATE(OMESH(NOM)%KKO_FC_R(SIZE_REC+DELTA_FC))
       OMESH(NOM)%KKO_FC_R(1:SIZE_REC)=KKO_FC_R_AUX(1:SIZE_REC)
       DEALLOCATE(OMESH(NOM)%AXS_FC_R); ALLOCATE(OMESH(NOM)%AXS_FC_R(SIZE_REC+DELTA_FC))
       OMESH(NOM)%AXS_FC_R(1:SIZE_REC)=AXS_FC_R_AUX(1:SIZE_REC)
       DEALLOCATE(IIO_FC_R_AUX,JJO_FC_R_AUX,KKO_FC_R_AUX,AXS_FC_R_AUX)
   ENDIF
   OMESH(NOM)%IIO_FC_R(OMESH(NOM)%NFCC_R(1)) = IIO
   OMESH(NOM)%JJO_FC_R(OMESH(NOM)%NFCC_R(1)) = JJO
   OMESH(NOM)%KKO_FC_R(OMESH(NOM)%NFCC_R(1)) = KKO
   OMESH(NOM)%AXS_FC_R(OMESH(NOM)%NFCC_R(1)) = X1AXIS
   IJKFACE2(LOW_IND:HIGH_IND,I,J,K,X1AXIS) = (/ NOM, OMESH(NOM)%NFCC_R(1) /)
ENDIF

RETURN
END SUBROUTINE ASSIGN_TO_FC_R

! ---------------------------- GET_INTSTENCILS_EP -------------------------------

SUBROUTINE GET_INTSTENCILS_EP(MASK_FLG,VIND,XYZ_PP,INTXN,NVEC,NPE_LIST_START,NPE_LIST_COUNT,&
INT_IJK,INT_COEF)

! This routine provides a set of interpolation points for an external normal point EP,
! located at position XYZE(IAXIS:KAXIS) = XYZ_PP(IAXIS:KAXIS) + INTXN*NVEC(IAXIS:KAXIS):
! The points will be face centered when:
!     VIND = IAXIS => X faces
!            JAXIS => Y faces
!            KAXIS => Z faces
! And cell centered when VIND = 0.
! The number of interpolation points is provided in variable NPE_LIST_COUNT.
! The IJK indexes on mesh of these points is defined in:
! INT_IJK(IAXIS:KAXIS,NPE_LIST_START+1:NPE_LIST_START+NPE_LIST_COUNT)
! INT_COEF(NPE_LIST_START+1:NPE_LIST_START+NPE_LIST_COUNT)

LOGICAL, INTENT(IN) :: MASK_FLG
INTEGER, INTENT(IN) :: VIND,NPE_LIST_START
INTEGER, INTENT(OUT):: NPE_LIST_COUNT
REAL(EB),INTENT(IN) :: XYZ_PP(IAXIS:KAXIS), NVEC(IAXIS:KAXIS), INTXN
INTEGER, INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: INT_IJK
REAL(EB), INTENT(INOUT),ALLOCATABLE, DIMENSION(:)   :: INT_COEF


! Local variables:
REAL(EB) :: XYZE(IAXIS:KAXIS)
LOGICAL  :: IS_FACE_X,IS_FACE_Y,IS_FACE_Z
INTEGER  :: INDX,INDY,INDZ,DIM_NPE,ILO,IHI,JLO,JHI,KLO,KHI,ILO_2,IHI_2,JLO_2,JHI_2,KLO_2,KHI_2
INTEGER  :: II,JJ,KK,DUMAXIS,COUNT
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MASK_IJK
REAL(EB),ALLOCATABLE, DIMENSION(:,:,:) :: RAW_COEF
REAL(EB) :: XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS:KAXIS),X_P(IAXIS:KAXIS),RED_COEF

LOGICAL, PARAMETER :: DO_TRILINEAR = .TRUE.

! Default number of interpolation stencil points:
NPE_LIST_COUNT = 0

! Define external point:
XYZE(IAXIS:KAXIS) = XYZ_PP(IAXIS:KAXIS) + INTXN*NVEC(IAXIS:KAXIS)

! Find closest point on mesh NM to point:
IS_FACE_X=.FALSE.;IS_FACE_Y=.FALSE.;IS_FACE_Z=.FALSE.
IF(VIND==IAXIS) IS_FACE_X=.TRUE.
IF(VIND==JAXIS) IS_FACE_Y=.TRUE.
IF(VIND==KAXIS) IS_FACE_Z=.TRUE.
CALL GET_X_IND(XYZE,IS_FACE_X,INDX)
CALL GET_Y_IND(XYZE,IS_FACE_Y,INDY)
CALL GET_Z_IND(XYZE,IS_FACE_Z,INDZ)

IF(INDX == INDEX_UNDEFINED) RETURN
IF(INDY == INDEX_UNDEFINED) RETURN
IF(INDZ == INDEX_UNDEFINED) RETURN

! Define stencil points:
DIM_NPE = SIZE(INT_IJK, DIM=2)
IF(NPE_LIST_START+MAX_INTERP_POINTS > DIM_NPE) THEN ! Reallocate size of INT_IJK, INT_COEF
   ALLOCATE(INT_IJK_AUX(IAXIS:KAXIS,DIM_NPE),INT_COEF_AUX(1:DIM_NPE))
   INT_IJK_AUX(IAXIS:KAXIS,1:DIM_NPE) = INT_IJK(IAXIS:KAXIS,1:DIM_NPE)
   INT_COEF_AUX(1:DIM_NPE)            = INT_COEF(1:DIM_NPE)
   DEALLOCATE(INT_IJK, INT_COEF)
   ALLOCATE(INT_IJK(IAXIS:KAXIS,NPE_LIST_START+MAX_INTERP_POINTS+DELTA_VERT)); INT_IJK = CC_UNDEFINED
   ALLOCATE(INT_COEF(1:NPE_LIST_START+MAX_INTERP_POINTS+DELTA_VERT)); INT_COEF = 0._EB
   INT_IJK(IAXIS:KAXIS,1:DIM_NPE) = INT_IJK_AUX(IAXIS:KAXIS,1:DIM_NPE)
   INT_COEF(1:DIM_NPE)            = INT_COEF_AUX(1:DIM_NPE)
   DEALLOCATE(INT_IJK_AUX,INT_COEF_AUX)
ENDIF

! Linear interpolation bounds:
ILO = INDX-1;  IHI = INDX
JLO = INDY-1;  JHI = INDY
KLO = INDZ-1;  KHI = INDZ
! Other interpolation bounds:
IF (STENCIL_INTERPOLATION /= CC_LINEAR_INTERPOLATION) THEN ! Either QUADRATIC_INTERPOLATION,WLS_INTERPOLATION.
   ILO_2 = ILO_CELL; IHI_2 = IHI_CELL
   JLO_2 = JLO_CELL; JHI_2 = JHI_CELL
   KLO_2 = KLO_CELL; KHI_2 = KHI_CELL
   SELECT CASE(VIND)
   CASE(IAXIS)
      ILO_2 = ILO_FACE; IHI_2 = IHI_FACE
   CASE(JAXIS)
      JLO_2 = JLO_FACE; JHI_2 = JHI_FACE
   CASE(KAXIS)
      KLO_2 = KLO_FACE; KHI_2 = KHI_FACE
   END SELECT
   IF(IHI == IHI_2+NGUARD ) THEN
      ILO = ILO - 1
   ELSEIF(ILO >= ILO_2-NGUARD ) THEN
      IHI = IHI + 1
   ENDIF
   IF(JHI == JHI_2+NGUARD ) THEN
      JLO = JLO - 1
   ELSEIF(JLO >= JLO_2-NGUARD ) THEN
      JHI = JHI + 1
   ENDIF
   IF(KHI == KHI_2+NGUARD ) THEN
      KLO = KLO - 1
   ELSEIF(KLO >= KLO_2-NGUARD ) THEN
      KHI = KHI + 1
   ENDIF

ENDIF

! Allocate stencil Allocatable arrays:
ALLOCATE(MASK_IJK(ILO:IHI,JLO:JHI,KLO:KHI)); MASK_IJK=0;
ALLOCATE(RAW_COEF(ILO:IHI,JLO:JHI,KLO:KHI)); RAW_COEF=0._EB;

! Add collocation points to interpolation stencil:
! Face vars:
IF(VIND > 0) THEN
   IF (MASK_FLG) THEN
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               IF(FCVAR(II,JJ,KK,CC_IDRC,VIND) /= CC_GASPHASE) CYCLE ! Cycle if facevar is masked by CC_IDRC.
               NPE_LIST_COUNT = NPE_LIST_COUNT + 1
               INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/ II, JJ, KK /)
               ! Coeff computation left for the end.
               MASK_IJK(II,JJ,KK) = 1
            ENDDO
         ENDDO
      ENDDO
   ELSE
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               IF(FCVAR(II,JJ,KK,CC_FGSC,VIND) == CC_SOLID) CYCLE ! Cycle solid faces.
               NPE_LIST_COUNT = NPE_LIST_COUNT + 1
               INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/ II, JJ, KK /)
               ! Coeff computation left for the end.
               MASK_IJK(II,JJ,KK) = 1
            ENDDO
         ENDDO
      ENDDO
   ENDIF ! MASK_FLG
ELSE ! Centered vars:
   DO KK = KLO,KHI
      DO JJ = JLO,JHI
         DO II = ILO,IHI
            IF(CCVAR(II,JJ,KK,CC_CGSC) == CC_SOLID) CYCLE ! Cycle solid cells.
            NPE_LIST_COUNT = NPE_LIST_COUNT + 1
            INT_IJK(IAXIS:KAXIS,NPE_LIST_START+NPE_LIST_COUNT) = (/ II, JJ, KK /)
            ! Coeff computation left for the end.
            MASK_IJK(II,JJ,KK) = 1
         ENDDO
      ENDDO
   ENDDO
ENDIF

! If NPE_LIST_COUNT == 0 return. Will use boundary value only on interpolated face:
IF (NPE_LIST_COUNT == 0) THEN
   DEALLOCATE(MASK_IJK,RAW_COEF)
   RETURN
ENDIF


! At this point we have the interpolation stencil points for EP and VIND mesh.
! Regarding the interpolation type chosen produce the interpolation coefficients in INT_COEFF.
! Define Bounding Box of the Stencil:
IF(VIND==IAXIS) THEN
   XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS) = (/ XFACE(ILO), XFACE(IHI) /)
ELSE
   XYZ_LOHI(LOW_IND:HIGH_IND,IAXIS) = (/ XCELL(ILO), XCELL(IHI) /)
ENDIF
IF(VIND==JAXIS) THEN
   XYZ_LOHI(LOW_IND:HIGH_IND,JAXIS) = (/ YFACE(JLO), YFACE(JHI) /)
ELSE
   XYZ_LOHI(LOW_IND:HIGH_IND,JAXIS) = (/ YCELL(JLO), YCELL(JHI) /)
ENDIF
IF(VIND==KAXIS) THEN
   XYZ_LOHI(LOW_IND:HIGH_IND,KAXIS) = (/ ZFACE(KLO), ZFACE(KHI) /)
ELSE
   XYZ_LOHI(LOW_IND:HIGH_IND,KAXIS) = (/ ZCELL(KLO), ZCELL(KHI) /)
ENDIF

! Masked Trilinear interpolation coefficients:
IF (STENCIL_INTERPOLATION == CC_LINEAR_INTERPOLATION) THEN
   DO DUMAXIS=IAXIS,KAXIS
      X_P(DUMAXIS) = (XYZE(DUMAXIS)-XYZ_LOHI(LOW_IND,DUMAXIS))/(XYZ_LOHI(HIGH_IND,DUMAXIS)-XYZ_LOHI(LOW_IND,DUMAXIS))
   ENDDO

   ! Case of Trilinear interpolation:
   DO_TRILINEAR_COND : IF (DO_TRILINEAR) THEN
      ! Masked trilinear interpolation:
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               RAW_COEF(II,JJ,KK) = (REAL(II-ILO,EB)*X_P(IAXIS)+REAL(IHI-II,EB)*(1._EB-X_P(IAXIS))) * &
                                    (REAL(JJ-JLO,EB)*X_P(JAXIS)+REAL(JHI-JJ,EB)*(1._EB-X_P(JAXIS))) * &
                                    (REAL(KK-KLO,EB)*X_P(KAXIS)+REAL(KHI-KK,EB)*(1._EB-X_P(KAXIS)))
            ENDDO
         ENDDO
      ENDDO
      ! Rescale remaining coefficients:
      RED_COEF = 0._EB
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               IF (MASK_IJK(II,JJ,KK) == 1) RED_COEF = RED_COEF + RAW_COEF(II,JJ,KK)
            ENDDO
         ENDDO
      ENDDO
      IF (ABS(RED_COEF) < TWENTY_EPSILON_EB) THEN
         NPE_LIST_COUNT = 0
         DEALLOCATE(MASK_IJK,RAW_COEF)
         RETURN
      ENDIF
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               IF (MASK_IJK(II,JJ,KK) == 1) THEN
                  RAW_COEF(II,JJ,KK) = RAW_COEF(II,JJ,KK)/RED_COEF
               ELSE
                  RAW_COEF(II,JJ,KK) = 0._EB
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ! Finally add coefficients to INT_COEF:
      COUNT = 0
      DO KK = KLO,KHI
         DO JJ = JLO,JHI
            DO II = ILO,IHI
               IF(MASK_IJK(II,JJ,KK) == 1) THEN
                  COUNT = COUNT + 1
                  INT_COEF(NPE_LIST_START+COUNT) = RAW_COEF(II,JJ,KK)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

   ! Case of Least Squares interpolation with up to 8 stencil points:
   ELSE
      ! To do.

   ENDIF DO_TRILINEAR_COND
   DEALLOCATE(MASK_IJK,RAW_COEF)
   RETURN
ENDIF

! Other interpolation schemes:
! To do.

DEALLOCATE(MASK_IJK,RAW_COEF)
RETURN
END SUBROUTINE GET_INTSTENCILS_EP

SUBROUTINE GET_X_IND(XYZE,IS_FACE,INDX)
REAL(EB),INTENT(IN) :: XYZE(IAXIS:KAXIS)
LOGICAL, INTENT(IN) :: IS_FACE
INTEGER, INTENT(OUT):: INDX
INTEGER :: IEP
INDX = INDEX_UNDEFINED
IF (IS_FACE) THEN ! X face.
   IF(XYZE(IAXIS) >= XFACE(ILO_FACE-NGUARD)) THEN
      DO IEP=ILO_FACE-CCGUARD,IHI_FACE+CCGUARD
         IF (XYZE(IAXIS)+GEOFCT*GEOMEPS < XFACE(IEP)) THEN ! First X index that XYZ(IAXIS) is smaller.
            INDX = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ELSE ! X center.
   IF(XYZE(IAXIS) >= XCELL(ILO_CELL-NGUARD)) THEN
      DO IEP=ILO_CELL-CCGUARD,IHI_CELL+CCGUARD
         IF (XYZE(IAXIS)+GEOFCT*GEOMEPS < XCELL(IEP)) THEN ! First X index that XYZ(IAXIS) is smaller.
            INDX = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ENDIF
END SUBROUTINE GET_X_IND

SUBROUTINE GET_Y_IND(XYZE,IS_FACE,INDY)
REAL(EB),INTENT(IN) :: XYZE(IAXIS:KAXIS)
LOGICAL, INTENT(IN) :: IS_FACE
INTEGER, INTENT(OUT):: INDY
INTEGER :: IEP
INDY = INDEX_UNDEFINED
IF (IS_FACE) THEN ! Y face.
   IF(XYZE(JAXIS) >= YFACE(JLO_FACE-NGUARD)) THEN
      DO IEP=JLO_FACE-CCGUARD,JHI_FACE+CCGUARD
         IF (XYZE(JAXIS)+GEOFCT*GEOMEPS < YFACE(IEP)) THEN ! First Y index that XYZ(JAXIS) is smaller.
            INDY = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ELSE ! Y center.
   IF(XYZE(JAXIS) >= YCELL(JLO_CELL-NGUARD)) THEN
      DO IEP=JLO_CELL-CCGUARD,JHI_CELL+CCGUARD
         IF (XYZE(JAXIS)+GEOFCT*GEOMEPS < YCELL(IEP)) THEN ! First Y index that XYZ(JAXIS) is smaller.
            INDY = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ENDIF
END SUBROUTINE GET_Y_IND

SUBROUTINE GET_Z_IND(XYZE,IS_FACE,INDZ)
REAL(EB),INTENT(IN) :: XYZE(IAXIS:KAXIS)
LOGICAL, INTENT(IN) :: IS_FACE
INTEGER, INTENT(OUT):: INDZ
INTEGER :: IEP
INDZ = INDEX_UNDEFINED
IF (IS_FACE) THEN ! Z face.
   IF(XYZE(KAXIS) >= ZFACE(KLO_FACE-NGUARD)) THEN
      DO IEP=KLO_FACE-CCGUARD,KHI_FACE+CCGUARD
         IF (XYZE(KAXIS)+GEOFCT*GEOMEPS < ZFACE(IEP)) THEN ! First Z index that XYZ(KAXIS) is smaller.
            INDZ = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ELSE ! Z center.
   IF(XYZE(KAXIS) >= ZCELL(KLO_CELL-NGUARD)) THEN
      DO IEP=KLO_CELL-CCGUARD,KHI_CELL+CCGUARD
         IF (XYZE(KAXIS)+GEOFCT*GEOMEPS < ZCELL(IEP)) THEN ! First Z index that XYZ(KAXIS) is smaller.
            INDZ = IEP; RETURN
         ENDIF
      ENDDO
   ENDIF
ENDIF
END SUBROUTINE GET_Z_IND


END SUBROUTINE GET_CRTCFCC_INT_STENCILS



! -------------------------------- FILL_IJKO_INTERP_STENCILS ----------------------------

SUBROUTINE FILL_IJKO_INTERP_STENCILS

USE MPI_F08

! Local Variables:
INTEGER :: NM,NOM,N,IERR
TYPE (OMESH_TYPE), POINTER :: M2,M3
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ0,REQ0DUM
INTEGER :: N_REQ0
REAL(EB) CPUTIME, CPUTIME_START
LOGICAL :: PROCESS_SENDREC

IF (N_MPI_PROCESSES>1) ALLOCATE(REQ0(NMESHES))

N_REQ0 = 0

IF(GET_CUTCELLS_VERBOSE) THEN
   WRITE(LU_SETCC,'(A)',advance='no') '   > First loop, CC info..'
   CALL CPU_TIME(CPUTIME_START)
ENDIF

! Exchange number of cut-cells information to be exchanged between MESH and OMESHES:
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      PROCESS_SENDREC = .FALSE.
      DO N=1,MESHES(NM)%N_NEIGHBORING_MESHES
         IF (NOM==MESHES(NM)%NEIGHBORING_MESH(N)) PROCESS_SENDREC = .TRUE.
      ENDDO
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK  .AND. PROCESS_SENDREC) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%NFCC_S(1),2,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > Second loop, Barrier and isend..'
ENDIF

! DEFINITION NCCC_S:   MESHES(NOM)%OMESH(NM)%NFCC_S   = MESHES(NM)%OMESH(NOM)%NFCC_R

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   DO NOM=1,NMESHES
      PROCESS_SENDREC = .FALSE.
      DO N=1,MESHES(NOM)%N_NEIGHBORING_MESHES
         IF (NM==MESHES(NOM)%NEIGHBORING_MESH(N)) PROCESS_SENDREC = .TRUE.
      ENDDO
      IF (.NOT.PROCESS_SENDREC) CYCLE
      M3 => MESHES(NM)%OMESH(NOM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%NFCC_R(1),2,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%NFCC_S(1:2) = M3%NFCC_R(1:2)
      ENDIF
   ENDDO
ENDDO

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > MPI_WAITALL and Alloc..'
ENDIF

IF ( (N_REQ0>0) .AND. (N_MPI_PROCESSES>1) ) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

! At this point values of M2%NFCC_S should have been received.

! Definition: MESHES(NOM)%OMESH(NM)%IIO_FC_S(:) = MESHES(NM)%OMESH(NOM)%IIO_FC_R(:)
!             MESHES(NOM)%OMESH(NM)%JJO_FC_S(:) = MESHES(NM)%OMESH(NOM)%JJO_FC_R(:)
!             MESHES(NOM)%OMESH(NM)%KKO_FC_S(:) = MESHES(NM)%OMESH(NOM)%KKO_FC_R(:)
!             MESHES(NOM)%OMESH(NM)%AXS_FC_S(:) = MESHES(NM)%OMESH(NOM)%AXS_FC_R(:)

! Exchange list of face and cutcells data:
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (M2%NFCC_S(1)>0) THEN
         ALLOCATE(M2%IIO_FC_S(M2%NFCC_S(1)))
         ALLOCATE(M2%JJO_FC_S(M2%NFCC_S(1)))
         ALLOCATE(M2%KKO_FC_S(M2%NFCC_S(1)))
         ALLOCATE(M2%AXS_FC_S(M2%NFCC_S(1)))
      ENDIF
      IF (M2%NFCC_S(2)>0) THEN
         ALLOCATE(M2%IIO_CC_S(M2%NFCC_S(2)))
         ALLOCATE(M2%JJO_CC_S(M2%NFCC_S(2)))
         ALLOCATE(M2%KKO_CC_S(M2%NFCC_S(2)))
      ENDIF
   ENDDO
ENDDO

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > Faces non-blocking send-receives..'
ENDIF

! Faces:
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NFCC_S(1)>0) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%IIO_FC_S(1),M2%NFCC_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%JJO_FC_S(1),M2%NFCC_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%KKO_FC_S(1),M2%NFCC_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%AXS_FC_S(1),M2%NFCC_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NFCC_R(1)<1) CYCLE
      IF (PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%IIO_FC_R(1),M3%NFCC_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%JJO_FC_R(1),M3%NFCC_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%KKO_FC_R(1),M3%NFCC_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%AXS_FC_R(1),M3%NFCC_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%IIO_FC_S(1:M2%NFCC_S(1)) = M3%IIO_FC_R(1:M3%NFCC_R(1))
         M2%JJO_FC_S(1:M2%NFCC_S(1)) = M3%JJO_FC_R(1:M3%NFCC_R(1))
         M2%KKO_FC_S(1:M2%NFCC_S(1)) = M3%KKO_FC_R(1:M3%NFCC_R(1))
         M2%AXS_FC_S(1:M2%NFCC_S(1)) = M3%AXS_FC_R(1:M3%NFCC_R(1))
      ENDIF
   ENDDO
ENDDO

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > MPI_WAITALL Faces..'
ENDIF

IF ( (N_REQ0>0) .AND. (N_MPI_PROCESSES>1) ) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > Cells non-blocking send-receives..'
ENDIF

! Cells:
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NFCC_S(2)>0) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%IIO_CC_S(1),M2%NFCC_S(2),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%JJO_CC_S(1),M2%NFCC_S(2),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%KKO_CC_S(1),M2%NFCC_S(2),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NFCC_R(2)<1) CYCLE
      IF (PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%IIO_CC_R(1),M3%NFCC_R(2),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%JJO_CC_R(1),M3%NFCC_R(2),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%KKO_CC_R(1),M3%NFCC_R(2),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%IIO_CC_S(1:M2%NFCC_S(2)) = M3%IIO_CC_R(1:M3%NFCC_R(2))
         M2%JJO_CC_S(1:M2%NFCC_S(2)) = M3%JJO_CC_R(1:M3%NFCC_R(2))
         M2%KKO_CC_S(1:M2%NFCC_S(2)) = M3%KKO_CC_R(1:M3%NFCC_R(2))
      ENDIF
   ENDDO
ENDDO

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
   CALL CPU_TIME(CPUTIME_START)
   WRITE(LU_SETCC,'(A)',advance='no') '   > MPI_WAITALL Cells..'
ENDIF

IF ( (N_REQ0>0) .AND. (N_MPI_PROCESSES>1) ) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

IF(GET_CUTCELLS_VERBOSE) THEN
   CALL CPU_TIME(CPUTIME)
   WRITE(LU_SETCC,'(A,F8.3,A)') ' done. Time taken : ',CPUTIME-CPUTIME_START,' sec.'
ENDIF

IF(ALLOCATED(REQ0)) DEALLOCATE(REQ0)

RETURN

CONTAINS
SUBROUTINE CHECK_REQ0_SIZE
IF(N_REQ0>SIZE(REQ0,DIM=1)) THEN
   ALLOCATE(REQ0DUM(SIZE(REQ0,DIM=1)+NMESHES))
   REQ0DUM(1:N_REQ0-1) = REQ0(1:N_REQ0-1)
   CALL MOVE_ALLOC(REQ0DUM,REQ0)
ENDIF
END SUBROUTINE CHECK_REQ0_SIZE

END SUBROUTINE FILL_IJKO_INTERP_STENCILS



! --------------------------- GET_CLSPT_INBCF -----------------------------------

SUBROUTINE GET_CLSPT_INBCF(NM,XYZ,INBFC,INBFC_LOC,XYZ_IP,DIST,FOUNDPT,INSEG,INSEG2)

INTEGER,  INTENT(IN) :: NM, INBFC, INBFC_LOC
REAL(EB), INTENT(IN) :: XYZ(MAX_DIM)
REAL(EB), INTENT(OUT):: XYZ_IP(MAX_DIM), DIST
LOGICAL,  INTENT(OUT):: FOUNDPT, INSEG
LOGICAL,  OPTIONAL, INTENT(OUT):: INSEG2

! Local Variables:
INTEGER :: BODTRI(1:2),VERT_CUTFACE
INTEGER, ALLOCATABLE, DIMENSION(:) :: CFELEM
INTEGER :: X1AXIS,X2AXIS,X3AXIS
INTEGER :: IBOD,IWSEL,NVFACE,IPT,NVERT
REAL(EB):: NVEC(MAX_DIM),ANVEC(MAX_DIM),P0(MAX_DIM),A,B,C,D,PROJ_COEFF,XYZ_P(MAX_DIM)
REAL(EB):: PTCEN(IAXIS:JAXIS) !,AREAI,V1(IAXIS:JAXIS),V2(IAXIS:JAXIS)
REAL(EB):: SQRDIST, SQRDISTI, X2X3_1(IAXIS:JAXIS), X2X3_2(IAXIS:JAXIS)
REAL(EB):: DP(IAXIS:JAXIS),PCM1(IAXIS:JAXIS),PCM2(IAXIS:JAXIS),X2X3_IP(IAXIS:JAXIS)
REAL(EB):: T,DPDOTDP,SLOC,ATEST
REAL(EB):: P(IAXIS:KAXIS),DPP(IAXIS:KAXIS)
LOGICAL :: IN_POLY

! Initialize:
XYZ_IP(IAXIS:KAXIS) = 0._EB
DIST    = 1._EB / GEOMEPS
FOUNDPT = .FALSE.
INSEG   = .FALSE.
IF(PRESENT(INSEG2)) INSEG2=.FALSE.

VERT_CUTFACE = SIZE(MESHES(NM)%CUT_FACE(INBFC)%CFELEM, DIM=1); ALLOCATE(CFELEM(1:VERT_CUTFACE+1))
CFELEM(1:VERT_CUTFACE)  = MESHES(NM)%CUT_FACE(INBFC)%CFELEM(1:VERT_CUTFACE,INBFC_LOC)
BODTRI(1:2)  = MESHES(NM)%CUT_FACE(INBFC)%BODTRI(1:2,INBFC_LOC)

! normal vector to boundary surface triangle:
IBOD    = BODTRI(1)
IWSEL   = BODTRI(2)
NVEC(IAXIS:KAXIS)    = GEOMETRY(IBOD)%FACES_NORMAL(IAXIS:KAXIS,IWSEL)
NVFACE  = CFELEM(1);   CFELEM(NVFACE+2)=CFELEM(2)

! Plane equation for INBOUNDARY cut-face plane:
! Location of first point in cf polygon is P0:
IPT = 1
P0(IAXIS:KAXIS) = MESHES(NM)%CUT_FACE(INBFC)%XYZVERT(IAXIS:KAXIS,CFELEM(IPT+1))
A = NVEC(IAXIS)
B = NVEC(JAXIS)
C = NVEC(KAXIS)
D = -(A*P0(IAXIS) + B*P0(JAXIS) + C*P0(KAXIS))

! Project xyz point into plane of cf polygon:
PROJ_COEFF = (A*XYZ(IAXIS)+B*XYZ(JAXIS)+C*XYZ(KAXIS)) + D ! /dot(n,n) = 1
XYZ_P(IAXIS:KAXIS) = XYZ(IAXIS:KAXIS) - PROJ_COEFF*NVEC(IAXIS:KAXIS)

! Which Cartesian plane we project to?
ANVEC(IAXIS) = ABS(NVEC(IAXIS)); ANVEC(JAXIS) = ABS(NVEC(JAXIS)); ANVEC(KAXIS) = ABS(NVEC(KAXIS))
IF ( MAX(ANVEC(IAXIS),MAX(ANVEC(JAXIS),ANVEC(KAXIS))) == ANVEC(IAXIS) ) THEN
   X1AXIS = IAXIS; X2AXIS = JAXIS; X3AXIS = KAXIS
ELSEIF ( MAX(ANVEC(IAXIS),MAX(ANVEC(JAXIS),ANVEC(KAXIS))) == ANVEC(JAXIS) ) THEN
   X1AXIS = JAXIS; X2AXIS = KAXIS; X3AXIS = IAXIS
ELSE
   X1AXIS = KAXIS; X2AXIS = IAXIS; X3AXIS = JAXIS
ENDIF

! Now find closest point in projected plane:
! First: Test if point is inside cf area: Compute area of triangles formed
! by projected point xyz_p in x2,x3 plane and cf XYZvert, resp to
! CUT_FACE.area*nvec(x1axis), i.e. cut-face area projected on the Cartesian
! plane orthogonal to x1axis:
PTCEN(IAXIS:JAXIS) = XYZ_P( (/ X2AXIS, X3AXIS /) )
NVERT = SIZE(MESHES(NM)%CUT_FACE(INBFC)%XYZVERT,DIM=2)
ATEST = MESHES(NM)%CUT_FACE(INBFC)%AREA(INBFC_LOC)*ANVEC(X1AXIS) ! Test Area is projected area into X1AXIS.
CALL POINT_IN_POLYGON(PTCEN,VERT_CUTFACE+1,CFELEM,NVERT,X2AXIS,X3AXIS,MESHES(NM)%CUT_FACE(INBFC)%XYZVERT,IN_POLY)

! Test if inside:
IF (IN_POLY) THEN
   ! Same areas, xyz_p inside INBOUNDARY cut-face:
   XYZ_IP(IAXIS:KAXIS) = XYZ_P(IAXIS:KAXIS)
   DIST   = SQRT( (XYZ(IAXIS)-XYZ_IP(IAXIS))**2._EB + &
                  (XYZ(JAXIS)-XYZ_IP(JAXIS))**2._EB + &
                  (XYZ(KAXIS)-XYZ_IP(KAXIS))**2._EB )
   FOUNDPT= .TRUE.
   ! Now check if point is in segment:
   IF (PRESENT(INSEG2)) THEN
      DO IPT=1,NVFACE
          P(IAXIS:KAXIS)  = MESHES(NM)%CUT_FACE(INBFC)%XYZVERT(IAXIS:KAXIS,CFELEM(IPT+1))
          DPP(IAXIS:KAXIS)= MESHES(NM)%CUT_FACE(INBFC)%XYZVERT(IAXIS:KAXIS,CFELEM(IPT+2))-P(IAXIS:KAXIS)
          IF (NORM2(DPP(IAXIS:KAXIS)) < TWENTY_EPSILON_EB) CYCLE
          DPP(IAXIS:KAXIS)=DPP(IAXIS:KAXIS)/NORM2(DPP(IAXIS:KAXIS))
          P = XYZ_IP - (P + DOT_PRODUCT(DPP,XYZ_IP-P)*DPP)
          IF (NORM2(P(IAXIS:KAXIS)) < GEOMEPS) THEN
             INSEG = .TRUE.
             INSEG2= .TRUE.
             EXIT
          ENDIF
      ENDDO
   ENDIF
   DEALLOCATE(CFELEM)
   RETURN
ENDIF

! Second, test against segments: Find closest point in segments in x2,x3 plane:
SQRDIST = 1._EB / GEOMEPS
DO IPT=1,NVFACE

    X2X3_1(IAXIS:JAXIS) = MESHES(NM)%CUT_FACE(INBFC)%XYZVERT((/ X2AXIS, X3AXIS /) ,CFELEM(IPT+1))
    X2X3_2(IAXIS:JAXIS) = MESHES(NM)%CUT_FACE(INBFC)%XYZVERT((/ X2AXIS, X3AXIS /) ,CFELEM(IPT+2))

    ! Smallest distance from point PC to segment x2x3_1-x2x3_2:
    DP(IAXIS:JAXIS)     = X2X3_2(IAXIS:JAXIS) - X2X3_1(IAXIS:JAXIS)
    PCM1(IAXIS:JAXIS)   =  PTCEN(IAXIS:JAXIS) - X2X3_1(IAXIS:JAXIS)
    T      = DP(IAXIS)*PCM1(IAXIS) + DP(JAXIS)*PCM1(JAXIS)
    DPDOTDP= DP(IAXIS)**2._EB + DP(JAXIS)**2._EB

    IF ( T < GEOMEPS ) THEN
        SQRDISTI = PCM1(IAXIS)**2._EB + PCM1(JAXIS)**2._EB ! x2x3_1 is closest pt.
        T = 0._EB
    ELSEIF ( T >= DPDOTDP ) THEN
        PCM2(IAXIS:JAXIS) = PTCEN(IAXIS:JAXIS) - X2X3_2(IAXIS:JAXIS)
        SQRDISTI = PCM2(IAXIS)**2._EB + PCM2(JAXIS)**2._EB ! x2x3_2 is closest pt.
        T = DPDOTDP
    ELSE
        SQRDISTI =(PCM1(IAXIS)**2._EB + PCM1(JAXIS)**2._EB) - T**2._EB/DPDOTDP
    ENDIF

    ! Test:
    IF ( SQRDISTI < SQRDIST ) THEN
        SQRDIST = SQRDISTI
        SLOC    = T/(DPDOTDP+TWENTY_EPSILON_EB)
        X2X3_IP(IAXIS:JAXIS) = X2X3_1(IAXIS:JAXIS) + SLOC * DP(IAXIS:JAXIS) ! intersection point in segment,
                                                                            ! plane x2,x3
        FOUNDPT= .TRUE.
        INSEG  = .TRUE.
    ENDIF
ENDDO

! Now pass x2x3 intersection point to 3D:
IF (FOUNDPT) THEN
    SELECT CASE(X1AXIS)
        CASE(IAXIS)
            XYZ_IP(JAXIS) = X2X3_IP(IAXIS)
            XYZ_IP(KAXIS) = X2X3_IP(JAXIS)
            XYZ_IP(IAXIS) = (-B*XYZ_IP(JAXIS) -C*XYZ_IP(KAXIS) - D)/A
        CASE(JAXIS)
            XYZ_IP(KAXIS) = X2X3_IP(IAXIS)
            XYZ_IP(IAXIS) = X2X3_IP(JAXIS)
            XYZ_IP(JAXIS) = (-A*XYZ_IP(IAXIS) -C*XYZ_IP(KAXIS) - D)/B
        CASE(KAXIS)
            XYZ_IP(IAXIS) = X2X3_IP(IAXIS)
            XYZ_IP(JAXIS) = X2X3_IP(JAXIS)
            XYZ_IP(KAXIS) = (-A*XYZ_IP(IAXIS) -B*XYZ_IP(JAXIS) - D)/C
    END SELECT
    DIST   = SQRT( (XYZ(IAXIS)-XYZ_IP(IAXIS))**2._EB + &
                   (XYZ(JAXIS)-XYZ_IP(JAXIS))**2._EB + &
                   (XYZ(KAXIS)-XYZ_IP(KAXIS))**2._EB )
ENDIF

DEALLOCATE(CFELEM)

RETURN
END SUBROUTINE GET_CLSPT_INBCF



! -------------------------- GET_CLOSEPT_CCVT -----------------------------------

SUBROUTINE GET_CLOSEPT_CCVT(NM,XYZ,ICC,XYZ_IP,DIST,FOUNDPT,IFCPT,IFCPT_LOC)

INTEGER,  INTENT(IN) :: NM, ICC
REAL(EB), INTENT(IN) :: XYZ(MAX_DIM)
REAL(EB), INTENT(OUT):: XYZ_IP(MAX_DIM), DIST
INTEGER,  INTENT(OUT):: IFCPT, IFCPT_LOC
LOGICAL,  INTENT(OUT):: FOUNDPT

! Local Variables:
INTEGER :: I,J,K,IJK(MAX_DIM),IJK_CELL(MAX_DIM),LOWHIGH,IND_ADD
INTEGER :: X1AXIS,X2AXIS,X3AXIS,XIAXIS,XJAXIS,XKAXIS,CEI,ICF,IX2,IX3
LOGICAL :: INLIST, ISCORN
INTEGER :: INDXI(MAX_DIM),IPT,IVERT,ICORN,INDI,INDJ,INDK
REAL(EB), POINTER, DIMENSION(:) :: X2FC,X3FC!,X1FC,X1CL,X2CL,X3CL,DX1FC,DX2FC,DX3FC,DX1CL,DX2CL,DX3CL
REAL(EB):: DV(MAX_DIM),XY1(1:4,IAXIS:JAXIS)
INTEGER :: JJ,KK,INDXI1(IAXIS:JAXIS),INDXI2(IAXIS:JAXIS),INDXI3(IAXIS:JAXIS),INDXI4(IAXIS:JAXIS)
LOGICAL :: CEIFLG

! Initialize:
XYZ_IP(IAXIS:KAXIS) = 0._EB
DIST    = 1._EB / GEOMEPS
FOUNDPT = .FALSE.
IFCPT   = 0; IFCPT_LOC = 0

! Here we need to look at Cartesian faces that are boundary of
! CUT_CELL(icc) (which has only regular or Gasphase cut-faces of regular
! size and find if a corner point (or sigular point) is type SOLID. The
! point found provides xyz_ip:
IJK_CELL(IAXIS:KAXIS) = MESHES(NM)%CUT_CELL(ICC)%IJK(IAXIS:KAXIS)

! Loop on different planes:
LOWHIGH_IND_LOOP : DO LOWHIGH=LOW_IND,HIGH_IND

   IND_ADD = LOWHIGH - LOW_IND  ! Index to add for face LOW-HIGH resp to cell.

   X1AXIS_LOOP : DO X1AXIS=IAXIS,KAXIS

      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         I = IJK_CELL(IAXIS)-1+IND_ADD
         J = IJK_CELL(JAXIS)
         K = IJK_CELL(KAXIS)
         X2AXIS = JAXIS; X3AXIS = KAXIS
         ! location in I,J,K of x1,x2,x3 axes:
         XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
         ! Centroid coordinates in x1,x2,x3 axes:
         ! X2CL => YCELL; DX2CL => DYCELL
         ! X3CL => ZCELL; DX3CL => DZCELL
         ! X1FC => XFACE
         X2FC => YFACE
         X3FC => ZFACE

      CASE(JAXIS)
         I = IJK_CELL(IAXIS)
         J = IJK_CELL(JAXIS)-1+IND_ADD
         K = IJK_CELL(KAXIS)
         X2AXIS = KAXIS; X3AXIS = IAXIS
         ! location in I,J,K of x1,x2,x3 axes:
         XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
         ! Centroid coordinates in x1,x2,x3 axes:
         ! X2CL => ZCELL; DX2CL => DZCELL
         ! X3CL => XCELL; DX3CL => DXCELL
         ! X1FC => YFACE;
         X2FC => ZFACE;
         X3FC => XFACE;
      CASE(KAXIS)

         I = IJK_CELL(IAXIS)
         J = IJK_CELL(JAXIS)
         K = IJK_CELL(KAXIS)-1+IND_ADD
         X2AXIS = IAXIS; X3AXIS = JAXIS
         ! location in I,J,K of x1,x2,x3 axes:
         XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
         ! Face coordinates in x1,x2,x3 axes:
         ! X2CL => XCELL; DX2CL => DXCELL
         ! X3CL => YCELL; DX3CL => DYCELL
         ! X1FC => ZFACE
         X2FC => XFACE
         X3FC => YFACE

      END SELECT

      ! Drop if face is regular GASPHASE or SOLID:
      IF ( MESHES(NM)%FCVAR(I,J,K,CC_FGSC,X1AXIS) /= CC_CUTCFE ) CYCLE

      ! Face IJK:
      IJK(IAXIS:KAXIS) = (/ I, J, K /)

      ! Cartesian Face centroid location x2-x3 plane:
      CEI = MESHES(NM)%FCVAR(I,J,K,CC_IDCE,X1AXIS)
      ICF = MESHES(NM)%FCVAR(I,J,K,CC_IDCF,X1AXIS)

      ! We might have single point CEIs:
      CEIFLG = .FALSE.
      IF (CEI <= 0) THEN
         CEIFLG = .TRUE.
      ELSEIF ( MESHES(NM)%CUT_EDGE(CEI)%NVERT == 1 ) THEN
         CEIFLG = .TRUE.
      ENDIF

      IF (CEIFLG) THEN ! Cut face is one regular face, with one SOLID vertex.

         ! Figure out which vertex is SOLID and location:
         INLIST = .FALSE.
         DO IX3=0,1
            DO IX2=0,1
               ! Vertex axes:
               INDXI(IAXIS:KAXIS) = (/ IJK(X1AXIS), IJK(X2AXIS)-1+IX2, IJK(X3AXIS)-1+IX3 /) ! x1,x2,x3
               INDI = INDXI(XIAXIS)
               INDJ = INDXI(XJAXIS)
               INDK = INDXI(XKAXIS)
               IF ( MESHES(NM)%VERTVAR(INDI,INDJ,INDK,CC_VGSC) == CC_SOLID ) THEN
                  INLIST = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (INLIST) THEN
               XYZ_IP(IAXIS:KAXIS)  = (/ XFACE(INDI), YFACE(INDJ), ZFACE(INDK) /)
               DIST   = SQRT( (XYZ(IAXIS)-XYZ_IP(IAXIS))**2._EB + &
                              (XYZ(JAXIS)-XYZ_IP(JAXIS))**2._EB + &
                              (XYZ(KAXIS)-XYZ_IP(KAXIS))**2._EB )
               IFCPT  = ICF
               ! Find local point:
               DO IPT=1,MESHES(NM)%CUT_FACE(ICF)%NVERT
                  DV = MESHES(NM)%CUT_FACE(ICF)%XYZVERT(IAXIS:KAXIS,IPT) - XYZ_IP(IAXIS:KAXIS)
                  IF( (ABS(DV(IAXIS))+ABS(DV(JAXIS))+ABS(DV(KAXIS))) < GEOMEPS ) THEN
                     IFCPT_LOC = IPT
                     EXIT
                  ENDIF
               ENDDO
               FOUNDPT = .TRUE.
               RETURN
            ENDIF
         ENDDO

         ! Check if there are more than 4 vertices, get vertex that is not in corners:
         IF ( MESHES(NM)%CUT_FACE(ICF)%NVERT > 4 ) THEN

            JJ = IJK(X2AXIS); KK = IJK(X3AXIS)
            ! Vertex at index jj-1,kk-1:
            INDXI1(IAXIS:JAXIS) = (/ JJ-1  , KK-1   /) ! Local x2,x3
            ! Vertex at index jj,kk-1:
            INDXI2(IAXIS:JAXIS) = (/ JJ    , KK-1   /) ! Local x2,x3
            ! Vertex at index jj,kk:
            INDXI3(IAXIS:JAXIS) = (/ JJ    , KK     /) ! Local x2,x3
            ! Vertex at index jj-1,kk:
            INDXI4(IAXIS:JAXIS) = (/ JJ-1  , KK     /) ! Local x2,x3

            XY1(1:4,IAXIS) = (/ X2FC(INDXI1(IAXIS)), X2FC(INDXI2(IAXIS)), &
                                X2FC(INDXI3(IAXIS)), X2FC(INDXI4(IAXIS)) /)
            XY1(1:4,JAXIS) = (/ X3FC(INDXI1(JAXIS)), X3FC(INDXI2(JAXIS)), &
                                X3FC(INDXI3(JAXIS)), X3FC(INDXI4(JAXIS)) /)

            ! Find vertex:
            DO IVERT=1,MESHES(NM)%CUT_FACE(ICF)%NVERT
               ISCORN = .FALSE.
               DO ICORN=1,4
                  IF( SQRT( (XY1(ICORN,IAXIS)-MESHES(NM)%CUT_FACE(ICF)%XYZVERT(X2AXIS,IVERT))**2._EB + &
                            (XY1(ICORN,JAXIS)-MESHES(NM)%CUT_FACE(ICF)%XYZVERT(X3AXIS,IVERT))**2._EB ) &
                            < GEOMEPS) THEN
                     ISCORN = .TRUE.
                     EXIT
                  ENDIF
               ENDDO
               IF (.NOT.ISCORN) THEN
                  XYZ_IP(IAXIS:KAXIS) = MESHES(NM)%CUT_FACE(ICF)%XYZVERT(IAXIS:KAXIS,IVERT)
                  INLIST = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (INLIST) THEN
               DIST   = SQRT( (XYZ(IAXIS)-XYZ_IP(IAXIS))**2._EB + &
                              (XYZ(JAXIS)-XYZ_IP(JAXIS))**2._EB + &
                              (XYZ(KAXIS)-XYZ_IP(KAXIS))**2._EB )
               IFCPT     = ICF
               IFCPT_LOC = IVERT
               FOUNDPT   = .TRUE.
               RETURN
            ENDIF

         ENDIF

      ENDIF ! CEI <= 0

      !NULLIFY(X2CL,X3CL,DX2CL,DX3CL,X1FC,X2FC,X3FC)
      NULLIFY(X2FC,X3FC)

   ENDDO X1AXIS_LOOP

ENDDO LOWHIGH_IND_LOOP


RETURN
END SUBROUTINE GET_CLOSEPT_CCVT



! ----------------------- SET_CC_MATVEC_DATA ---------------------------------

SUBROUTINE SET_CC_MATVEC_DATA

USE MPI_F08

! Local variables:
INTEGER :: NM, I, IPROC, IERR

! 1. Define unknown numbers for Scalars:
CALL GET_LINKED_MATRIX_INDEXES_Z

! 2. For each CC_GASPHASE (cut or regular) face, find global numeration of the volumes
! that share it, store a list of areas and centroids for diffussion operator in FV form.
! 3. Get CC_GASPHASE regular faces data, for scalars Z:
CALL GET_GASPHASE_REGRCFACES_DATA

! 4. Get CC_GASPHASE cut-faces data:
CALL GET_GASPHASE_CUTFACES_DATA ! Here there is no need to populate CELL_LIST on CUT_FACE,
                                ! list of low/high cut-cell volumes that share the cut-face, as
                                ! this has been done before calling SET_CC_MATVEC_DATA, when calling
                                ! GET_CRTCFCC_INT_STENCILS.

! 6. Exchange information at block boundaries for RC_FACE, CUT_FACE
! fields on each mesh:
CALL DEFINE_SHARED_FACES
CALL GET_RCEDGE_FACE_LIST ! This routine makes use of RCFACES index defined in ECVAR(:,:,:,CC_IDRC,:) in the previous.

! 6.5 Link faces for momentum transport if unstructured projection:
CALL GET_LINKED_FACE_INDEXES_F

! 7. Get nonzeros graph of the scalar diffusion/advection matrix, defined as:
!    - NNZ_D_MAT_Z(1:NUNKZ_LOCAL) Number of nonzeros on per matrix row.
!    - JD_MAT_Z(1:NNZ_ROW_Z,1:NUNKZ_LOCAL) Column location of nonzeros, global numeration.
NUNKZ_LOCAL = sum(NUNKZ_LOC(1:NMESHES)) ! Filled in GET_MATRIX_INDEXES, only nonzeros are for meshes
                                        ! that belong to this process.
NUNKZ_TOTAL = sum(NUNKZ_TOT(1:NMESHES))

IF (GET_CUTCELLS_VERBOSE) THEN
   IF (MY_RANK==0) THEN
      WRITE(LU_ERR,*) ' '
      WRITE(LU_ERR,'(A)') ' Cut-cell region scalar transport advanced explicitly.'
      WRITE(LU_ERR,'(A)') ' List of Scalar unknown numbers per proc:'
   ENDIF
   DO IPROC=0,N_MPI_PROCESSES-1
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      IF(MY_RANK==IPROC) WRITE(LU_ERR,'(A,I8,A,I8)') ' MY_RANK=',MY_RANK,', NUNKZ_LOCAL=',NUNKZ_LOCAL
   ENDDO
ENDIF

! Allocate NNZ_D_MAT_Z, JD_MAT_Z:
ALLOCATE( NNZ_D_MAT_Z(1:NUNKZ_LOCAL) )
ALLOCATE( JD_MAT_Z(1:NNZ_ROW_Z,1:NUNKZ_LOCAL) ) ! Contains on first index nonzeros per local row.
NNZ_D_MAT_Z(:) = 0
JD_MAT_Z(:,:)  = HUGE(I)

! Find NM_START: first mesh that belongs to the processor.
NM_START = CC_UNDEFINED
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   NM_START = NM
   EXIT
ENDDO

! 8. Build Mass (volumes) matrix for scalars:
CALL GET_MMATRIX_SCALAR_3D

! Allocate rhs and solution arrays for species:
ALLOCATE( F_Z(1:NUNKZ_LOCAL) , F_Z0(1:NUNKZ_LOCAL,1:N_TOTAL_SCALARS) , RZ_Z(1:NUNKZ_LOCAL) , RZ_ZS(1:NUNKZ_LOCAL) )
ALLOCATE( RZ_Z0(1:NUNKZ_LOCAL,1:N_TOTAL_SCALARS) )
ALLOCATE( P_0_CV(1:NUNKZ_LOCAL), TMP_0_CV(1:NUNKZ_LOCAL), RHO_0_CV(1:NUNKZ_LOCAL), ZCEN_CV(1:NUNKZ_LOCAL) )

RETURN
END SUBROUTINE SET_CC_MATVEC_DATA



! ------------------------- GET_RCEDGE_FACE_LIST --------------------------------

SUBROUTINE GET_RCEDGE_FACE_LIST

INTEGER :: NM,ICF,JCF,I,J,K,X1AXIS,EAXIS,LOHI_AXIS,IRC,ISIDE,HILO
TYPE(MESH_TYPE), POINTER :: M

MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M=>MESHES(NM)
   ! X axis edges:
   EAXIS = IAXIS
   DO K=0,M%KBAR
      DO J=0,M%JBAR
         DO I=1,M%IBAR
            IRC = M%ECVAR(I,J,K,CC_IDCE,EAXIS)
            IF (M%ECVAR(I,J,K,CC_EGSC,EAXIS)/=CC_GASPHASE .OR. IRC<1) CYCLE
            ALLOCATE(M%CC_RCEDGE(IRC)%FACE_LIST(1:3,-2:2)); M%CC_RCEDGE(IRC)%FACE_LIST = CC_UNDEFINED
            ! Faces -1 and 1 : KAXIS, I,J,K and KAXIS I,J+1,K
            X1AXIS = KAXIS; LOHI_AXIS = JAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I,J+ISIDE,K,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I,J+ISIDE,K,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I,J+ISIDE,K,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I,J+ISIDE,K,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
            ! Faces -2 and 2 : JAXIS, I,J,K and JAXIS I,J,K+1
            X1AXIS = JAXIS; LOHI_AXIS = KAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I,J,K+ISIDE,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I,J,K+ISIDE,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I,J,K+ISIDE,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I,J,K+ISIDE,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Y axis edges:
   EAXIS = JAXIS
   DO K=0,M%KBAR
      DO J=1,M%JBAR
         DO I=0,M%IBAR
            IRC = M%ECVAR(I,J,K,CC_IDCE,EAXIS)
            IF (M%ECVAR(I,J,K,CC_EGSC,EAXIS)/=CC_GASPHASE .OR. IRC<1) CYCLE
            ALLOCATE(M%CC_RCEDGE(IRC)%FACE_LIST(1:3,-2:2)); M%CC_RCEDGE(IRC)%FACE_LIST = CC_UNDEFINED
            ! Faces -1 and 1 : IAXIS, I,J,K and IAXIS I,J,K+1
            X1AXIS = IAXIS; LOHI_AXIS = KAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I,J,K+ISIDE,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I,J,K+ISIDE,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I,J,K+ISIDE,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I,J,K+ISIDE,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
            ! Faces -2 and 2 : KAXIS, I,J,K and KAXIS I+1,J,K
            X1AXIS = KAXIS; LOHI_AXIS = IAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I+ISIDE,J,K,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I+ISIDE,J,K,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I+ISIDE,J,K,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I+ISIDE,J,K,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! Z axis edges:
   EAXIS = KAXIS
   DO K=1,M%KBAR
      DO J=0,M%JBAR
         DO I=0,M%IBAR
            IRC = M%ECVAR(I,J,K,CC_IDCE,EAXIS)
            IF (M%ECVAR(I,J,K,CC_EGSC,EAXIS)/=CC_GASPHASE .OR. IRC<1) CYCLE
            ALLOCATE(M%CC_RCEDGE(IRC)%FACE_LIST(1:3,-2:2)); M%CC_RCEDGE(IRC)%FACE_LIST = CC_UNDEFINED
            ! Faces -1 and 1 : JAXIS, I,J,K and JAXIS I+1,J,K
            X1AXIS = JAXIS; LOHI_AXIS = IAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I+ISIDE,J,K,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I+ISIDE,J,K,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I+ISIDE,J,K,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I+ISIDE,J,K,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,2*ISIDE-1) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
            ! Faces -2 and 2 : IAXIS, I,J,K and IAXIS I,J+1,K
            X1AXIS = IAXIS; LOHI_AXIS = JAXIS
            DO ISIDE=0,1
               IF(M%FCVAR(I,J+ISIDE,K,CC_IDCF,X1AXIS)>0) THEN ! CUT_FACE
                  ICF = M%FCVAR(I,J+ISIDE,K,CC_IDCF,X1AXIS)
                  ! Find which cut-face in ICF entry has the RCEDGE(IRC) as boundary:
                  HILO=2-ISIDE ! 2=HIGH_IND side reg edge for face; 1=LOW_IND reg side edge.
                  CALL GET_RCEDGE_CUTFACE(NM,LOHI_AXIS,HILO,ICF,JCF)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_CFGAS, ICF, JCF /)
               ELSEIF (M%FCVAR(I,J+ISIDE,K,CC_IDRC,X1AXIS)>0) THEN ! RCFACE
                  ICF = M%FCVAR(I,J+ISIDE,K,CC_IDRC,X1AXIS)
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RCGAS, ICF,   0 /)
               ELSE ! Regular gas face. There can be no solid faces adjacent to an RCEDGE.
                  M%CC_RCEDGE(IRC)%FACE_LIST(1:3,4*ISIDE-2) = (/ CC_FTYPE_RGGAS,   0,   0 /)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! DO IRC=1,M%CC_NRCEDGE_Z
   !    WRITE(LU_ERR,*) ' '
   !    WRITE(LU_ERR,*) 'IRC=',IRC,M%CC_RCEDGE(IRC)%IJK(1:4)
   !    WRITE(LU_ERR,*) '-2 =',M%CC_RCEDGE(IRC)%FACE_LIST(1:3,-2)
   !    WRITE(LU_ERR,*) '-1 =',M%CC_RCEDGE(IRC)%FACE_LIST(1:3,-1)
   !    WRITE(LU_ERR,*) ' 1 =',M%CC_RCEDGE(IRC)%FACE_LIST(1:3, 1)
   !    WRITE(LU_ERR,*) ' 2 =',M%CC_RCEDGE(IRC)%FACE_LIST(1:3, 2)
   ! ENDDO

ENDDO MESH_LOOP

END SUBROUTINE GET_RCEDGE_FACE_LIST



! ------------------------------ GET_RCEDGE_CUTFACE -----------------------------------

SUBROUTINE GET_RCEDGE_CUTFACE(NM,AXIS,SIDE,ICF,JCF)

INTEGER, INTENT(IN) :: NM,AXIS,SIDE,ICF
INTEGER, INTENT(OUT):: JCF

INTEGER :: IED,IEDGE
TYPE(CC_CUTFACE_TYPE), POINTER :: CF

CF=> MESHES(NM)%CUT_FACE(ICF)
DO JCF=1,CF%NFACE
   DO IED=2,CF%CEDGES(1,JCF)+1
      IEDGE=CF%CEDGES(IED,JCF)
      IF(CF%EDGE_LIST(1,IEDGE)/=CC_ETYPE_RGGAS) CYCLE  ! Regular gas edge, RCEDGES defined like this in
                                                        ! GET_CARTFACE_CUTFACES.
      IF(CF%EDGE_LIST(2,IEDGE)/=SIDE)            CYCLE  ! Edge in low or high side of face.
      IF(CF%EDGE_LIST(3,IEDGE)/=AXIS)            CYCLE  ! Direction in which edge is low or high resp to face.
      RETURN
   ENDDO
ENDDO

END SUBROUTINE GET_RCEDGE_CUTFACE



! ------------------------- DEFINE_SHARED_FACES --------------------------------
SUBROUTINE DEFINE_SHARED_FACES

! Once cells have been linked, regular and RC faces shared by linked cells are defined as SHARED=.TRUE.

INTEGER :: NM,ICF,JCF,I,J,K,X1AXIS

DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   ! First cut-faces:
   DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
      IF(CUT_FACE(ICF)%STATUS/=CC_GASPHASE) CYCLE
      CUT_FACE(ICF)%SHARED = .FALSE.
      DO JCF=1,CUT_FACE(ICF)%NFACE
         IF(CUT_FACE(ICF)%UNKZ(LOW_IND,JCF)==CUT_FACE(ICF)%UNKZ(HIGH_IND,JCF)) &
         CUT_FACE(ICF)%SHARED(JCF) = .TRUE. ! This face is shared by two linked cut-cells.
      ENDDO
   ENDDO
   ! Then RC faces:
   FCVAR(:,:,:,CC_IDRC,:) = 0
   DO ICF=1,MESHES(NM)%CC_NRCFACE_Z
      ! Add ICF position in FCVAR(I,J,K,CC_IDRC,X1AXIS):
      I      = RC_FACE(ICF)%IJK(IAXIS)
      J      = RC_FACE(ICF)%IJK(JAXIS)
      K      = RC_FACE(ICF)%IJK(KAXIS)
      X1AXIS = RC_FACE(ICF)%IJK(KAXIS+1)
      FCVAR(I,J,K,CC_IDRC,X1AXIS) = ICF
      ! Then set if SHARED:
      IF(RC_FACE(ICF)%UNKZ(LOW_IND)==RC_FACE(ICF)%UNKZ(HIGH_IND)) RC_FACE(ICF)%SHAREDZ = .TRUE.
   ENDDO

ENDDO

RETURN
END SUBROUTINE DEFINE_SHARED_FACES


! ----------------------- FILL_UNKZ_GUARDCELLS ---------------------------------
SUBROUTINE FILL_UNKZ_GUARDCELLS

USE MPI_F08

! Local Variables:
INTEGER :: NM,NOM,IERR
TYPE (MESH_TYPE), POINTER :: M
TYPE (OMESH_TYPE), POINTER :: M2,M3
TYPE (MPI_REQUEST), ALLOCATABLE, DIMENSION(:) :: REQ0,REQ0DUM
INTEGER :: N_REQ0, NICC_R, ICC, ICC1, NCELL, JCC, NICF_R, ICF
INTEGER, ALLOCATABLE, DIMENSION(:) :: NCC_SV
INTEGER :: ISTR,IEND,JSTR,JEND,KSTR,KEND,IIO,JJO,KKO,IOR,IW,N_INT,IIOF,JJOF,KKOF,X1AXIS
LOGICAL :: ALL_FLG
INTEGER, ALLOCATABLE, DIMENSION(:) :: INT1D


! First allocate buffers to receive UNKZ information:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      IF (MESHES(NM)%OMESH(NOM)%NIC_R>0) THEN
         M3 => MESHES(NM)%OMESH(NOM)
         ALLOCATE(M3%UNKZ_CT_R(M3%NIC_R))
      ENDIF
      IF (MESHES(NM)%OMESH(NOM)%NICC_R(1)>0) THEN
         M3 => MESHES(NM)%OMESH(NOM)
         ALLOCATE(M3%ICC_UNKZ_CC_R(3*M3%NICC_R(1))); M3%ICC_UNKZ_CC_R = CC_UNDEFINED
         ALLOCATE(M3%UNKZ_CC_R(M3%NICC_R(2)))

         ! Dump cut-cell indexes on sending mesh NOM, whose info will be received:
         NICC_R = 0
         CALL POINT_TO_MESH(NM)
         ! Loop over cut-cells:
         EXTERNAL_WALL_LOOP_1A : DO IW=1,N_EXTERNAL_WALL_CELLS
            WC=>WALL(IW)
            IF (.NOT.ANY(WC%BOUNDARY_TYPE==(/INTERPOLATED_BOUNDARY,MIRROR_BOUNDARY/))) CYCLE EXTERNAL_WALL_LOOP_1A
            EWC=>EXTERNAL_WALL(IW); BC=>BOUNDARY_COORD(WC%BC_INDEX)
            IF ( WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY ) THEN
               IF (EWC%NOM/=NOM) CYCLE EXTERNAL_WALL_LOOP_1A
               IF (CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1A
               DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
                  DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                     DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                       ICC   = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
                       IF (ICC > 0) THEN
                          NICC_R = NICC_R + 1
                          M3%ICC_UNKZ_CC_R(3*NICC_R-2:3*NICC_R) = (/ IIO, JJO, KKO /) ! Note : This ICC index refers to NOM mesh.
                       ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ELSEIF ( WC%BOUNDARY_TYPE==MIRROR_BOUNDARY ) THEN
               IF (NM/=NOM) CYCLE EXTERNAL_WALL_LOOP_1A
               IF (CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1A
               IIO = BC%IIG; JJO = BC%JJG; KKO = BC%KKG; IOR = BC%IOR
               ! CYCLE if OBJECT face is in the Mirror Boundary, normal out into ghost-cell:
               SELECT CASE(IOR)
               CASE( IAXIS); IF(FCVAR(IIO-1,JJO  ,KKO  ,CC_FGSC,IAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               CASE(-IAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,IAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               CASE( JAXIS); IF(FCVAR(IIO  ,JJO-1,KKO  ,CC_FGSC,JAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               CASE(-JAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,JAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               CASE( KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO-1,CC_FGSC,KAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               CASE(-KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,KAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_1A
               END SELECT
               ICC = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC); IF (ICC<1) CYCLE
               NICC_R = NICC_R + 1
               M3%ICC_UNKZ_CC_R(3*NICC_R-2:3*NICC_R) = (/ IIO, JJO, KKO /) ! Note : This ICC index refers to NOM==NM mesh.
            ENDIF
         ENDDO EXTERNAL_WALL_LOOP_1A
      ENDIF

      ! Count and add interpolated boundary cut-faces:
      MESHES(NM)%OMESH(NOM)%NICF_R(1) = 0
      MESHES(NM)%OMESH(NOM)%NICF_R(2) = 0
      CALL POINT_TO_MESH(NM)
      ! Loop over cut-cells:
      DO IW=1,N_EXTERNAL_WALL_CELLS
         WC=>WALL(IW); IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
         EWC=>EXTERNAL_WALL(IW); IF (EWC%NOM/=NOM) CYCLE
         BC=>BOUNDARY_COORD(WC%BC_INDEX); IIO = BC%IIG; JJO = BC%JJG; KKO = BC%KKG
         SELECT CASE(BC%IOR)
         CASE( IAXIS); IF(FCVAR(IIO-1,JJO  ,KKO  ,CC_FGSC,IAXIS) /= CC_CUTCFE) CYCLE
         CASE(-IAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,IAXIS) /= CC_CUTCFE) CYCLE
         CASE( JAXIS); IF(FCVAR(IIO  ,JJO-1,KKO  ,CC_FGSC,JAXIS) /= CC_CUTCFE) CYCLE
         CASE(-JAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,JAXIS) /= CC_CUTCFE) CYCLE
         CASE( KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO-1,CC_FGSC,KAXIS) /= CC_CUTCFE) CYCLE
         CASE(-KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,KAXIS) /= CC_CUTCFE) CYCLE
         END SELECT
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                SELECT CASE(-BC%IOR)
                CASE( IAXIS); ICF=MESHES(NOM)%FCVAR(IIO-1,JJO  ,KKO  ,CC_IDCF,IAXIS)
                CASE(-IAXIS); ICF=MESHES(NOM)%FCVAR(IIO  ,JJO  ,KKO  ,CC_IDCF,IAXIS)
                CASE( JAXIS); ICF=MESHES(NOM)%FCVAR(IIO  ,JJO-1,KKO  ,CC_IDCF,JAXIS)
                CASE(-JAXIS); ICF=MESHES(NOM)%FCVAR(IIO  ,JJO  ,KKO  ,CC_IDCF,JAXIS)
                CASE( KAXIS); ICF=MESHES(NOM)%FCVAR(IIO  ,JJO  ,KKO-1,CC_IDCF,KAXIS)
                CASE(-KAXIS); ICF=MESHES(NOM)%FCVAR(IIO  ,JJO  ,KKO  ,CC_IDCF,KAXIS)
                END SELECT
                IF (ICF > 0) THEN
                   MESHES(NM)%OMESH(NOM)%NICF_R(1) = MESHES(NM)%OMESH(NOM)%NICF_R(1) + 1
                   MESHES(NM)%OMESH(NOM)%NICF_R(2) = MESHES(NM)%OMESH(NOM)%NICF_R(2) + &
                                                     MESHES(NOM)%CUT_FACE(ICF)%NFACE
                ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (MESHES(NM)%OMESH(NOM)%NICF_R(1)>0) THEN
         M3 => MESHES(NM)%OMESH(NOM)
         ALLOCATE(M3%ICF_UFFB_CF_R(4*M3%NICF_R(1)))
         ! Dump cut-cell indexes on sending mesh NOM, whose info will be received:
         NICF_R = 0
         CALL POINT_TO_MESH(NM)
         ! Loop over cut-cells:
         EXTERNAL_WALL_LOOP_1B : DO IW=1,N_EXTERNAL_WALL_CELLS
            WC=>WALL(IW)
            IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE EXTERNAL_WALL_LOOP_1B
            EWC=>EXTERNAL_WALL(IW); IF (EWC%NOM/=NOM) CYCLE EXTERNAL_WALL_LOOP_1B
            BC=>BOUNDARY_COORD(WC%BC_INDEX)
            IIO = BC%IIG; JJO = BC%JJG; KKO = BC%KKG
            SELECT CASE(BC%IOR)
            CASE( IAXIS); IF(FCVAR(IIO-1,JJO  ,KKO  ,CC_FGSC,IAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            CASE(-IAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,IAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            CASE( JAXIS); IF(FCVAR(IIO  ,JJO-1,KKO  ,CC_FGSC,JAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            CASE(-JAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,JAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            CASE( KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO-1,CC_FGSC,KAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            CASE(-KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,KAXIS) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_1B
            END SELECT
            DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
               DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
                  DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                   IIOF=IIO; JJOF=JJO; KKOF=KKO
                   SELECT CASE(-BC%IOR)
                   CASE( IAXIS); IIOF=IIO-1
                   CASE( JAXIS); JJOF=JJO-1
                   CASE( KAXIS); KKOF=KKO-1
                   END SELECT
                   ICF=MESHES(NOM)%FCVAR(IIOF,JJOF,KKOF,CC_IDCF,ABS(BC%IOR))
                   IF (ICF > 0) THEN
                      NICF_R = NICF_R + 1
                      M3%ICF_UFFB_CF_R(4*NICF_R-3:4*NICF_R) = (/IIOF,JJOF,KKOF,ABS(BC%IOR)/) ! Note:ICF index refers to NOM mesh
                      CF => MESHES(NOM)%CUT_FACE(ICF)
                      IF(ALLOCATED(CF%VELS_OMESH)) DEALLOCATE(CF%VELS_OMESH)
                      IF(ALLOCATED(CF% VEL_OMESH)) DEALLOCATE(CF% VEL_OMESH)
                      IF(ALLOCATED(CF% VEL_LNK)) DEALLOCATE(CF% VEL_LNK) ! For special cases of grid refinement.
                      IF(ALLOCATED(CF% VEL_LNK_OMESH)) DEALLOCATE(CF% VEL_LNK_OMESH)
                      IF(ALLOCATED(CF%  FN_OMESH)) DEALLOCATE(CF%  FN_OMESH)
                      ALLOCATE(CF%VELS_OMESH(1:CF%NFACE));     CF%   VELS_OMESH = 0._EB
                      ALLOCATE(CF% VEL_OMESH(1:CF%NFACE));     CF%    VEL_OMESH = 0._EB
                      ALLOCATE(CF% VEL_LNK(1:CF%NFACE));       CF%      VEL_LNK = 0._EB
                      ALLOCATE(CF% VEL_LNK_OMESH(1:CF%NFACE)); CF%VEL_LNK_OMESH = 0._EB
                      ALLOCATE(CF%  FN_OMESH(1:CF%NFACE));     CF%     FN_OMESH = 0._EB
                   ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO EXTERNAL_WALL_LOOP_1B
      ENDIF

      ! Here allocate in M3 U_LNK,V_LNK,W_LNK if necessary:
      ! Here test should be done to figure out if any of the faces exchanged is actually a linked face.
      ! For now assume there are:
      IF (MESHES(NM)%OMESH(NOM)%NIC_R > 0) THEN ! There are velocity variables to receive from NOM.
         M3 => MESHES(NM)%OMESH(NOM)
         M3%NLKF_R = (M3%I_MAX_R-M3%I_MIN_R+1)*(M3%J_MAX_R-M3%J_MIN_R+1)*(M3%K_MAX_R-M3%K_MIN_R+1)
         IF(ALLOCATED(M3%U_LNK)) DEALLOCATE(M3%U_LNK)
         IF(ALLOCATED(M3%V_LNK)) DEALLOCATE(M3%V_LNK)
         IF(ALLOCATED(M3%W_LNK)) DEALLOCATE(M3%W_LNK)
         ALLOCATE(M3%U_LNK(M3%I_MIN_R:M3%I_MAX_R,M3%J_MIN_R:M3%J_MAX_R,M3%K_MIN_R:M3%K_MAX_R))
         ALLOCATE(M3%V_LNK(M3%I_MIN_R:M3%I_MAX_R,M3%J_MIN_R:M3%J_MAX_R,M3%K_MIN_R:M3%K_MAX_R))
         ALLOCATE(M3%W_LNK(M3%I_MIN_R:M3%I_MAX_R,M3%J_MIN_R:M3%J_MAX_R,M3%K_MIN_R:M3%K_MAX_R))
         M3%U_LNK = 0._EB; M3%V_LNK = 0._EB; M3%W_LNK = 0._EB
      ENDIF
      IF (MESHES(NM)%OMESH(NOM)%NIC_S>0) THEN
         M3 => MESHES(NM)%OMESH(NOM)
         M3%NLKF_S = (M3%I_MAX_S-M3%I_MIN_S+1)*(M3%J_MAX_S-M3%J_MIN_S+1)*(M3%K_MAX_S-M3%K_MIN_S+1)
      ENDIF
   ENDDO
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

IF (N_MPI_PROCESSES>1) ALLOCATE(REQ0(NMESHES))
! Exchange number of cut-cells information to be exchanged between MESH and OMESHES:
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. MESHES(NOM)%CONNECTED_MESH(NM)) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%NICC_S(1),2,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
! DEFINITION NCC_S:   MESHES(NOM)%OMESH(NM)%NCC_S   = MESHES(NM)%OMESH(NOM)%NCC_R
DO NM=1,NMESHES
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)  ! This call orders the sending mesh by mesh.
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   DO NOM=1,NMESHES
      IF (NM/=NOM .AND. .NOT.MESHES(NM)%CONNECTED_MESH(NOM)) CYCLE
      M3 => MESHES(NM)%OMESH(NOM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NOM)/=MY_RANK .AND. MESHES(NM)%CONNECTED_MESH(NOM)) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%NICC_R(1),2,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         ! 2D, NM/NOM, CONNECTED_MESH=F and Several MPI process run.
         IF(.NOT. ALLOCATED(MESHES(NOM)%OMESH)) CYCLE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%NICC_S(1:2) = M3%NICC_R(1:2)
      ENDIF
   ENDDO
ENDDO
IF ((N_REQ0>0) .AND. (N_MPI_PROCESSES>1)) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)
! At this point values of M2%NICC_S should have been received.
! Definition: MESHES(NOM)%OMESH(NM)%UNKZ_CT_S(:) = MESHES(NM)%OMESH(NOM)%UNKZ_CT_R(:)
!             MESHES(NOM)%OMESH(NM)%UNKZ_CC_S(:) = MESHES(NM)%OMESH(NOM)%UNKZ_CC_R(:)
! Now allocate buffers to send UNKZ information:
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (MESHES(NOM)%OMESH(NM)%NIC_S>0) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         ALLOCATE(M2%UNKZ_CT_S(M2%NIC_S))
      ENDIF

      IF (MESHES(NOM)%OMESH(NM)%NICC_S(1)>0) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         ALLOCATE(M2%ICC_UNKZ_CC_S(3*M2%NICC_S(1))); M2%ICC_UNKZ_CC_S = CC_UNDEFINED
         ALLOCATE(M2%UNKZ_CC_S(M2%NICC_S(2)))
      ENDIF
   ENDDO
ENDDO

! Exchange list of cutcells in ICC_UNKZ_CC_S/R:
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NICC_S(1)>0) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%ICC_UNKZ_CC_S(1),3*M2%NICC_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICC_R(1)<1) CYCLE
      IF (PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%ICC_UNKZ_CC_R(1),3*M3%NICC_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%ICC_UNKZ_CC_S(1:3*M2%NICC_S(1)) = M3%ICC_UNKZ_CC_R(1:3*M3%NICC_R(1))
         ! Here write the ICC into ICC_UNKZ_CC_S for OMESH NM:
         M => MESHES(NOM); DO ICC1=1,M2%NICC_S(1)
            IIO=M2%ICC_UNKZ_CC_S(3*ICC1-2)
            JJO=M2%ICC_UNKZ_CC_S(3*ICC1-1)
            KKO=M2%ICC_UNKZ_CC_S(3*ICC1  )
            M2%ICC_UNKZ_CC_S(ICC1) = M%CCVAR(IIO,JJO,KKO,CC_IDCC)
         ENDDO
         ALLOCATE(INT1D(1:M2%NICC_S(1))); INT1D(1:M2%NICC_S(1)) = M2%ICC_UNKZ_CC_S(1:M2%NICC_S(1))
         CALL MOVE_ALLOC(FROM=INT1D,TO=M2%ICC_UNKZ_CC_S)
      ENDIF
   ENDDO
ENDDO
IF ((N_REQ0>0) .AND. (N_MPI_PROCESSES>1)) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NICC_S(1)>0) THEN
         ! Here write the ICC into ICC_UNKZ_CC_S for OMESH NM:
         M => MESHES(NOM); DO ICC1=1,M2%NICC_S(1)
            IIO=M2%ICC_UNKZ_CC_S(3*ICC1-2)
            JJO=M2%ICC_UNKZ_CC_S(3*ICC1-1)
            KKO=M2%ICC_UNKZ_CC_S(3*ICC1  )
            M2%ICC_UNKZ_CC_S(ICC1) = M%CCVAR(IIO,JJO,KKO,CC_IDCC)
         ENDDO
         ALLOCATE(INT1D(1:M2%NICC_S(1))); INT1D(1:M2%NICC_S(1)) = M2%ICC_UNKZ_CC_S(1:M2%NICC_S(1))
         CALL MOVE_ALLOC(FROM=INT1D,TO=M2%ICC_UNKZ_CC_S)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICC_R(1)<1) CYCLE
      M => MESHES(NOM); DO ICC1=1,M3%NICC_R(1)
         IIO=M3%ICC_UNKZ_CC_R(3*ICC1-2)
         JJO=M3%ICC_UNKZ_CC_R(3*ICC1-1)
         KKO=M3%ICC_UNKZ_CC_R(3*ICC1  )
         M3%ICC_UNKZ_CC_R(ICC1) = M%CCVAR(IIO,JJO,KKO,CC_IDCC)
      ENDDO
      ALLOCATE(INT1D(1:M3%NICC_R(1))); INT1D(1:M3%NICC_R(1)) = M3%ICC_UNKZ_CC_R(1:M3%NICC_R(1))
      CALL MOVE_ALLOC(FROM=INT1D,TO=M3%ICC_UNKZ_CC_R)
   ENDDO
ENDDO

! Then exchange cut-face ICF values:
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. MESHES(NOM)%CONNECTED_MESH(NM)) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%NICF_S(1),2,MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
DO NM=1,NMESHES
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)  ! This call orders the sending mesh by mesh.
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   DO NOM=1,NMESHES
      IF (NM/=NOM .AND. .NOT.MESHES(NM)%CONNECTED_MESH(NOM)) CYCLE
      M3 => MESHES(NM)%OMESH(NOM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NOM)/=MY_RANK .AND. MESHES(NM)%CONNECTED_MESH(NOM)) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%NICF_R(1),2,MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         ! 2D, NM/NOM, CONNECTED_MESH=F and Several MPI process run.
         IF(.NOT. ALLOCATED(MESHES(NOM)%OMESH)) CYCLE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%NICF_S(1:2) = M3%NICF_R(1:2)
      ENDIF
   ENDDO
ENDDO
IF ((N_REQ0>0) .AND. (N_MPI_PROCESSES>1)) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

! Now allocate buffers to send ICF information:
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      IF (MESHES(NOM)%OMESH(NM)%NICF_S(1)>0) THEN
         M2 => MESHES(NOM)%OMESH(NM)
         ALLOCATE(M2%ICF_UFFB_CF_S(4*M2%NICF_S(1)))
      ENDIF
   ENDDO
ENDDO
N_REQ0 = 0
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NICF_S(1)>0) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M2%ICF_UFFB_CF_S(1),4*M2%NICF_S(1),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICF_R(1)<1) CYCLE
      IF (PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M3%ICF_UFFB_CF_R(1),4*M3%NICF_R(1),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M2%ICF_UFFB_CF_S(1:4*M2%NICF_S(1)) = M3%ICF_UFFB_CF_R(1:4*M3%NICF_R(1))
         ! Here write the ICF into ICF_UFFB_CF_S for OMESH NM:
         M => MESHES(NOM); DO ICF=1,M2%NICF_S(1)
            IIOF= M2%ICF_UFFB_CF_S(4*ICF-3)
            JJOF= M2%ICF_UFFB_CF_S(4*ICF-2)
            KKOF= M2%ICF_UFFB_CF_S(4*ICF-1)
            X1AXIS=M2%ICF_UFFB_CF_S(4*ICF)
            M2%ICF_UFFB_CF_S(ICF) = M%FCVAR(IIOF,JJOF,KKOF,CC_IDCF,X1AXIS)
         ENDDO
         ALLOCATE(INT1D(1:M2%NICF_S(1))); INT1D(1:M2%NICF_S(1)) = M2%ICF_UFFB_CF_S(1:M2%NICF_S(1))
         CALL MOVE_ALLOC(FROM=INT1D,TO=M2%ICF_UFFB_CF_S)
      ENDIF
   ENDDO
ENDDO
IF ((N_REQ0>0) .AND. (N_MPI_PROCESSES>1)) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NICF_S(1)>0) THEN
         ! Here write the ICF into ICF_UFFB_CF_S for OMESH NM:
         M => MESHES(NOM); DO ICF=1,M2%NICF_S(1)
            IIOF= M2%ICF_UFFB_CF_S(4*ICF-3)
            JJOF= M2%ICF_UFFB_CF_S(4*ICF-2)
            KKOF= M2%ICF_UFFB_CF_S(4*ICF-1)
            X1AXIS=M2%ICF_UFFB_CF_S(4*ICF)
            M2%ICF_UFFB_CF_S(ICF) = M%FCVAR(IIOF,JJOF,KKOF,CC_IDCF,X1AXIS)
         ENDDO
         ALLOCATE(INT1D(1:M2%NICF_S(1))); INT1D(1:M2%NICF_S(1)) = M2%ICF_UFFB_CF_S(1:M2%NICF_S(1))
         CALL MOVE_ALLOC(FROM=INT1D,TO=M2%ICF_UFFB_CF_S)
      ENDIF
   ENDDO
ENDDO
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICF_R(1)<1) CYCLE
      M => MESHES(NOM); DO ICF=1,M3%NICF_R(1)
         IIOF= M3%ICF_UFFB_CF_R(4*ICF-3)
         JJOF= M3%ICF_UFFB_CF_R(4*ICF-2)
         KKOF= M3%ICF_UFFB_CF_R(4*ICF-1)
         X1AXIS=M3%ICF_UFFB_CF_R(4*ICF)
         M3%ICF_UFFB_CF_R(ICF) = M%FCVAR(IIOF,JJOF,KKOF,CC_IDCF,X1AXIS)
      ENDDO
      ALLOCATE(INT1D(1:M3%NICF_R(1))); INT1D(1:M3%NICF_R(1)) = M3%ICF_UFFB_CF_R(1:M3%NICF_R(1))
      CALL MOVE_ALLOC(FROM=INT1D,TO=M3%ICF_UFFB_CF_R)
   ENDDO
ENDDO

! Allocate VEL_LNK on Sending cut-faces if not yet allocated:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICF_S(1)>0) THEN
         DO IW=1,M3%NICF_S(1)
            ICF=  M3%ICF_UFFB_CF_S(IW)
            IF (ICF<1) CYCLE
            CF => MESHES(NM)%CUT_FACE(ICF)
            IF(.NOT.ALLOCATED(MESHES(NM)%CUT_FACE(ICF)%VEL_LNK)) ALLOCATE(MESHES(NM)%CUT_FACE(ICF)%VEL_LNK(1:CF%NFACE))
            MESHES(NM)%CUT_FACE(ICF)%VEL_LNK = 0._EB
         ENDDO
      ENDIF
   ENDDO
ENDDO

! Senders populate UNKZ_CC_S with computed UNKZ values:
ALLOCATE(NCC_SV(1:NMESHES));
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (M2%NICC_S(1)<1) CYCLE
      M => MESHES(NOM)
      NCC_SV(NOM) = 0
      DO ICC1=1,M2%NICC_S(1)
         ICC = M2%ICC_UNKZ_CC_S(ICC1)
         NCELL=M%CUT_CELL(ICC)%NCELL
         DO JCC=1,NCELL
            NCC_SV(NOM) = NCC_SV(NOM) + 1
            M2%UNKZ_CC_S(NCC_SV(NOM)) = M%CUT_CELL(ICC)%UNKZ(JCC)
         ENDDO
      ENDDO
   ENDDO
ENDDO
DEALLOCATE(NCC_SV)

! Finally exchange UNKZ values:
N_REQ0 = 0
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   DO NOM=1,NMESHES
      M3 => MESHES(NM)%OMESH(NOM)
      IF (M3%NICC_R(1)<1) CYCLE
      IF (PROCESS(NOM)/=MY_RANK) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_IRECV(M3%UNKZ_CC_R(1),M3%NICC_R(2),MPI_INTEGER,PROCESS(NOM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ELSE
         M2 => MESHES(NOM)%OMESH(NM)
         M3%UNKZ_CC_R(1:M3%NICC_R(2)) = M2%UNKZ_CC_S(1:M2%NICC_S(2))
      ENDIF
   ENDDO
ENDDO
DO NM=1,NMESHES
   DO NOM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      M2 => MESHES(NOM)%OMESH(NM)
      IF (N_MPI_PROCESSES>1 .AND. NM/=NOM .AND. PROCESS(NM)/=MY_RANK .AND. M2%NICC_S(1)>0) THEN
         N_REQ0 = N_REQ0 + 1; CALL CHECK_REQ0_SIZE
         CALL MPI_ISEND(M2%UNKZ_CC_S(1),M2%NICC_S(2),MPI_INTEGER,PROCESS(NM),NM,MPI_COMM_WORLD,REQ0(N_REQ0),IERR)
      ENDIF
   ENDDO
ENDDO
IF ( (N_REQ0>0) .AND. (N_MPI_PROCESSES>1) ) CALL MPI_WAITALL(N_REQ0,REQ0(1:N_REQ0),MPI_STATUSES_IGNORE,IERR)

! Copy to guard-cell cut-cells:
ALLOCATE(NCC_SV(1:NMESHES));
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   M => MESHES(NM)
   NCC_SV(:)=0
   CALL POINT_TO_MESH(NM)
   ! Loop over cut-cells:
   EXTERNAL_WALL_LOOP_2 : DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      EWC=>EXTERNAL_WALL(IW)
      BC=>BOUNDARY_COORD(WC%BC_INDEX)
      IF (.NOT.(WC%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY .OR. &
                WC%BOUNDARY_TYPE == MIRROR_BOUNDARY) ) CYCLE EXTERNAL_WALL_LOOP_2
      IF ( WC%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY ) THEN
         IF (CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_2
         NOM = EWC%NOM
         DO KKO=EWC%KKO_MIN,EWC%KKO_MAX
            DO JJO=EWC%JJO_MIN,EWC%JJO_MAX
               DO IIO=EWC%IIO_MIN,EWC%IIO_MAX
                 ICC   = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
                 IF (ICC > 0) THEN
                    DO JCC=1,MESHES(NOM)%CUT_CELL(ICC)%NCELL
                       NCC_SV(NOM)=NCC_SV(NOM)+1
                       MESHES(NOM)%CUT_CELL(ICC)%UNKZ(JCC) = M%OMESH(NOM)%UNKZ_CC_R(NCC_SV(NOM))
                    ENDDO
                 ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSEIF ( WC%BOUNDARY_TYPE==MIRROR_BOUNDARY ) THEN
         IIO = BC%IIG; JJO = BC%JJG; KKO = BC%KKG
         IF (CCVAR(BC%II,BC%JJ,BC%KK,CC_CGSC) /= CC_CUTCFE) CYCLE EXTERNAL_WALL_LOOP_2
         ! CYCLE if OBJECT face is in the Mirror Boundary, normal out into ghost-cell:
         SELECT CASE(BC%IOR)
         CASE( IAXIS); IF(FCVAR(IIO-1,JJO  ,KKO  ,CC_FGSC,IAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         CASE(-IAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,IAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         CASE( JAXIS); IF(FCVAR(IIO  ,JJO-1,KKO  ,CC_FGSC,JAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         CASE(-JAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,JAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         CASE( KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO-1,CC_FGSC,KAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         CASE(-KAXIS); IF(FCVAR(IIO  ,JJO  ,KKO  ,CC_FGSC,KAXIS) == CC_SOLID) CYCLE EXTERNAL_WALL_LOOP_2
         END SELECT
         NOM = NM; ICC = MESHES(NOM)%CCVAR(IIO,JJO,KKO,CC_IDCC)
         IF (ICC > 0) THEN
            DO JCC=1,MESHES(NOM)%CUT_CELL(ICC)%NCELL
               NCC_SV(NOM)=NCC_SV(NOM)+1
               MESHES(NOM)%CUT_CELL(ICC)%UNKZ(JCC) = M%OMESH(NOM)%UNKZ_CC_R(NCC_SV(NOM))
            ENDDO
         ENDIF
      ENDIF
   ENDDO EXTERNAL_WALL_LOOP_2
ENDDO
DEALLOCATE(NCC_SV)

! Finally Exchange Cartesian cell UNKZ:
IF (N_MPI_PROCESSES>1) THEN
   DO NM=1,NMESHES
      IF (MPI_COMM_NEIGHBORS(NM)==MPI_COMM_NULL) CYCLE
      M => MESHES(NM)
      ! X direction bounds:
      ILO_FACE = 0                    ! Low mesh boundary face index.
      IHI_FACE = M%IBAR               ! High mesh boundary face index.
      ISTR     = ILO_FACE - NGUARD    ! Allocation start x arrays.
      IEND     = IHI_FACE + NGUARD    ! Allocation end x arrays.

      ! Y direction bounds:
      JLO_FACE = 0                    ! Low mesh boundary face index.
      JHI_FACE = M%JBAR               ! High mesh boundary face index.
      JSTR     = JLO_FACE - NGUARD    ! Allocation start y arrays.
      JEND     = JHI_FACE + NGUARD    ! Allocation end y arrays.

      ! Z direction bounds:
      KLO_FACE = 0                    ! Low mesh boundary face index.
      KHI_FACE = M%KBAR               ! High mesh boundary face index.
      KSTR     = KLO_FACE - NGUARD    ! Allocation start z arrays.
      KEND     = KHI_FACE + NGUARD    ! Allocation end z arrays.

      ALL_FLG = .FALSE.
      IF (.NOT.ALLOCATED(M%CCVAR)) THEN; ALL_FLG=.TRUE.; ALLOCATE(M%CCVAR(ISTR:IEND,JSTR:JEND,KSTR:KEND,CC_UNKZ:CC_UNKZ)); ENDIF
      N_INT = (IEND-ISTR+1)*(JEND-JSTR+1)*(KEND-KSTR+1)
      CALL MPI_BCAST(M%CCVAR(ISTR,JSTR,KSTR,CC_UNKZ),N_INT,MPI_INTEGER,MPI_COMM_NEIGHBORS_ROOT(NM),MPI_COMM_NEIGHBORS(NM),IERR)
      IF (ALL_FLG) DEALLOCATE(M%CCVAR)
   ENDDO
ENDIF

IF (N_MPI_PROCESSES>1) THEN
   DEALLOCATE(REQ0)
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
ENDIF

RETURN

CONTAINS
SUBROUTINE CHECK_REQ0_SIZE
IF(N_REQ0>SIZE(REQ0,DIM=1)) THEN
   ALLOCATE(REQ0DUM(SIZE(REQ0,DIM=1)+NMESHES))
   REQ0DUM(1:N_REQ0-1) = REQ0(1:N_REQ0-1)
   CALL MOVE_ALLOC(REQ0DUM,REQ0)
ENDIF
END SUBROUTINE CHECK_REQ0_SIZE

END SUBROUTINE FILL_UNKZ_GUARDCELLS



! --------------------------- GET_MMATRIX_SCALAR_3D -------------------------------

SUBROUTINE GET_MMATRIX_SCALAR_3D


! Local Variables:
INTEGER :: NM
INTEGER :: I,J,K,IROW,IROW_LOC,ICC,ICC2

! INTEGER :: ILK

! Allocate mass matrix: Diagonal containing cell volumes on implicit region:
ALLOCATE(  M_MAT_Z(1:NUNKZ_LOCAL) );  M_MAT_Z = 0._EB
ALLOCATE( JM_MAT_Z(1:NUNKZ_LOCAL) ); JM_MAT_Z = 0 ! local index of diagonal entry in JD_MAT_Z

! Mesh Loop:
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! 1. Number Regular GASPHASE cells:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_UNKZ) <= 0) CYCLE ! Either explicit region or solid cell.
            IROW     = CCVAR(I,J,K,CC_UNKZ)
            IROW_LOC = IROW - UNKZ_IND(NM_START)
            M_MAT_Z(IROW_LOC) = M_MAT_Z(IROW_LOC) + DX(I)*DY(J)*DZ(K)
         ENDDO
      ENDDO
   ENDDO

   ! 2. Now Cut cells:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      ! Drop cut-cells inside an OBST:
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO ICC2 = 1,CC%NCELL
         IROW     = CC%UNKZ(ICC2)
         IROW_LOC = IROW - UNKZ_IND(NM_START)
         M_MAT_Z(IROW_LOC) = M_MAT_Z(IROW_LOC) + CC%VOLUME(ICC2)
      ENDDO
   ENDDO

ENDDO MESH_LOOP

RETURN
END SUBROUTINE GET_MMATRIX_SCALAR_3D


! ---------------------- GET_GASPHASE_CUTFACES_DATA -----------------------------

SUBROUTINE GET_GASPHASE_CUTFACES_DATA

USE MPI_F08

! Local variables:
INTEGER :: NM
INTEGER :: NCELL,ICC,JCC,IFC,IFACE,LOWHIGH,ICF1,ICF2
INTEGER :: IW,II,JJ,KK,IIF,JJF,KKF,IOR,IIG,JJG,KKG,LOWHIGH_TEST,LOWHIGH_TEST_G,X1AXIS
INTEGER :: IERR

! Mesh loop:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   CALL POINT_TO_MESH(NM)

   ! First Scalars:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      NCELL = CUT_CELL(ICC)%NCELL
      DO JCC=1,NCELL
         DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
            IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
            ! If face type in face_list is not CC_FTYPE_CFGAS, drop:
            IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_CFGAS ) CYCLE
            ! Which face?
            LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
            ICF1    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
            ICF2    = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
            SELECT CASE(LOWHIGH)
            CASE( LOW_IND); CUT_FACE(ICF1)%UNKZ(HIGH_IND,ICF2) = CUT_CELL(ICC)%UNKZ(JCC) !Cutface on low side of cutcell
            CASE(HIGH_IND); CUT_FACE(ICF1)%UNKZ( LOW_IND,ICF2) = CUT_CELL(ICC)%UNKZ(JCC) !CF on high side of CC.
            END SELECT
         ENDDO
      ENDDO
   ENDDO

   ! Now Apply wall cell loop for cut cells:
   WALL_CUTFACE_LOOP :  DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      II  = BC%II
      JJ  = BC%JJ
      KK  = BC%KK
      IOR = BC%IOR

      ! Drop if face is not of type CC_CUTCFE:
      X1AXIS=ABS(IOR)
      SELECT CASE(IOR)                                           ! V Face on high side of Guard-Cell
      CASE( IAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST = HIGH_IND; LOWHIGH_TEST_G =  LOW_IND
      CASE(-IAXIS); IIF=II-1; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST =  LOW_IND; LOWHIGH_TEST_G = HIGH_IND
      CASE( JAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST = HIGH_IND; LOWHIGH_TEST_G =  LOW_IND
      CASE(-JAXIS); IIF=II  ; JJF=JJ-1; KKF=KK  ; LOWHIGH_TEST =  LOW_IND; LOWHIGH_TEST_G = HIGH_IND
      CASE( KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST = HIGH_IND; LOWHIGH_TEST_G =  LOW_IND
      CASE(-KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK-1; LOWHIGH_TEST =  LOW_IND; LOWHIGH_TEST_G = HIGH_IND
      END SELECT

      IF (FCVAR(IIF,JJF,KKF,CC_FGSC,X1AXIS) /= CC_CUTCFE) CYCLE WALL_CUTFACE_LOOP

      BC => BOUNDARY_COORD(WC%BC_INDEX)
      IIG  = BC%IIG
      JJG  = BC%JJG
      KKG  = BC%KKG

      ! Add IWC field to CUT_FACE from internal cell:
      ICC = MESHES(NM)%CCVAR(IIG,JJG,KKG,CC_IDCC)
      DO JCC=1,CUT_CELL(ICC)%NCELL
         ! Loop faces and test:
         DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
            IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
            ! Which face ?
            LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
            IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_CFGAS) CYCLE ! Must Be gasphase cut-face
            IF ( LOWHIGH                          /=  LOWHIGH_TEST_G) CYCLE ! In same side as EWC from internal cell
            IF ( CUT_CELL(ICC)%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE ! Normal to same axis as EWC
            ICF1    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
            CUT_FACE(ICF1)%IWC = IW   ! Rest of info from internal cut-cell has been filled in previous ICC loop.
            WC%CUT_FACE_INDEX  = ICF1 ! Add CUT_FACE index in WALL(:) array.
         ENDDO
      ENDDO

      ! Now CCVAR(II,JJ,KK,CC_CGSC) from guard cell:
      IF (CELL(CELL_INDEX(II,JJ,KK))%SOLID) CYCLE WALL_CUTFACE_LOOP
      ICC = MESHES(NM)%CCVAR(II,JJ,KK,CC_IDCC)
      DO JCC=1,CUT_CELL(ICC)%NCELL
         ! Loop faces and test:
         DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
            IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
            ! Which face ?
            LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
            IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_CFGAS) CYCLE ! Must Be gasphase cut-face
            IF ( LOWHIGH                              /= LOWHIGH_TEST) CYCLE ! In same side as EWC from guard-cell
            IF ( CUT_CELL(ICC)%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE ! Normal to same axis as EWC

            ICF1    = CUT_CELL(ICC)%FACE_LIST(4,IFACE)
            ICF2    = CUT_CELL(ICC)%FACE_LIST(5,IFACE)
            SELECT CASE(LOWHIGH)
            CASE( LOW_IND); CUT_FACE(ICF1)%UNKZ(HIGH_IND,ICF2) = CUT_CELL(ICC)%UNKZ(JCC) !Cutface on low side of cutcell
            CASE(HIGH_IND); CUT_FACE(ICF1)%UNKZ( LOW_IND,ICF2) = CUT_CELL(ICC)%UNKZ(JCC)
            END SELECT
         ENDDO
      ENDDO

   ENDDO WALL_CUTFACE_LOOP

ENDDO MAIN_MESH_LOOP

IF (DEBUG_MATVEC_DATA) THEN
   DBG_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      IF(MY_RANK==PROCESS(NM)) THEN
      CALL POINT_TO_MESH(NM)
      WRITE(LU_ERR,*) ' '
      WRITE(LU_ERR,*) 'MY_RANK, NM, N_BBCUTFACE_MESH : ',MY_RANK,NM,MESHES(NM)%N_BBCUTFACE_MESH
      DO IFC=1,MESHES(NM)%N_BBCUTFACE_MESH
         IF(CUT_FACE(IFC)%STATUS/=CC_GASPHASE) CYCLE
         WRITE(LU_ERR,*) 'BB CUT_FACE, IFC, IWC=',IFC,CUT_FACE(IFC)%IWC,CUT_FACE(IFC)%STATUS
      ENDDO
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
   ENDDO DBG_MESH_LOOP
ENDIF

RETURN
END SUBROUTINE GET_GASPHASE_CUTFACES_DATA



! ---------------------- GET_GASPHASE_REGRCFACES_DATA -----------------------------

SUBROUTINE GET_GASPHASE_REGRCFACES_DATA

USE MPI_F08

! Local variables:
INTEGER :: NM
INTEGER :: ILO,IHI,JLO,JHI,KLO,KHI
INTEGER :: I,J,K,II,IREG,IRC,IIFC,X1AXIS,X2AXIS,X3AXIS
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJKBUFFER
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOHIBUFF
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: IJKFACE
INTEGER :: ICC,JCC,IJK(MAX_DIM),IFC,IFACE,LOWHIGH
INTEGER :: XIAXIS,XJAXIS,XKAXIS,INDXI1(MAX_DIM),INCELL,JNCELL,KNCELL,INFACE,JNFACE,KNFACE
INTEGER :: ISTR, IEND, JSTR, JEND, KSTR, KEND
LOGICAL :: INLIST

INTEGER :: IW,JJ,KK,IIF,JJF,KKF,IOR,LOWHIGH_TEST,IIG,JJG,KKG
INTEGER :: IBNDINT,IC,IC2
LOGICAL :: FLGIN
INTEGER, PARAMETER :: OZPOS=0, ICPOS=1, JCPOS=2, IFPOS=3
INTEGER :: IERR
TYPE(MESH_TYPE), POINTER :: M

! Mesh loop:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M => MESHES(NM)
   CALL POINT_TO_MESH(NM)

   ! Mesh sizes:
   NXB=IBAR
   NYB=JBAR
   NZB=KBAR

   ! X direction bounds:
   ILO_FACE = 0                    ! Low mesh boundary face index.
   IHI_FACE = IBAR                 ! High mesh boundary face index.
   ILO_CELL = ILO_FACE + 1     ! First internal cell index. See notes.
   IHI_CELL = IHI_FACE ! Last internal cell index.
   ISTR     = ILO_FACE - NGUARD    ! Allocation start x arrays.
   IEND     = IHI_FACE + NGUARD    ! Allocation end x arrays.


   ! Y direction bounds:
   JLO_FACE = 0                    ! Low mesh boundary face index.
   JHI_FACE = JBAR                 ! High mesh boundary face index.
   JLO_CELL = JLO_FACE + 1     ! First internal cell index. See notes.
   JHI_CELL = JHI_FACE ! Last internal cell index.
   JSTR     = JLO_FACE - NGUARD    ! Allocation start y arrays.
   JEND     = JHI_FACE + NGUARD    ! Allocation end y arrays.

   ! Z direction bounds:
   KLO_FACE = 0                    ! Low mesh boundary face index.
   KHI_FACE = KBAR                 ! High mesh boundary face index.
   KLO_CELL = KLO_FACE + 1     ! First internal cell index. See notes.
   KHI_CELL = KHI_FACE ! Last internal cell index.
   KSTR     = KLO_FACE - NGUARD    ! Allocation start z arrays.
   KEND     = KHI_FACE + NGUARD    ! Allocation end z arrays.

   ! Define grid arrays for this mesh:
   ! Populate position and cell size arrays: Uniform grid implementation.
   ! X direction:
   ALLOCATE(XCELL(ISTR:IEND));  XCELL = 1._EB/GEOMEPS ! Initialize huge.
   XCELL(ILO_CELL-1:IHI_CELL+1) = M%XC(ILO_CELL-1:IHI_CELL+1)

   ! Y direction:
   ALLOCATE(YCELL(JSTR:JEND));  YCELL = 1._EB/GEOMEPS ! Initialize huge.
   YCELL(JLO_CELL-1:JHI_CELL+1) = M%YC(JLO_CELL-1:JHI_CELL+1)

   ! Z direction:
   ALLOCATE(ZCELL(KSTR:KEND));  ZCELL = 1._EB/GEOMEPS ! Initialize huge.
   ZCELL(KLO_CELL-1:KHI_CELL+1) = M%ZC(KLO_CELL-1:KHI_CELL+1)

   ! Set starting number of regular faces for NM to zero:
   M%CC_NREGFACE_Z(IAXIS:KAXIS) = 0

   ! 1. Regular GASPHASE faces connected to Gasphase cells:
   ALLOCATE(IJKBUFFER(IAXIS:KAXIS+1,1:(NXB+1)*(NYB+1)*(NZB+1)))
   ALLOCATE(LOHIBUFF(LOW_IND:HIGH_IND,1:(NXB+1)*(NYB+1)*(NZB+1)))

   ! First Scalars:
   ! axis = IAXIS:
   X1AXIS = IAXIS
   IJKBUFFER(:,:)=0; LOHIBUFF(:,:)=.FALSE.; IREG = 0
   ! First Reg Faces in mesh block boundaries, then inside mesh. Count for allocation:
   IBNDINT_LOOP_X : DO IBNDINT=1,3
      IF(IBNDINT==3) M%CC_NBBREGFACE_Z(X1AXIS)=IREG ! REG faces in block boundaries.
      SELECT CASE(IBNDINT)
      CASE(1)
         ILO = ILO_FACE; IHI = ILO_FACE
         JLO = JLO_CELL; JHI = JHI_CELL
         KLO = KLO_CELL; KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I  ,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I+1,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I+1,J,K,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I  ,J,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I+1,J,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I+1  ,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX(-X1AXIS) /) ! Face on low side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         ILO = IHI_FACE; IHI = IHI_FACE
         JLO = JLO_CELL; JHI = JHI_CELL
         KLO = KLO_CELL; KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I  ,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I+1,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I  ,J,K,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I  ,J,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I+1,J,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /) ! Face on high side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         ILO = ILO_FACE+1; IHI = IHI_FACE-1
         JLO = JLO_CELL;   JHI = JHI_CELL
         KLO = KLO_CELL;   KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I  ,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I+1,J,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( (CCVAR(I,J,K,CC_UNKZ)<=0) .AND. (CCVAR(I+1,J,K,CC_UNKZ)<=0) ) CYCLE
                  IREG = IREG + 1
                  IF ( CCVAR(I  ,J,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I+1,J,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /)
                  IC2 =  CELL_INDEX(I+1,J,K)
                  IF(CELL(IC)%WALL_INDEX( X1AXIS)>0 .AND. CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IREG = IREG + 1
                     LOHIBUFF(LOW_IND,IREG) = LOHIBUFF(LOW_IND,IREG-1)
                     LOHIBUFF(HIGH_IND,IREG)= LOHIBUFF(HIGH_IND,IREG-1)
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ELSEIF(CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      END SELECT
   ENDDO IBNDINT_LOOP_X
   M%CC_NREGFACE_Z(X1AXIS) = IREG
   IF(ALLOCATED(M%CC_REGFACE_IAXIS_Z)) DEALLOCATE(M%CC_REGFACE_IAXIS_Z)
   ALLOCATE(M%CC_REGFACE_IAXIS_Z(IREG))
   DO II=1,IREG
      M%CC_REGFACE_IAXIS_Z(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
      M%CC_REGFACE_IAXIS_Z(II)%IWC              = IJKBUFFER(KAXIS+1,II)
      M%CC_REGFACE_IAXIS_Z(II)%DO_LO_IND        = LOHIBUFF(LOW_IND,II)
      M%CC_REGFACE_IAXIS_Z(II)%DO_HI_IND        = LOHIBUFF(HIGH_IND,II)
      ALLOCATE(M%CC_REGFACE_IAXIS_Z(II)%RHOZZ_U(1:N_TOTAL_SCALARS)   ,&
               M%CC_REGFACE_IAXIS_Z(II)%FN_ZZ(1:N_TOTAL_SCALARS)     ,&
               M%CC_REGFACE_IAXIS_Z(II)%RHO_D_DZDN(1:N_TOTAL_SCALARS),&
               M%CC_REGFACE_IAXIS_Z(II)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))
      M%CC_REGFACE_IAXIS_Z(II)%RHOZZ_U      = 0._EB
      M%CC_REGFACE_IAXIS_Z(II)%FN_ZZ        = 0._EB
      M%CC_REGFACE_IAXIS_Z(II)%RHO_D_DZDN   = 0._EB
      M%CC_REGFACE_IAXIS_Z(II)%H_RHO_D_DZDN = 0._EB
   ENDDO

   ! axis = JAXIS:
   X1AXIS = JAXIS
   IJKBUFFER(:,:)=0; LOHIBUFF(:,:)=.FALSE.; IREG = 0
   ! First Reg Faces in mesh block boundaries, then inside mesh. Count for allocation:
   IBNDINT_LOOP_Y : DO IBNDINT=1,3
      IF(IBNDINT==3) M%CC_NBBREGFACE_Z(X1AXIS)=IREG ! REG faces in block boundaries.
      SELECT CASE(IBNDINT)
      CASE(1)
         ILO = ILO_CELL; IHI = IHI_CELL
         JLO = JLO_FACE; JHI = JLO_FACE
         KLO = KLO_CELL; KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J  ,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J+1,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J+1,K,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I,J  ,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J+1,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J+1  ,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX(-X1AXIS) /) ! Face on low side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         ILO = ILO_CELL; IHI = IHI_CELL
         JLO = JHI_FACE; JHI = JHI_FACE
         KLO = KLO_CELL; KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J  ,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J+1,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J  ,K,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I,J  ,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J+1,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /) ! Face on high side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         ILO = ILO_CELL;   IHI = IHI_CELL
         JLO = JLO_FACE+1; JHI = JHI_FACE-1
         KLO = KLO_CELL;   KHI = KHI_CELL
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J  ,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J+1,K,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( (CCVAR(I,J,K,CC_UNKZ)<=0) .AND. (CCVAR(I,J+1,K,CC_UNKZ)<=0) ) CYCLE
                  IREG = IREG + 1
                  IF ( CCVAR(I,J  ,K,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J+1,K,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /)
                  IC2 =  CELL_INDEX(I,J+1,K)
                  IF(CELL(IC)%WALL_INDEX( X1AXIS)>0 .AND. CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IREG = IREG + 1
                     LOHIBUFF(LOW_IND,IREG) = LOHIBUFF(LOW_IND,IREG-1)
                     LOHIBUFF(HIGH_IND,IREG)= LOHIBUFF(HIGH_IND,IREG-1)
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ELSEIF(CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      END SELECT
   ENDDO IBNDINT_LOOP_Y
   M%CC_NREGFACE_Z(X1AXIS) = IREG
   IF(ALLOCATED(M%CC_REGFACE_JAXIS_Z)) DEALLOCATE(M%CC_REGFACE_JAXIS_Z)
   ALLOCATE(M%CC_REGFACE_JAXIS_Z(IREG))
   DO II=1,IREG
      M%CC_REGFACE_JAXIS_Z(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
      M%CC_REGFACE_JAXIS_Z(II)%IWC              = IJKBUFFER(KAXIS+1,II)
      M%CC_REGFACE_JAXIS_Z(II)%DO_LO_IND        = LOHIBUFF(LOW_IND,II)
      M%CC_REGFACE_JAXIS_Z(II)%DO_HI_IND        = LOHIBUFF(HIGH_IND,II)
      ALLOCATE(M%CC_REGFACE_JAXIS_Z(II)%RHOZZ_U(1:N_TOTAL_SCALARS)   ,&
               M%CC_REGFACE_JAXIS_Z(II)%FN_ZZ(1:N_TOTAL_SCALARS)     ,&
               M%CC_REGFACE_JAXIS_Z(II)%RHO_D_DZDN(1:N_TOTAL_SCALARS),&
               M%CC_REGFACE_JAXIS_Z(II)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))
      M%CC_REGFACE_JAXIS_Z(II)%RHOZZ_U      = 0._EB
      M%CC_REGFACE_JAXIS_Z(II)%FN_ZZ        = 0._EB
      M%CC_REGFACE_JAXIS_Z(II)%RHO_D_DZDN   = 0._EB
      M%CC_REGFACE_JAXIS_Z(II)%H_RHO_D_DZDN = 0._EB
   ENDDO

   ! axis = KAXIS:
   X1AXIS = KAXIS
   IJKBUFFER(:,:)=0; LOHIBUFF(:,:)=.FALSE.; IREG = 0
   ! First Reg Faces in mesh block boundaries, then inside mesh.
   IBNDINT_LOOP_Z : DO IBNDINT=1,3
      IF(IBNDINT==3) M%CC_NBBREGFACE_Z(X1AXIS)=IREG ! REG faces in block boundaries.
      SELECT CASE(IBNDINT)
      CASE(1)
         ILO = ILO_CELL; IHI = IHI_CELL
         JLO = JLO_CELL; JHI = JHI_CELL
         KLO = KLO_FACE; KHI = KLO_FACE
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J,K  ,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J,K+1,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J,K+1,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I,J,K  ,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J,K+1,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J,K+1 )
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX(-X1AXIS) /)  ! Face on low side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(2)
         ILO = ILO_CELL; IHI = IHI_CELL
         JLO = JLO_CELL; JHI = JHI_CELL
         KLO = KHI_FACE; KHI = KHI_FACE
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J,K  ,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J,K+1,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J,K  ,CC_UNKZ)<=0 ) CYCLE ! Either face out of CCREGION or EXIM face with CCREGION outside of mesh.
                  IREG = IREG + 1
                  IF ( CCVAR(I,J,K  ,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J,K+1,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /)  ! Face on high side of cell.
               ENDDO
            ENDDO
         ENDDO
      CASE(3)
         ILO = ILO_CELL;   IHI = IHI_CELL
         JLO = JLO_CELL;   JHI = JHI_CELL
         KLO = KLO_FACE+1; KHI = KHI_FACE-1
         DO K=KLO,KHI
            DO J=JLO,JHI
               DO I=ILO,IHI
                  IF ( CCVAR(I,J,K  ,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( CCVAR(I,J,K+1,CC_CGSC) /= CC_GASPHASE ) CYCLE
                  IF ( (CCVAR(I,J,K,CC_UNKZ)<=0) .AND. (CCVAR(I,J,K+1,CC_UNKZ)<=0) ) CYCLE
                  IREG = IREG + 1
                  IF ( CCVAR(I,J,K  ,CC_UNKZ)>0 ) LOHIBUFF(LOW_IND,IREG) = .TRUE.
                  IF ( CCVAR(I,J,K+1,CC_UNKZ)>0 ) LOHIBUFF(HIGH_IND,IREG)= .TRUE.
                  IC   = CELL_INDEX(I,J,K)
                  IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC)%WALL_INDEX( X1AXIS) /)
                  IC2 =  CELL_INDEX(I,J,K+1)
                  IF(CELL(IC)%WALL_INDEX( X1AXIS)>0 .AND. CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IREG = IREG + 1
                     LOHIBUFF(LOW_IND,IREG) = LOHIBUFF(LOW_IND,IREG-1)
                     LOHIBUFF(HIGH_IND,IREG)= LOHIBUFF(HIGH_IND,IREG-1)
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ELSEIF(CELL(IC2)%WALL_INDEX(-X1AXIS)>0) THEN
                     IJKBUFFER(IAXIS:KAXIS+1,IREG) = (/ I, J, K, CELL(IC2)%WALL_INDEX(-X1AXIS) /)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      END SELECT
   ENDDO IBNDINT_LOOP_Z
   M%CC_NREGFACE_Z(X1AXIS) = IREG
   IF(ALLOCATED(M%CC_REGFACE_KAXIS_Z)) DEALLOCATE(M%CC_REGFACE_KAXIS_Z)
   ALLOCATE(M%CC_REGFACE_KAXIS_Z(IREG))
   DO II=1,IREG
      M%CC_REGFACE_KAXIS_Z(II)%IJK(IAXIS:KAXIS) = IJKBUFFER(IAXIS:KAXIS,II)
      M%CC_REGFACE_KAXIS_Z(II)%IWC              = IJKBUFFER(KAXIS+1,II)
      M%CC_REGFACE_KAXIS_Z(II)%DO_LO_IND        = LOHIBUFF(LOW_IND,II)
      M%CC_REGFACE_KAXIS_Z(II)%DO_HI_IND        = LOHIBUFF(HIGH_IND,II)
      ALLOCATE(M%CC_REGFACE_KAXIS_Z(II)%RHOZZ_U(1:N_TOTAL_SCALARS)   ,&
               M%CC_REGFACE_KAXIS_Z(II)%FN_ZZ(1:N_TOTAL_SCALARS)     ,&
               M%CC_REGFACE_KAXIS_Z(II)%RHO_D_DZDN(1:N_TOTAL_SCALARS),&
               M%CC_REGFACE_KAXIS_Z(II)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))
      M%CC_REGFACE_KAXIS_Z(II)%RHOZZ_U      = 0._EB
      M%CC_REGFACE_KAXIS_Z(II)%FN_ZZ        = 0._EB
      M%CC_REGFACE_KAXIS_Z(II)%RHO_D_DZDN   = 0._EB
      M%CC_REGFACE_KAXIS_Z(II)%H_RHO_D_DZDN = 0._EB
   ENDDO

   ! 2. Lists of Regular Gasphase faces, connected to one regular gasphase and one cut-cell:
   ! First count for allocation:
   ALLOCATE( IJKFACE(ILO_FACE:IHI_FACE,JLO_FACE:JHI_FACE,KLO_FACE:KHI_FACE,IAXIS:KAXIS,OZPOS:IFPOS) )
   IJKFACE(:,:,:,:,:) = 0
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IJK(IAXIS:KAXIS) = CC%IJK(IAXIS:KAXIS)
      DO JCC=1,CC%NCELL
         IF(CC%UNKZ(JCC) < 1) CYCLE
         ! Loop faces and test:
         DO IFC=1,CC%CCELEM(1,JCC)
            IFACE = CC%CCELEM(IFC+1,JCC)
            ! If face type in face_list is not CC_FTYPE_RCGAS, drop:
            IF(CC%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE
            ! Which face?
            LOWHIGH = CC%FACE_LIST(2,IFACE)
            X1AXIS  = CC%FACE_LIST(3,IFACE)
            SELECT CASE(X1AXIS)
            CASE(IAXIS)
               X2AXIS = JAXIS; X3AXIS = KAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
            CASE(JAXIS)
               X2AXIS = KAXIS; X3AXIS = IAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
            CASE(KAXIS)
               X2AXIS = IAXIS; X3AXIS = JAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
            END SELECT

            ! Face indexes:
            INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1+(LOWHIGH-1), IJK(X2AXIS), IJK(X3AXIS) /)
            INFACE = INDXI1(XIAXIS)
            JNFACE = INDXI1(XJAXIS)
            KNFACE = INDXI1(XKAXIS)

            IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,OZPOS) = 1
            IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,ICPOS) = ICC
            IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,JCPOS) = JCC
            IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,IFPOS) = IFACE

         ENDDO
      ENDDO
   ENDDO

   ! Check for RC_FACE on the boundary of the domain, where the cut-cell is in the guard-cell region.
   ! Now Apply external wall cell loop for guard-cell cut cells:
   GUARD_CUT_CELL_LOOP_1A :  DO IW=1,N_EXTERNAL_WALL_CELLS
      WC=>WALL(IW)
      BC => BOUNDARY_COORD(WC%BC_INDEX)
      II = BC%II
      JJ = BC%JJ
      KK = BC%KK
      IOR = BC%IOR

      ! Which face:
      X1AXIS=ABS(IOR)
      SELECT CASE(IOR)
      CASE( IAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK
         LOWHIGH_TEST=HIGH_IND ! Face on high side of Guard-Cell
      CASE(-IAXIS)
         IIF=II-1; JJF=JJ  ; KKF=KK
         LOWHIGH_TEST=LOW_IND
      CASE( JAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK
         LOWHIGH_TEST=HIGH_IND
      CASE(-JAXIS)
         IIF=II  ; JJF=JJ-1; KKF=KK
         LOWHIGH_TEST=LOW_IND
      CASE( KAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK
         LOWHIGH_TEST=HIGH_IND
      CASE(-KAXIS)
         IIF=II  ; JJF=JJ  ; KKF=KK-1
         LOWHIGH_TEST=LOW_IND
      END SELECT

      ! Drop if FACE is not type CC_GASPHASE
      IF (FCVAR(IIF,JJF,KKF,CC_FGSC,X1AXIS) /= CC_GASPHASE) CYCLE GUARD_CUT_CELL_LOOP_1A

      IIG  = BC%IIG
      JJG  = BC%JJG
      KKG  = BC%KKG

      ! Is this an actual RCFACE laying on the mesh boundary, where the cut-cell is in the guard-cell region?
      FLGIN = (CCVAR(II,JJ,KK,CC_CGSC)==CC_CUTCFE) ! Note that this will overrride the Cut-cell to cut-cell case in
                                                     ! the block boundary, picked up in previous loop. Thats fine.

      IF(.NOT.FLGIN) CYCLE GUARD_CUT_CELL_LOOP_1A

      ICC=CCVAR(II,JJ,KK,CC_IDCC)
      DO JCC=1,CUT_CELL(ICC)%NCELL
         ! Loop faces and test:
         DO IFC=1,CUT_CELL(ICC)%CCELEM(1,JCC)
            IFACE = CUT_CELL(ICC)%CCELEM(IFC+1,JCC)
            ! Which face ?
            LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)
            IF ( CUT_CELL(ICC)%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE ! Must Be gasphase cut-face
            IF ( LOWHIGH                             /= LOWHIGH_TEST) CYCLE ! In same side as EWC from guard-cell
            IF ( CUT_CELL(ICC)%FACE_LIST(3,IFACE) /= X1AXIS) CYCLE ! Normal to same axis as EWC

            IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) = 1
            IJKFACE(IIF,JJF,KKF,X1AXIS,ICPOS) = ICC
            IJKFACE(IIF,JJF,KKF,X1AXIS,JCPOS) = JCC
            IJKFACE(IIF,JJF,KKF,X1AXIS,IFPOS) = IFACE

            CYCLE GUARD_CUT_CELL_LOOP_1A
         ENDDO
      ENDDO

   ENDDO GUARD_CUT_CELL_LOOP_1A

   IRC = SUM(IJKFACE(:,:,:,:,OZPOS))
   IF (IRC == 0) THEN
      DEALLOCATE(XCELL,YCELL,ZCELL)
      DEALLOCATE(IJKBUFFER,LOHIBUFF,IJKFACE)
      CYCLE
   ENDIF

   ! Now actual computation for Scalars:
   ALLOCATE( M%RC_FACE(IRC) )
   DO II=1,IRC
      ALLOCATE(M%RC_FACE(II)%ZZ_FACE(1:N_TOTAL_SCALARS)   ,&
               M%RC_FACE(II)%RHO_D_DZDN(1:N_TOTAL_SCALARS),&
               M%RC_FACE(II)%H_RHO_D_DZDN(1:N_TOTAL_SCALARS))
      M%RC_FACE(II)%ZZ_FACE      = 0._EB
      M%RC_FACE(II)%RHO_D_DZDN   = 0._EB
      M%RC_FACE(II)%H_RHO_D_DZDN = 0._EB
   ENDDO
   IRC = 0
   ! Start with external wall cells:
   GUARD_CUT_CELL_LOOP_1B :  DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      WC=>WALL(IW); BC => BOUNDARY_COORD(WC%BC_INDEX); II = BC%II; JJ = BC%JJ; KK = BC%KK; IOR = BC%IOR
      ! Which face:
      X1AXIS=ABS(IOR)
      SELECT CASE(IOR)
      CASE( IAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND ! Face on high side of Guard-Cell
      CASE(-IAXIS); IIF=II-1; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST= LOW_IND
      CASE( JAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND
      CASE(-JAXIS); IIF=II  ; JJF=JJ-1; KKF=KK  ; LOWHIGH_TEST= LOW_IND
      CASE( KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK  ; LOWHIGH_TEST=HIGH_IND
      CASE(-KAXIS); IIF=II  ; JJF=JJ  ; KKF=KK-1; LOWHIGH_TEST= LOW_IND
      END SELECT

      IF(IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) < 1) CYCLE GUARD_CUT_CELL_LOOP_1B ! Not an RC_FACE.

      ! None of the following defined RCFACES in block boundaries have been added to RC_FACE before:
      ICC   = IJKFACE(IIF,JJF,KKF,X1AXIS,ICPOS)
      JCC   = IJKFACE(IIF,JJF,KKF,X1AXIS,JCPOS)
      IFACE = IJKFACE(IIF,JJF,KKF,X1AXIS,IFPOS)
      IJK(IAXIS:KAXIS) = CUT_CELL(ICC)%IJK(IAXIS:KAXIS)
      ! Which face?
      LOWHIGH = CUT_CELL(ICC)%FACE_LIST(2,IFACE)

      ! Add face to RC face:
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         X2AXIS = JAXIS
         X3AXIS = KAXIS
         ! location in I,J,K od x2,x2,x3 axes:
         XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
      CASE(JAXIS)
         X2AXIS = KAXIS
         X3AXIS = IAXIS
         ! location in I,J,K od x2,x2,x3 axes:
         XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
      CASE(KAXIS)
         X2AXIS = IAXIS
         X3AXIS = JAXIS
         ! location in I,J,K od x2,x2,x3 axes:
         XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
      END SELECT

      IF_LOW_HIGH_1B : IF (LOWHIGH == LOW_IND) THEN

         ! Face indexes:
         INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
         INFACE = INDXI1(XIAXIS)
         JNFACE = INDXI1(XJAXIS)
         KNFACE = INDXI1(XKAXIS)

         ! Location of next Cartesian cell:
         INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
         INCELL = INDXI1(XIAXIS)
         JNCELL = INDXI1(XJAXIS)
         KNCELL = INDXI1(XKAXIS)

         ! Scalar:
         IF (CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_GASPHASE ) THEN ! next cell is reg-cell:

            ! Set OZPOS to 2, to be used in next cycle:
            IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) = 2
            ! Add face to RC_FACE data structure:
            IRC = IRC + 1
            M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)
            M%RC_FACE(IRC)%IWC = IW ! Locate WALL CELL for boundary M%RC_FACE(IRC).

            ! Add all info required for matrix build:
            ! Cell at i-1, i.e. regular GASPHASE:
            M%RC_FACE(IRC)%UNKZ(LOW_IND)                     = CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)         = (/ XCELL(INCELL), YCELL(JNCELL), ZCELL(KNCELL) /)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1, LOW_IND) = (/ CC_FTYPE_RGGAS, INCELL, JNCELL, KNCELL /)

            ! Cell at i+1, i.e. cut-cell:
            M%RC_FACE(IRC)%UNKZ(HIGH_IND)                    = CUT_CELL(ICC)%UNKZ(JCC)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = CUT_CELL(ICC)%XYZCEN(IAXIS:KAXIS,JCC)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)

            ! Modify FACE_LIST for the given cut-cell:
            CUT_CELL(ICC)%FACE_LIST(4,IFACE) = IRC
         ELSEIF(CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE) THEN ! next cell is cc:

            ! Set OZPOS to 2, to be used in next cycle:
            IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) = 2
            ! Add face to RC_FACE  data structure:
            IRC = IRC + 1
            M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)
            M%RC_FACE(IRC)%IWC = IW ! Locate WALL CELL for boundary M%RC_FACE(IRC).

            ! Add all info required for matrix build:
            ! Cell at i+1, i.e. cut-cell:
            M%RC_FACE(IRC)%UNKZ(HIGH_IND)                    = CUT_CELL(ICC)%UNKZ(JCC)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = CUT_CELL(ICC)%XYZCEN(IAXIS:KAXIS,JCC)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)

            ! Modify FACE_LIST for the given cut-cell:
            CUT_CELL(ICC)%FACE_LIST(4,IFACE) = IRC
         ELSE
            WRITE(LU_ERR,*) 'MISSING BOUNDARY RCFACE',IIF,JJF,KKF,X1AXIS
         ENDIF

      ELSE ! IF_LOW_HIGH : HIGH_IND

         ! Face indexes:
         INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS), IJK(X2AXIS), IJK(X3AXIS) /)
         INFACE = INDXI1(XIAXIS)
         JNFACE = INDXI1(XJAXIS)
         KNFACE = INDXI1(XKAXIS)

         ! Location of next Cartesian cell:
         INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)+1, IJK(X2AXIS), IJK(X3AXIS) /)
         INCELL = INDXI1(XIAXIS)
         JNCELL = INDXI1(XJAXIS)
         KNCELL = INDXI1(XKAXIS)

         IF (CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_GASPHASE ) THEN

            ! Set OZPOS to 2, to be used in next cycle:
            IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) = 2
            ! Add face to RC_FACE data structure:
            IRC = IRC + 1
            M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)
            M%RC_FACE(IRC)%IWC = IW ! Locate WALL CELL for boundary M%RC_FACE(IRC).

            ! Add all info required for matrix build:
            ! Cell at i-1, i.e. cut-cell:
            M%RC_FACE(IRC)%UNKZ(LOW_IND)                     = CUT_CELL(ICC)%UNKZ(JCC)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)         = CUT_CELL(ICC)%XYZCEN(IAXIS:KAXIS,JCC)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1, LOW_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)
            ! Modify FACE_LIST for the given cut-cell:
            CUT_CELL(ICC)%FACE_LIST(4,IFACE) = IRC

            ! Cell at i+1, i.e. regular GASPHASE:
            M%RC_FACE(IRC)%UNKZ(HIGH_IND)                    = CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = (/ XCELL(INCELL), YCELL(JNCELL), ZCELL(KNCELL) /)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_RGGAS, INCELL, JNCELL, KNCELL /)

         ELSEIF(CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE) THEN ! next cell is cc:

            ! Set OZPOS to 2, to be used in next cycle:
            IJKFACE(IIF,JJF,KKF,X1AXIS,OZPOS) = 2
            ! Add face to RC_FACE data structure:
            IRC = IRC + 1
            M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)
            M%RC_FACE(IRC)%IWC = IW ! Locate WALL CELL for boundary M%RC_FACE(IRC).

            ! Add high cell info required for matrix build:
            ! Cell at i-1, i.e. cut-cell:
            M%RC_FACE(IRC)%UNKZ(LOW_IND)                     = CUT_CELL(ICC)%UNKZ(JCC)
            M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)         = CUT_CELL(ICC)%XYZCEN(IAXIS:KAXIS,JCC)
            M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1, LOW_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)

            ! Modify FACE_LIST for the given cut-cell:
            CUT_CELL(ICC)%FACE_LIST(4,IFACE) = IRC
           ELSE
              WRITE(LU_ERR,*) 'MISSING BOUNDARY RCFACE',IIF,JJF,KKF,X1AXIS
         ENDIF
      ENDIF IF_LOW_HIGH_1B

   ENDDO GUARD_CUT_CELL_LOOP_1B

   ! Number of RC faces defined in block boundaries:
   M%CC_NBBRCFACE_Z = IRC

   ! Now run regular cut-cell loop to define internal RCFACES:
   DO ICC=1,M%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); IJK(IAXIS:KAXIS) = CC%IJK(IAXIS:KAXIS)
      DO JCC=1,CC%NCELL
         IF(CC%UNKZ(JCC) < 1) CYCLE
         ! Loop faces and test:
         IFC_LOOP : DO IFC=1,CC%CCELEM(1,JCC)
            IFACE = CC%CCELEM(IFC+1,JCC)
            ! If face type in face_list is not CC_FTYPE_RCGAS, drop:
            IF(CC%FACE_LIST(1,IFACE) /= CC_FTYPE_RCGAS) CYCLE IFC_LOOP

            ! Which face?
            LOWHIGH = CC%FACE_LIST(2,IFACE)
            X1AXIS  = CC%FACE_LIST(3,IFACE)

            SELECT CASE(X1AXIS)
            CASE(IAXIS)
               X2AXIS = JAXIS
               X3AXIS = KAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = IAXIS; XJAXIS = JAXIS; XKAXIS = KAXIS
            CASE(JAXIS)
               X2AXIS = KAXIS
               X3AXIS = IAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = KAXIS; XJAXIS = IAXIS; XKAXIS = JAXIS
            CASE(KAXIS)
               X2AXIS = IAXIS
               X3AXIS = JAXIS
               ! location in I,J,K od x2,x2,x3 axes:
               XIAXIS = JAXIS; XJAXIS = KAXIS; XKAXIS = IAXIS
            END SELECT

            IF_LOW_HIGH : IF (LOWHIGH == LOW_IND) THEN

               ! Face indexes:
               INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
               INFACE = INDXI1(XIAXIS)
               JNFACE = INDXI1(XJAXIS)
               KNFACE = INDXI1(XKAXIS)

               ! Location of next Cartesian cell:
               INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)-1, IJK(X2AXIS), IJK(X3AXIS) /)
               INCELL = INDXI1(XIAXIS)
               JNCELL = INDXI1(XJAXIS)
               KNCELL = INDXI1(XKAXIS)

               ! Scalar:
               IF (CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ) > 0 ) THEN ! next cell is reg-cell:

                  IF(IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,OZPOS) == 2) CYCLE IFC_LOOP ! CC-REG face already counted in
                                                                                     ! external boundary loop.

                  ! Add face to RC_FACE data structure:
                  IRC = IRC + 1
                  M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)

                  ! Add all info required for matrix build:
                  ! Cell at i-1, i.e. regular GASPHASE:
                  M%RC_FACE(IRC)%UNKZ(LOW_IND)                     = CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)         = (/ XCELL(INCELL), YCELL(JNCELL), ZCELL(KNCELL) /)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,LOW_IND)  = (/ CC_FTYPE_RGGAS, INCELL, JNCELL, KNCELL /)

                  ! Cell at i+1, i.e. cut-cell:
                  M%RC_FACE(IRC)%UNKZ(HIGH_IND)                    = CC%UNKZ(JCC)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = CC%XYZCEN(IAXIS:KAXIS,JCC)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)

                  ! Modify FACE_LIST for the given cut-cell:
                  CC%FACE_LIST(4,IFACE) = IRC
               ELSEIF(CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE) THEN ! next cell is cc:

                  ! Test that Cut-cell to Cut-cell reg face hasn't been added before:
                  INLIST = .FALSE.
                  DO IIFC=1,IRC
                     IF ( M%RC_FACE(IIFC)%IJK(IAXIS)   /= INFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(JAXIS)   /= JNFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(KAXIS)   /= KNFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(KAXIS+1) /= X1AXIS ) CYCLE
                     INLIST = .TRUE.
                     EXIT
                  ENDDO
                  IF (INLIST) THEN
                     ! This cut-cell is on the high side of face iifc:
                     ! Cell at i+1, i.e. cut-cell:
                     M%RC_FACE(IIFC)%UNKZ(HIGH_IND)                    = CC%UNKZ(JCC)
                     M%RC_FACE(IIFC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = CC%XYZCEN(IAXIS:KAXIS,JCC)
                     M%RC_FACE(IIFC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)
                     ! Modify FACE_LIST for the given cut-cell:
                     CC%FACE_LIST(4,IFACE) = IIFC
                     CYCLE IFC_LOOP
                  ENDIF

                  ! Add face to RC_FACE  data structure:
                  IRC = IRC + 1
                  M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)

                  ! Add all info required for matrix build:
                  ! Cell at i+1, i.e. cut-cell:
                  M%RC_FACE(IRC)%UNKZ(HIGH_IND)                        = CC%UNKZ(JCC)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)            = CC%XYZCEN(IAXIS:KAXIS,JCC)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND)     = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)

                  ! Modify FACE_LIST for the given cut-cell:
                  CC%FACE_LIST(4,IFACE) = IRC
               ENDIF

            ELSE ! IF_LOW_HIGH : HIGH_IND

               ! Face indexes:
               INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS), IJK(X2AXIS), IJK(X3AXIS) /)
               INFACE = INDXI1(XIAXIS)
               JNFACE = INDXI1(XJAXIS)
               KNFACE = INDXI1(XKAXIS)

               ! Location of next Cartesian cell:
               INDXI1(IAXIS:KAXIS) = (/ IJK(X1AXIS)+1, IJK(X2AXIS), IJK(X3AXIS) /)
               INCELL = INDXI1(XIAXIS)
               JNCELL = INDXI1(XJAXIS)
               KNCELL = INDXI1(XKAXIS)

               IF (CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ) > 0 ) THEN

                  IF(IJKFACE(INFACE,JNFACE,KNFACE,X1AXIS,OZPOS) == 2) CYCLE IFC_LOOP ! CC-REG face already counted in
                                                                                     ! external boundary loop.

                  ! Add face to RC_FACE data structure:
                  IRC = IRC + 1
                  M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)

                  ! Add all info required for matrix build:
                  ! Cell at i-1, i.e. cut-cell:
                  M%RC_FACE(IRC)%UNKZ(LOW_IND)                     = CC%UNKZ(JCC)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)         = CC%XYZCEN(IAXIS:KAXIS,JCC)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,LOW_IND)  = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)
                  ! Modify FACE_LIST for the given cut-cell:
                  CC%FACE_LIST(4,IFACE) = IRC

                  ! Cell at i+1, i.e. regular GASPHASE:
                  M%RC_FACE(IRC)%UNKZ(HIGH_IND)                    = CCVAR(INCELL,JNCELL,KNCELL,CC_UNKZ)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,HIGH_IND)        = (/ XCELL(INCELL), YCELL(JNCELL), ZCELL(KNCELL) /)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1,HIGH_IND) = (/ CC_FTYPE_RGGAS, INCELL, JNCELL, KNCELL /)

               ELSEIF(CCVAR(INCELL,JNCELL,KNCELL,CC_CGSC) == CC_CUTCFE) THEN ! next cell is cc:

                  ! Test that Cut-cell to Cut-cell reg face hasn't been added before:
                  INLIST = .FALSE.
                  DO IIFC=1,IRC
                     IF ( M%RC_FACE(IIFC)%IJK(IAXIS)   /= INFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(JAXIS)   /= JNFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(KAXIS)   /= KNFACE ) CYCLE
                     IF ( M%RC_FACE(IIFC)%IJK(KAXIS+1) /= X1AXIS ) CYCLE
                     INLIST = .TRUE.
                     EXIT
                  ENDDO
                  IF (INLIST) THEN
                     ! This cut-cell is on the high side of face iifc:
                     ! Cell at i-1, i.e. cut-cell:
                     M%RC_FACE(IIFC)%UNKZ(LOW_IND)                     = CC%UNKZ(JCC)
                     M%RC_FACE(IIFC)%XCEN(IAXIS:KAXIS,LOW_IND)         = CC%XYZCEN(IAXIS:KAXIS,JCC)
                     M%RC_FACE(IIFC)%CELL_LIST(IAXIS:KAXIS+1, LOW_IND) = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)
                     ! Modify FACE_LIST for the given cut-cell:
                     CC%FACE_LIST(4,IFACE) = IIFC
                     CYCLE IFC_LOOP
                  ENDIF

                  ! Add face to RC_FACE data structure:
                  IRC = IRC + 1
                  M%RC_FACE(IRC)%IJK(IAXIS:KAXIS+1) = (/ INFACE, JNFACE, KNFACE, X1AXIS/)

                  ! Add high cell info required for matrix build:
                  ! Cell at i-1, i.e. cut-cell:
                  M%RC_FACE(IRC)%UNKZ(LOW_IND)                         = CC%UNKZ(JCC)
                  M%RC_FACE(IRC)%XCEN(IAXIS:KAXIS,LOW_IND)             = CC%XYZCEN(IAXIS:KAXIS,JCC)
                  M%RC_FACE(IRC)%CELL_LIST(IAXIS:KAXIS+1, LOW_IND)     = (/ CC_FTYPE_CFGAS, ICC, JCC, IFC /)
                  ! Modify FACE_LIST for the given cut-cell:
                  CC%FACE_LIST(4,IFACE) = IRC
               ENDIF
            ENDIF IF_LOW_HIGH

         ENDDO IFC_LOOP

      ENDDO
   ENDDO

   ! Final number of RC faces:
   M%CC_NRCFACE_Z = IRC

   ! Note WALL_CELLs for internal RC_FACES:
   DO IRC=M%CC_NBBRCFACE_Z+1,M%CC_NRCFACE_Z
      RCF => M%RC_FACE(IRC);
      I = RCF%IJK(IAXIS); J = RCF%IJK(JAXIS); K = RCF%IJK(KAXIS); X1AXIS = RCF%IJK(KAXIS+1)
      ! Don't count cut-faces inside an OBST, or don't lay on a WALL_CELL:
      SELECT CASE(X1AXIS)
      CASE(IAXIS)
         IF (ALL(CELL(CELL_INDEX(I:I+1,J,K))%SOLID) .OR. ALL(.NOT.CELL(CELL_INDEX(I:I+1,J,K))%SOLID)) THEN; CYCLE
         ELSEIF(    CELL(CELL_INDEX(I,J,K))%SOLID .AND. .NOT.CELL(CELL_INDEX(I+1,J,K))%SOLID) THEN
            IW = CELL(CELL_INDEX(I+1,J,K))%WALL_INDEX(-X1AXIS) ! Low face of I+1 cell.
         ELSEIF(.NOT.CELL(CELL_INDEX(I,J,K))%SOLID .AND.     CELL(CELL_INDEX(I+1,J,K))%SOLID) THEN
            IW = CELL(CELL_INDEX(I  ,J,K))%WALL_INDEX( X1AXIS) ! High face of I cell.
         ENDIF
      CASE(JAXIS)
         IF (ALL(CELL(CELL_INDEX(I,J:J+1,K))%SOLID) .OR. ALL(.NOT.CELL(CELL_INDEX(I,J:J+1,K))%SOLID)) THEN; CYCLE
         ELSEIF(    CELL(CELL_INDEX(I,J,K))%SOLID .AND. .NOT.CELL(CELL_INDEX(I,J+1,K))%SOLID) THEN
            IW = CELL(CELL_INDEX(I,J+1,K))%WALL_INDEX(-X1AXIS) ! Low face of J+1 cell.
         ELSEIF(.NOT.CELL(CELL_INDEX(I,J,K))%SOLID .AND.     CELL(CELL_INDEX(I,J+1,K))%SOLID) THEN
            IW = CELL(CELL_INDEX(I,J  ,K))%WALL_INDEX( X1AXIS) ! High face of J cell.
         ENDIF
      CASE(KAXIS)
         IF (ALL(CELL(CELL_INDEX(I,J,K:K+1))%SOLID) .OR. ALL(.NOT.CELL(CELL_INDEX(I,J,K:K+1))%SOLID)) THEN; CYCLE
         ELSEIF(    CELL(CELL_INDEX(I,J,K))%SOLID .AND. .NOT.CELL(CELL_INDEX(I,J,K+1))%SOLID) THEN
            IW = CELL(CELL_INDEX(I,J,K+1))%WALL_INDEX(-X1AXIS) ! Low face of K+1 cell.
         ELSEIF(.NOT.CELL(CELL_INDEX(I,J,K))%SOLID .AND.     CELL(CELL_INDEX(I,J,K+1))%SOLID) THEN
            IW = CELL(CELL_INDEX(I,J,K  ))%WALL_INDEX( X1AXIS) ! High face of K cell.
         ENDIF
      END SELECT
      IF(IW>0) RCF%IWC=IW
   ENDDO

   ! Cell centered positions and cell sizes:
   IF (ALLOCATED(XCELL)) DEALLOCATE(XCELL)
   IF (ALLOCATED(YCELL)) DEALLOCATE(YCELL)
   IF (ALLOCATED(ZCELL)) DEALLOCATE(ZCELL)
   IF (ALLOCATED(IJKBUFFER)) DEALLOCATE(IJKBUFFER)
   IF (ALLOCATED(LOHIBUFF))  DEALLOCATE(LOHIBUFF)
   IF (ALLOCATED(IJKFACE))   DEALLOCATE(IJKFACE)

ENDDO MAIN_MESH_LOOP


IF (DEBUG_MATVEC_DATA) THEN
   DBG_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      IF(MY_RANK/=PROCESS(NM)) CYCLE DBG_MESH_LOOP
      CALL POINT_TO_MESH(NM)
      WRITE(LU_ERR,*) ' '
      WRITE(LU_ERR,*) 'MY_RANK, NM : ',MY_RANK,NM
      WRITE(LU_ERR,*) 'CC_NBBREGFACE(1:3) : ',M%CC_NBBREGFACE_Z(IAXIS:KAXIS)
      WRITE(LU_ERR,*) 'CC_NREGFACE(1:3)   : ',M%CC_NREGFACE_Z(IAXIS:KAXIS)
      WRITE(LU_ERR,*) 'CC_NRCFACE_Z, CC_NBBRCFACE_Z : ', &
                       M%CC_NRCFACE_Z,M%CC_NBBRCFACE_Z
   ENDDO DBG_MESH_LOOP
ENDIF

RETURN
END SUBROUTINE GET_GASPHASE_REGRCFACES_DATA



! ----------------------------- GET_LINKED_MATRIX_INDEXES_Z ---------------------------------
SUBROUTINE GET_LINKED_MATRIX_INDEXES_Z
USE MPI_F08

! Local variables:
INTEGER :: NM
INTEGER :: X1AXIS,I,J,K,ICC,JCC,ICC2,JCC2,ILEV,INGH,JNGH,KNGH,IERR

! Linking variables associated data:
INTEGER, ALLOCATABLE, DIMENSION(:) :: CELLPUNKZ, INDUNKZ
INTEGER :: COUNT, ICF, CF_STATUS

! Define local number of cut-cell:
IF (ALLOCATED(NUNKZ_LOC)) DEALLOCATE(NUNKZ_LOC)
ALLOCATE(NUNKZ_LOC(1:NMESHES)); NUNKZ_LOC = 0

! Cell numbers for Scalar equations:
MAIN_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)

   ! Reset UNKZ to CC_UNDEFINED:
   CCVAR(:,:,:,CC_UNKZ) = CC_UNDEFINED
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CUT_CELL(ICC)%UNKZ(:) = CC_UNDEFINED
   ENDDO

   ! 1. Number regular GASPHASE cells:
   IF (PERIODIC_TEST==103 .OR. PERIODIC_TEST==11 .OR. PERIODIC_TEST==7) THEN
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               ! If regular cell centroid is outside the test box + DELTA -> drop:
               IF(XC(I) < (VAL_TESTX_LOW-DX(I) +GEOMEPS)) CYCLE; IF(XC(I) > (VAL_TESTX_HIGH+DX(I)-GEOMEPS)) CYCLE
               IF(YC(J) < (VAL_TESTY_LOW-DY(J) +GEOMEPS)) CYCLE; IF(YC(J) > (VAL_TESTY_HIGH+DY(J)-GEOMEPS)) CYCLE
               IF(ZC(K) < (VAL_TESTZ_LOW-DZ(K) +GEOMEPS)) CYCLE; IF(ZC(K) > (VAL_TESTZ_HIGH+DZ(K)-GEOMEPS)) CYCLE
               IF(CCVAR(I,J,K,CC_CGSC) /= CC_GASPHASE) CYCLE
               NUNKZ_LOC(NM) = NUNKZ_LOC(NM) + 1
               CCVAR(I,J,K,CC_UNKZ) = NUNKZ_LOC(NM)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Loop on cut-cells and surrounding cartesian cells, number unknowns for cells in LINK_LEV=0:
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            ! Drop if cartesian cell is not type CC_CUTCFE:
            IF (  CCVAR(I,J,K,CC_CGSC) /= CC_CUTCFE ) CYCLE
            ! First Add the Cut-Cell:
            ICC  = CCVAR(I,J,K,CC_IDCC) ! The following test excludes GC cut-cells from numbering.
            IF (ICC <= MESHES(NM)%N_CUTCELL_MESH .AND. .NOT.CELL(CELL_INDEX(I,J,K))%SOLID ) THEN
               DO JCC=1,CUT_CELL(ICC)%NCELL
                  IF ( CUT_CELL(ICC)%LINK_LEV(JCC)/=0) CYCLE ! Linked cell, dealt with later.
                  NUNKZ_LOC(NM) = NUNKZ_LOC(NM) + 1
                  CUT_CELL(ICC)%UNKZ(JCC) = NUNKZ_LOC(NM)
               ENDDO
            ENDIF
            ! Surrounding regular cells:
            DO KNGH=K-1,K+1
               IF ( (KNGH < 1) .OR. (KNGH > KBAR) ) CYCLE
               DO JNGH=J-1,J+1
                  IF ( (JNGH < 1) .OR. (JNGH > JBAR) ) CYCLE
                  DO INGH=I-1,I+1
                     ! Either not GASPHASE or already counted:
                     IF ((CCVAR(INGH,JNGH,KNGH,CC_CGSC)/=CC_GASPHASE) .OR. (CCVAR(INGH,JNGH,KNGH,CC_UNKZ)>0)) CYCLE
                     IF ( (INGH < 1) .OR. (INGH > IBAR) ) CYCLE
                     IF (CELL(CELL_INDEX(INGH,JNGH,KNGH))%SOLID) CYCLE
                     ! Add Scalar unknown:
                     NUNKZ_LOC(NM) = NUNKZ_LOC(NM) + 1
                     CCVAR(INGH,JNGH,KNGH,CC_UNKZ) = NUNKZ_LOC(NM)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ! The do all link tree levels from -1 to FINEST_LINK_LEV:
   LINK_LEV_DO : DO ILEV=-1,MESHES(NM)%FINEST_LINK_LEV,-1
      DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
         DO JCC=1,CUT_CELL(ICC)%NCELL
            IF (CUT_CELL(ICC)%LINK_LEV(JCC) /= ILEV) CYCLE

            ! Find master cell for this CC:
            I = CUT_CELL(ICC)%IJK_LINK(2,JCC); J = CUT_CELL(ICC)%IJK_LINK(3,JCC); K = CUT_CELL(ICC)%IJK_LINK(4,JCC)
            SELECT CASE(CUT_CELL(ICC)%IJK_LINK(1,JCC))
            CASE(CC_GASPHASE)
               CUT_CELL(ICC)%UNKZ(JCC) = CCVAR(I,J,K,CC_UNKZ)
            CASE(CC_CUTCFE)
               ICC2 = CCVAR(I,J,K,CC_IDCC); JCC2 = CUT_CELL(ICC)%IJK_LINK(5,JCC)
               CUT_CELL(ICC)%UNKZ(JCC) = CUT_CELL(ICC2)%UNKZ(JCC2)
            END SELECT

         ENDDO
      ENDDO
   ENDDO LINK_LEV_DO

ENDDO MAIN_MESH_LOOP

! After fixing cut-cell unkz for a given Cartesian cells there might be UNKZ values that haven't been assigned.
! Condense:
REIND_MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   IF(NUNKZ_LOC(NM) == 0) CYCLE
   CALL POINT_TO_MESH(NM)
   ALLOCATE(CELLPUNKZ(1:NUNKZ_LOC(NM)), INDUNKZ(1:NUNKZ_LOC(NM))); CELLPUNKZ = 0; INDUNKZ = 0;
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CCVAR(I,J,K,CC_UNKZ) > 0) CELLPUNKZ(CCVAR(I,J,K,CC_UNKZ)) = CELLPUNKZ(CCVAR(I,J,K,CC_UNKZ)) + 1
         ENDDO
      ENDDO
   ENDDO
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      DO JCC=1,CUT_CELL(ICC)%NCELL
         IF (CUT_CELL(ICC)%UNKZ(JCC) > 0) CELLPUNKZ(CUT_CELL(ICC)%UNKZ(JCC)) = CELLPUNKZ(CUT_CELL(ICC)%UNKZ(JCC)) + 1
      ENDDO
   ENDDO
   ! Now re-index:
   COUNT=0
   DO I=1,NUNKZ_LOC(NM)
      IF(CELLPUNKZ(I) == 0) CYCLE ! This UNKZ_LOC value has no cells associated to it.
      COUNT = COUNT + 1; INDUNKZ(I) = COUNT
   ENDDO
   NUNKZ_LOC(NM) = COUNT
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF(CCVAR(I,J,K,CC_UNKZ) > 0) CCVAR(I,J,K,CC_UNKZ) = INDUNKZ(CCVAR(I,J,K,CC_UNKZ)) ! Condensed value.
         ENDDO
      ENDDO
   ENDDO
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      DO JCC=1,CUT_CELL(ICC)%NCELL
         IF (CUT_CELL(ICC)%UNKZ(JCC) > 0) CUT_CELL(ICC)%UNKZ(JCC) = INDUNKZ(CUT_CELL(ICC)%UNKZ(JCC)) ! Condensed value.
      ENDDO
   ENDDO
   DEALLOCATE(CELLPUNKZ,INDUNKZ)
ENDDO REIND_MESH_LOOP

! Define total number of unknowns and global unknown index start per MESH:
IF (ALLOCATED(NUNKZ_TOT)) DEALLOCATE(NUNKZ_TOT)
ALLOCATE(NUNKZ_TOT(1:NMESHES)); NUNKZ_TOT = 0
IF (N_MPI_PROCESSES > 1) THEN
   CALL MPI_ALLREDUCE(NUNKZ_LOC(1), NUNKZ_TOT(1), NMESHES, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERR)
ELSE
   NUNKZ_TOT = NUNKZ_LOC
ENDIF
! Define global start indexes for each mesh:
IF (ALLOCATED(UNKZ_ILC)) DEALLOCATE(UNKZ_ILC)
ALLOCATE(UNKZ_ILC(1:NMESHES)); UNKZ_ILC(1:NMESHES) = 0
IF (ALLOCATED(UNKZ_IND)) DEALLOCATE(UNKZ_IND)
ALLOCATE(UNKZ_IND(1:NMESHES)); UNKZ_IND(1:NMESHES) = 0
DO NM=2,NMESHES
   UNKZ_ILC(NM) = UNKZ_ILC(NM-1) + NUNKZ_LOC(NM-1)
   UNKZ_IND(NM) = UNKZ_IND(NM-1) + NUNKZ_TOT(NM-1)
ENDDO

! Cell numbers for Scalar equations in global numeration:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   ! 1. Number regular GASPHASE cells within the implicit region:
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (CCVAR(I,J,K,CC_CGSC) /= CC_GASPHASE .OR. CCVAR(I,J,K,CC_UNKZ) <= 0 ) CYCLE
            CCVAR(I,J,K,CC_UNKZ) = CCVAR(I,J,K,CC_UNKZ) + UNKZ_IND(NM)
         ENDDO
      ENDDO
   ENDDO
   ! 2. Number cut-cells:
   DO ICC=1,MESHES(NM)%N_CUTCELL_MESH
      CC => CUT_CELL(ICC); I = CC%IJK(IAXIS); J = CC%IJK(JAXIS); K = CC%IJK(KAXIS)
      IF (CELL(CELL_INDEX(I,J,K))%SOLID) CYCLE
      DO JCC=1,CC%NCELL; CC%UNKZ(JCC) = CC%UNKZ(JCC) + UNKZ_IND(NM); ENDDO
   ENDDO
ENDDO

! Exchange Guardcell + guard cc information on CC_UNKZ:
CALL FILL_UNKZ_GUARDCELLS

! Finally set to solid Gasphase cut-faces which have a surrounding cut-cell inside an OBST:
DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
   CALL POINT_TO_MESH(NM)
   DO ICF=1,MESHES(NM)%N_CUTFACE_MESH
      IF (CUT_FACE(ICF)%STATUS /= CC_GASPHASE) CYCLE
      I=CUT_FACE(ICF)%IJK(IAXIS); J=CUT_FACE(ICF)%IJK(JAXIS); K=CUT_FACE(ICF)%IJK(KAXIS)
      X1AXIS=CUT_FACE(ICF)%IJK(KAXIS+1)
      CF_STATUS = CC_GASPHASE
      SELECT CASE(X1AXIS)
      CASE(IAXIS); IF ( CELL(CELL_INDEX(I,J,K))%SOLID .AND. CELL(CELL_INDEX(I+1,J,K))%SOLID ) CF_STATUS = CC_SOLID
      CASE(JAXIS); IF ( CELL(CELL_INDEX(I,J,K))%SOLID .AND. CELL(CELL_INDEX(I,J+1,K))%SOLID ) CF_STATUS = CC_SOLID
      CASE(KAXIS); IF ( CELL(CELL_INDEX(I,J,K))%SOLID .AND. CELL(CELL_INDEX(I,J,K+1))%SOLID ) CF_STATUS = CC_SOLID
      END SELECT
      CUT_FACE(ICF)%STATUS = CF_STATUS
   ENDDO
ENDDO

RETURN
END SUBROUTINE GET_LINKED_MATRIX_INDEXES_Z



SUBROUTINE GET_LINKED_FACE_INDEXES_F

INTEGER :: NM,I,J,K,X1AXIS,X2AXIS,X3AXIS,ICF,JCF,LO_UNKZ,HI_UNKZ,IEC,IEDGE,LOHI,II,JJ,KK,IIO,JJO,KKO,IERC,&
           OFACE(3),OLO_UNKZ,OHI_UNKZ,OICF,OJCF,IECE,JECE,ILINK,ICC,JCC,IW,COUNT
REAL(EB):: ACRT,CCVOL_THRES
LOGICAL :: ALL_FLG,CC_LINKED
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: ILOC,SIZE_FACE
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: FACE_LIST,FACELAUX
REAL(EB), ALLOCATABLE,DIMENSION(:)   :: FACE_AREA,FACEARAUX
TYPE(CC_CUTFACE_TYPE),       POINTER :: CF2

LOGICAL, PARAMETER :: NO_FACE_LINKING = .FALSE.

SIZE_FACE = 20; ALLOCATE(FACE_LIST(4,SIZE_FACE),FACE_AREA(SIZE_FACE))

! Define Face Linking:
! Important approximation: As we do not compute cut-face volumes from the computational geometry engine (4 times cost),
! face volumes are assumed to be composed of Area*DXN; DXN is the cartesian cell size in the normal face direction.
MESH_LOOP : DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX

   M => MESHES(NM)

   ! 1. Run across gasphase cut-faces and number those which have Area*DXN > CCVOL_LINK*DX(I)*DX(2)*DX(3):
   M%NUNK_F = 0
   FACE_LINK_IF : IF (NO_FACE_LINKING) THEN

      ICF_LOOP_0 : DO ICF=1,M%N_CUTFACE_MESH
         CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE ICF_LOOP_0
         IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE ICF_LOOP_0
         DO JCF=1,CF%NFACE; M%NUNK_F = M%NUNK_F+1; CF%UNKF(JCF)=M%NUNK_F; CF%LINK_LEV(JCF)=0; ENDDO
      ENDDO ICF_LOOP_0

   ELSE FACE_LINK_IF

   ICF_LOOP_1 : DO ICF=1,M%N_CUTFACE_MESH
      CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE ICF_LOOP_1
      IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE ICF_LOOP_1
      ! Give the mesh boundary cut-faces with boundary not periodic or interpolated their own unlinked UNKF.
      IF(CF%IWC>0) THEN
         IF( .NOT.ANY(M%WALL(CF%IWC)%BOUNDARY_TYPE == (/PERIODIC_BOUNDARY,INTERPOLATED_BOUNDARY,OPEN_BOUNDARY/)) ) THEN
            DO JCF=1,CF%NFACE; M%NUNK_F = M%NUNK_F+1; CF%UNKF(JCF)=M%NUNK_F; CF%LINK_LEV(JCF)=0; ENDDO
            CYCLE ICF_LOOP_1
         ENDIF
      ENDIF
      I  = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS)
      SELECT CASE(CF%IJK(KAXIS+1))
      CASE(IAXIS); ACRT = M%DY(J)*M%DZ(K)
      CASE(JAXIS); ACRT = M%DX(I)*M%DZ(K)
      CASE(KAXIS); ACRT = M%DY(J)*M%DX(I)
      END SELECT
      DO JCF=1,CF%NFACE
         ! Test if any of surrounding cut-cells are linked. If any link attempt to link cut-face.
         CC_LINKED = .FALSE.
         DO LOHI=LOW_IND,HIGH_IND
            ICC = CF%CELL_LIST(2,LOHI,JCF); JCC = CF%CELL_LIST(3,LOHI,JCF); CC => M%CUT_CELL(ICC)
            CCVOL_THRES = CCVOL_LINK * (M%DX(CC%IJK(IAXIS))*M%DY(CC%IJK(JAXIS))*M%DZ(CC%IJK(KAXIS)))
            IF ( CC%VOLUME(JCC) <= CCVOL_THRES ) THEN; CC_LINKED = .TRUE.; EXIT; ENDIF
         ENDDO
         ! If surrounding cut-cells not linked, face UNKF not numbered and Area over threshold, add to NUNK_F:
         IF(CF%UNKF(JCF)<1 .AND. CF%AREA(JCF)>CCVOL_LINK*ACRT .AND. .NOT.CC_LINKED) THEN
            M%NUNK_F = M%NUNK_F + 1
            CF%UNKF(JCF) = M%NUNK_F
            CF%LINK_LEV(JCF) = 0
         ENDIF
      ENDDO
   ENDDO ICF_LOOP_1

   ! 2. Link small faces with neighbor faces that share CV, and then Link small faces with neighbor faces:
   ALL_FLG=.FALSE.
   LINK_ITER : DO ILINK=1,N_LINK_ATTMP_F
      IF(ILINK>3) ALL_FLG=.TRUE.
      CALL LINK_FACES_TO_NEIGHBORS(ALL_FLG)
      ! Test For remaining unlinked faces:
      DO ICF=1,M%N_CUTFACE_MESH
         CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
         IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE
         DO JCF=1,CF%NFACE
            IF(CF%UNKF(JCF)<1) CYCLE LINK_ITER
         ENDDO
      ENDDO
      EXIT LINK_ITER
   ENDDO LINK_ITER

   ! Finally force link small unlinked faces or set their UNKF, LINK_LEV=0.
   DO ICF=1,M%N_CUTFACE_MESH
      CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
      IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE
      I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
      DO JCF=1,CF%NFACE
         IF(CF%UNKF(JCF)>0) CYCLE
         NULLIFY(CF2)
         SELECT CASE(X1AXIS)
         CASE(IAXIS)
            IF(I>0) THEN
               IF(M%FCVAR(I-1,J,K,CC_IDCF,IAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I-1,J,K,CC_IDCF,IAXIS))
               ELSEIF(M%FCVAR(I-1,J,K,CC_FGSC,IAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I-1,J,K,CC_UNKF,IAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I-1,J,K,CC_UNKF,IAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I-1,J,K,CC_UNKF,IAXIS); CYCLE
               ENDIF
            ENDIF
            IF(.NOT.ASSOCIATED(CF2) .AND. I<M%IBAR) THEN
               IF(M%FCVAR(I+1,J,K,CC_IDCF,IAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I+1,J,K,CC_IDCF,IAXIS))
               ELSEIF(M%FCVAR(I+1,J,K,CC_FGSC,IAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I+1,J,K,CC_UNKF,IAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I+1,J,K,CC_UNKF,IAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I+1,J,K,CC_UNKF,IAXIS); CYCLE
               ENDIF
            ENDIF
         CASE(JAXIS)
            IF(J>0) THEN
               IF(M%FCVAR(I,J-1,K,CC_IDCF,JAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I,J-1,K,CC_IDCF,JAXIS))
               ELSEIF(M%FCVAR(I,J-1,K,CC_FGSC,JAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I,J-1,K,CC_UNKF,JAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I,J-1,K,CC_UNKF,JAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I,J-1,K,CC_UNKF,JAXIS); CYCLE
               ENDIF
            ENDIF
            IF(.NOT.ASSOCIATED(CF2) .AND. J<M%JBAR) THEN
               IF(M%FCVAR(I,J+1,K,CC_IDCF,JAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I,J+1,K,CC_IDCF,JAXIS))
               ELSEIF(M%FCVAR(I,J+1,K,CC_FGSC,JAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I,J+1,K,CC_UNKF,JAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I,J+1,K,CC_UNKF,JAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I,J+1,K,CC_UNKF,JAXIS); CYCLE
               ENDIF
            ENDIF
         CASE(KAXIS)
            IF(K>0) THEN
               IF(M%FCVAR(I,J,K-1,CC_IDCF,KAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I,J,K-1,CC_IDCF,KAXIS))
               ELSEIF(M%FCVAR(I,J,K-1,CC_FGSC,KAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I,J,K-1,CC_UNKF,KAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I,J,K-1,CC_UNKF,KAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I,J,K-1,CC_UNKF,KAXIS); CYCLE
               ENDIF
            ENDIF
            IF(.NOT.ASSOCIATED(CF2) .AND. K<M%KBAR) THEN
               IF(M%FCVAR(I,J,K+1,CC_IDCF,KAXIS)>0) THEN
                  CF2=>M%CUT_FACE(M%FCVAR(I,J,K+1,CC_IDCF,KAXIS))
               ELSEIF(M%FCVAR(I,J,K+1,CC_FGSC,KAXIS)==CC_GASPHASE) THEN
                  IF(M%FCVAR(I,J,K+1,CC_UNKF,KAXIS)<1)THEN; M%NUNK_F=M%NUNK_F+1; M%FCVAR(I,J,K+1,CC_UNKF,KAXIS)=M%NUNK_F; ENDIF
                  CF%UNKF(JCF)=M%FCVAR(I,J,K+1,CC_UNKF,KAXIS); CYCLE
               ENDIF
            ENDIF
         END SELECT
         IF(ASSOCIATED(CF2)) THEN
            IF(CF2%UNKF(1)>0) THEN; CF%UNKF(JCF)=CF2%UNKF(1); CF%LINK_LEV(JCF)=CF2%LINK_LEV(1)-1; ENDIF
         ENDIF
         IF(CF%UNKF(JCF)<1) THEN; CF%UNKF(JCF)=0; CF%LINK_LEV(JCF)=0; ENDIF
      ENDDO
   ENDDO


   ! I = 0; J = 0; K = 0
   ! DO ICF=1,M%N_CUTFACE_MESH
   !    CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
   !    IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE
   !    DO JCF=1,CF%NFACE
   !       IF(CF%UNKF(JCF)>0) THEN
   !          IF(CF%IJK(KAXIS+1)==IAXIS) I=I+1
   !          IF(CF%IJK(KAXIS+1)==JAXIS) J=J+1
   !          IF(CF%IJK(KAXIS+1)==KAXIS) K=K+1
   !          WRITE(LU_ERR,*) ICF,JCF,CF%IJK(1:4),CF%UNKF(JCF)
   !       ENDIF
   !    ENDDO
   ! ENDDO
   ! WRITE(LU_ERR,*) 'N_CUTFACE_MESH,Unassigned CFs,NUNK_F=',ILINK,':',M%N_CUTFACE_MESH,I,J,K,I+J+K,';',M%NUNK_F
   !
   ! DO K=1,M%KBAR
   !    DO J=1,M%JBAR
   !       DO I=0,M%IBAR
   !          IF(M%FCVAR(I,J,K,CC_UNKF,IAXIS)>0) WRITE(LU_ERR,*) 'CRT=',I,J,K,IAXIS,M%FCVAR(I,J,K,CC_UNKF,X1AXIS)
   !       ENDDO
   !    ENDDO
   ! ENDDO

  ENDIF FACE_LINK_IF

  ! Allocate linked face velocity arrays:
  IF(ALLOCATED(M%EWC_UN_LNK)) DEALLOCATE(M%EWC_UN_LNK)
  IF(ALLOCATED(M%UN_LNK)) DEALLOCATE(M%UN_LNK)
  ALLOCATE(M%EWC_UN_LNK(0:M%N_EXTERNAL_WALL_CELLS)); M%EWC_UN_LNK = 0._EB
  ALLOCATE(M%UN_LNK(0:M%NUNK_F)); M%UN_LNK = 0._EB
  DO IW=1,M%N_EXTERNAL_WALL_CELLS
     IF(M%WALL(IW)%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY .OR. M%WALL(IW)%CUT_FACE_INDEX<1) CYCLE
     CF => M%CUT_FACE(M%WALL(IW)%CUT_FACE_INDEX)
     IF(.NOT.ALLOCATED(CF%VEL_LNK)) ALLOCATE(CF%VEL_LNK(1:CF%NFACE)); CF%VEL_LNK = 0._EB
  ENDDO

  ! Allocate UN_ULNK:
  COUNT=0
  DO K=0,M%KBAR
     DO J=0,M%JBAR
        DO I=0,M%IBAR
           IF (M%FCVAR(I,J,K,CC_UNKF,IAXIS)>0) COUNT = COUNT+1
           IF (M%FCVAR(I,J,K,CC_UNKF,JAXIS)>0 .AND. .NOT.TWO_D) COUNT = COUNT+1
           IF (M%FCVAR(I,J,K,CC_UNKF,KAXIS)>0) COUNT = COUNT+1
        ENDDO
     ENDDO
  ENDDO
  DO ICF=1,M%CC_NRCFACE_Z
     IF(M%RC_FACE(ICF)%UNKF<1) CYCLE
     I = M%RC_FACE(ICF)%IJK(IAXIS); J = M%RC_FACE(ICF)%IJK(JAXIS); K = M%RC_FACE(ICF)%IJK(KAXIS)
     X1AXIS = M%RC_FACE(ICF)%IJK(KAXIS+1)
     SELECT CASE(X1AXIS)
     CASE(IAXIS); COUNT = COUNT+1
     CASE(JAXIS); IF(.NOT.TWO_D) COUNT = COUNT+1
     CASE(KAXIS); COUNT = COUNT+1
     END SELECT
  ENDDO
  DO ICF=1,M%N_CUTFACE_MESH
     CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE
     IF(TWO_D .AND. CF%IJK(KAXIS+1)==JAXIS) CYCLE
     I = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
     DO JCF=1,CF%NFACE
        IF (CF%UNKF(JCF)<1) CYCLE
        COUNT = COUNT+1
     ENDDO
  ENDDO
  DO IW=1,M%N_EXTERNAL_WALL_CELLS
     WC=>M%WALL(IW)
     IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
     EWC=>M%EXTERNAL_WALL(IW)
     BC =>M%BOUNDARY_COORD(WC%BC_INDEX)
     I  = BC%II; J = BC%JJ; K = BC%KK; X1AXIS = ABS(BC%IOR)
     SELECT CASE(BC%IOR)
     CASE(-IAXIS); I=I-1
     CASE(-JAXIS); J=J-1
     CASE(-KAXIS); K=K-1
     END SELECT
     IF(M%FCVAR(I,J,K,CC_IDCF,X1AXIS)>0) THEN ! Cut-face.
        ICF=M%FCVAR(I,J,K,CC_IDCF,X1AXIS)
        COUNT=COUNT+M%CUT_FACE(ICF)%NFACE
     ELSE ! All other reg faces.
        COUNT = COUNT + 1
     ENDIF
  ENDDO
  IF(ALLOCATED(M%UN_ULNK)) DEALLOCATE(M%UN_ULNK)
  ALLOCATE(M%UN_ULNK(COUNT)); M%UN_ULNK = 0._EB

ENDDO MESH_LOOP
DEALLOCATE(FACE_LIST,FACE_AREA)

RETURN
CONTAINS

SUBROUTINE LINK_FACES_TO_NEIGHBORS(ALL_FLG)

LOGICAL, INTENT(IN) :: ALL_FLG

LOGICAL, PARAMETER  :: CANDIDATE_LIST=.TRUE.

IF(.NOT.CANDIDATE_LIST) THEN

ICF_LOOP_2 : DO ICF=1,M%N_CUTFACE_MESH
   CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE ICF_LOOP_2
   I  = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   SELECT CASE(X1AXIS)
   CASE(IAXIS) ! Face in X axis:
       ACRT = M%DY(J)*M%DZ(K)
       IAXIS_JCF_LOOP : DO JCF=1,CF%NFACE
          IF(CF%UNKF(JCF)<1) THEN ! .AND. CF%AREA(JCF) < CCVOL_LINK*ACRT+TWENTY_EPSILON_EB) THEN
             LO_UNKZ = CF%UNKZ(LOW_IND, JCF)
             HI_UNKZ = CF%UNKZ(HIGH_IND,JCF)
             ! Loop edges to find next face in X2AXIS, X3AXIS plane:
             IAXIS_IEC_LOOP : DO IEC=2,CF%CEDGES(1,JCF)+1
                IEDGE = CF%CEDGES(IEC,JCF)
                IAXIS_IEC_SELECT : SELECT CASE(CF%EDGE_LIST(1,IEDGE))
                CASE(CC_ETYPE_RGGAS) IAXIS_IEC_SELECT ! Cut-faces Regular Gas Edge.
                   LOHI   = CF%EDGE_LIST(2,IEDGE)
                   X2AXIS = CF%EDGE_LIST(3,IEDGE)
                   IF (X2AXIS==JAXIS) THEN ! Edge pointed in X3AXIS=KAXIS direction.
                      X3AXIS = KAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) JJ=J-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE IAXIS_IEC_LOOP ! Index in CC_RCEDGE
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ J-1
                         IIO=I; JJO=J-1; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-2) ! X Face in low J
                      ELSE ! Indexes of other Cartesian face @ J+1
                         IIO=I; JJO=J+1; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 2) ! X Face in high J
                      ENDIF
                      IF(JJO<1 .OR. JJO>M%JBAR) CYCLE IAXIS_IEC_LOOP
                   ELSEIF(X2AXIS==KAXIS) THEN ! Edge pointed in X3AXIS=JAXIS direction.
                      X3AXIS = JAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) KK=K-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE IAXIS_IEC_LOOP
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ K-1
                         IIO=I; JJO=J; KKO=K-1; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-1) ! Face in low K
                      ELSE ! Indexes of other Cartesian face @ K+1
                         IIO=I; JJO=J; KKO=K+1; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 1) ! Face in high K
                      ENDIF
                      IF(KKO<1 .OR. KKO>M%KBAR) CYCLE IAXIS_IEC_LOOP
                   ENDIF
                   ! Now look for other potential face and link CF%UNKF(JCF):
                   IAXIS_OFACE_SELECT : SELECT CASE (OFACE(1))
                   CASE(CC_FTYPE_RGGAS) IAXIS_OFACE_SELECT
                      OLO_UNKZ = M%CCVAR(IIO  ,JJO  ,KKO  ,CC_UNKZ)
                      OHI_UNKZ = M%CCVAR(IIO+1,JJO  ,KKO  ,CC_UNKZ)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to regular face:
                         IF(M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS) = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)
                         CF%LINK_LEV(JCF) = -1
                         CYCLE IAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_RCGAS) IAXIS_OFACE_SELECT
                      OICF     = OFACE(2)
                      OLO_UNKZ = M%RC_FACE(OICF)%UNKZ(LOW_IND)
                      OHI_UNKZ = M%RC_FACE(OICF)%UNKZ(HIGH_IND)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to RC face:
                         IF(M%RC_FACE(OICF)%UNKF<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%RC_FACE(OICF)%UNKF = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%RC_FACE(OICF)%UNKF
                         CF%LINK_LEV(JCF) = -1
                         CYCLE IAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_CFGAS) IAXIS_OFACE_SELECT
                      OICF     = OFACE(2); OJCF     = OFACE(3)
                      OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                      OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                          IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                             CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                             CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                             CYCLE IAXIS_JCF_LOOP
                          ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                             CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                             CYCLE IAXIS_JCF_LOOP
                          ENDIF
                      ENDIF
                   END SELECT IAXIS_OFACE_SELECT

                CASE(CC_ETYPE_CFGAS) IAXIS_IEC_SELECT ! Gas cut-edge.
                   IECE = CF%EDGE_LIST(2,IEDGE); JECE = CF%EDGE_LIST(3,IEDGE)
                   II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                   X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                   IF (X3AXIS==KAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = JAXIS
                      IF(JJ==J-1) THEN ! LOWER cut-face in the J direction.
                         IIO = I; JJO = J-1; KKO = K ! Indexes of other Cartesian face @ J-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-2, JECE) ! X Face in low J
                      ELSEIF(JJ==J) THEN ! UPPER cut-face in the J direction.
                         IIO = I; JJO = J+1; KKO = K ! Indexes of other Cartesian face @ J+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 2, JECE) ! X Face in high J
                      ENDIF
                      IF(JJO<1 .OR. JJO>M%JBAR) CYCLE IAXIS_IEC_LOOP
                   ELSEIF (X3AXIS==JAXIS) THEN ! X2AXIS = KAXIS
                      IF(KK==K-1) THEN ! LOWER cut-face in the K direction.
                         IIO = I; JJO = J; KKO = K-1 ! Indexes of other Cartesian face @ K-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-1, JECE) ! X Face in low K
                      ELSEIF(KK==K) THEN ! UPPER cut-face in the K direction.
                         IIO = I; JJO = J; KKO = K+1 ! Indexes of other Cartesian face @ K+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 1, JECE) ! X Face in high K
                      ENDIF
                      IF(KKO<1 .OR. KKO>M%KBAR) CYCLE IAXIS_IEC_LOOP
                   ENDIF
                   OICF     = OFACE(2); OJCF     = OFACE(3)
                   IF(OICF<1 .OR. OJCF<1) CYCLE IAXIS_IEC_LOOP
                   OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                   OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                   IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                       IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                          CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                          CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                          CYCLE IAXIS_JCF_LOOP
                       ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                          M%NUNK_F = M%NUNK_F + 1
                          CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                          CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                          CYCLE IAXIS_JCF_LOOP
                       ENDIF
                   ENDIF
                END SELECT IAXIS_IEC_SELECT
             ENDDO IAXIS_IEC_LOOP
          ENDIF
       ENDDO IAXIS_JCF_LOOP
   CASE(JAXIS)
       IF (TWO_D) CYCLE ICF_LOOP_2
       ACRT = M%DX(I)*M%DZ(K)
       JAXIS_JCF_LOOP : DO JCF=1,CF%NFACE
          IF(CF%UNKF(JCF)<1) THEN ! .AND. CF%AREA(JCF) < CCVOL_LINK*ACRT+TWENTY_EPSILON_EB) THEN
             LO_UNKZ = CF%UNKZ(LOW_IND, JCF)
             HI_UNKZ = CF%UNKZ(HIGH_IND,JCF)
             ! Loop edges to find next face in X2AXIS, X3AXIS plane:
             JAXIS_IEC_LOOP : DO IEC=2,CF%CEDGES(1,JCF)+1
                IEDGE = CF%CEDGES(IEC,JCF)
                JAXIS_IEC_SELECT : SELECT CASE(CF%EDGE_LIST(1,IEDGE))
                CASE(CC_ETYPE_RGGAS) JAXIS_IEC_SELECT ! Cut-faces Regular Gas Edge.
                   LOHI   = CF%EDGE_LIST(2,IEDGE)
                   X2AXIS = CF%EDGE_LIST(3,IEDGE)
                   IF (X2AXIS==IAXIS) THEN ! Edge pointed in X3AXIS=KAXIS direction.
                      X3AXIS = KAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) II=I-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE JAXIS_IEC_LOOP ! Index in CC_RCEDGE
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ I-1
                         IIO=I-1; JJO=J; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-1) ! Y Face in low I
                      ELSE ! Indexes of other Cartesian face @ I+1
                         IIO=I+1; JJO=J; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 1) ! Y Face in high I
                      ENDIF
                      IF(IIO<1 .OR. IIO>M%IBAR) CYCLE JAXIS_IEC_LOOP
                   ELSEIF(X2AXIS==KAXIS) THEN ! Edge pointed in X3AXIS=IAXIS direction.
                      X3AXIS = IAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) KK=K-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE JAXIS_IEC_LOOP
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ K-1
                         IIO=I; JJO=J; KKO=K-1; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-2) ! Y Face in low K
                      ELSE ! Indexes of other Cartesian face @ K+1
                         IIO=I; JJO=J; KKO=K+1; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 2) ! Y Face in high K
                      ENDIF
                      IF(KKO<1 .OR. KKO>M%KBAR) CYCLE JAXIS_IEC_LOOP
                   ENDIF
                   ! Now look for other potential face and link CF%UNKF(JCF):
                   JAXIS_OFACE_SELECT : SELECT CASE (OFACE(1))
                   CASE(CC_FTYPE_RGGAS) JAXIS_OFACE_SELECT
                      OLO_UNKZ = M%CCVAR(IIO  ,JJO  ,KKO  ,CC_UNKZ)
                      OHI_UNKZ = M%CCVAR(IIO  ,JJO+1,KKO  ,CC_UNKZ)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to regular face:
                         IF(M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS) = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)
                         CF%LINK_LEV(JCF) = -1
                         ! IF(I==35 .AND. J== 1 .AND. K==15) &
                         ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, RGGAS J= 1=',I,J,K,IIO,JJO,KKO
                         ! IF(I==35 .AND. J==30 .AND. K==15) &
                         ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, RGGAS J=30=',I,J,K,IIO,JJO,KKO
                         CYCLE JAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_RCGAS) JAXIS_OFACE_SELECT
                      OICF     = OFACE(2)
                      OLO_UNKZ = M%RC_FACE(OICF)%UNKZ(LOW_IND)
                      OHI_UNKZ = M%RC_FACE(OICF)%UNKZ(HIGH_IND)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to RC face:
                         IF(M%RC_FACE(OICF)%UNKF<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%RC_FACE(OICF)%UNKF = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%RC_FACE(OICF)%UNKF
                         CF%LINK_LEV(JCF) = -1
                         ! IF(I==35 .AND. J== 1 .AND. K==15) &
                         ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, RCGAS J= 1=',I,J,K,IIO,JJO,KKO
                         ! IF(I==35 .AND. J==30 .AND. K==15) &
                         ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, RCGAS J=30=',I,J,K,IIO,JJO,KKO
                         CYCLE JAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_CFGAS) JAXIS_OFACE_SELECT
                      OICF     = OFACE(2); OJCF     = OFACE(3)
                      OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                      OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                          IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                             CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                             CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                             ! IF(I==35 .AND. J== 1 .AND. K==15) &
                             ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS1 J= 1=',I,J,K,IIO,JJO,KKO
                             ! IF(I==35 .AND. J==30 .AND. K==15) &
                             ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS1 J=30=',I,J,K,IIO,JJO,KKO
                             CYCLE JAXIS_JCF_LOOP
                          ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                             CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                             ! IF(I==35 .AND. J== 1 .AND. K==15) &
                             ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS2 J= 1=',I,J,K,IIO,JJO,KKO
                             ! IF(I==35 .AND. J==30 .AND. K==15) &
                             ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS2 J=30=',I,J,K,IIO,JJO,KKO
                             CYCLE JAXIS_JCF_LOOP
                          ENDIF
                      ENDIF
                   END SELECT JAXIS_OFACE_SELECT

                CASE(CC_ETYPE_CFGAS) JAXIS_IEC_SELECT ! Gas cut-edge.
                   IECE = CF%EDGE_LIST(2,IEDGE);  JECE = CF%EDGE_LIST(3,IEDGE)
                   II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                   X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                   IF (X3AXIS==KAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = IAXIS
                      IF(II==I-1) THEN ! LOWER cut-face in the I direction.
                         IIO = I-1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-1, JECE) ! Y Face in low I
                      ELSEIF(II==I) THEN ! UPPER cut-face in the I direction.
                         IIO = I+1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 1, JECE) ! Y Face in high I
                      ENDIF
                      IF(IIO<1 .OR. IIO>M%IBAR) CYCLE JAXIS_IEC_LOOP
                   ELSEIF (X3AXIS==IAXIS) THEN ! X2AXIS = KAXIS
                      IF(KK==K-1) THEN ! LOWER cut-face in the K direction.
                         IIO = I; JJO = J; KKO = K-1 ! Indexes of other Cartesian face @ K-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-2, JECE) ! Y Face in low K
                      ELSEIF(KK==K) THEN ! UPPER cut-face in the K direction.
                         IIO = I; JJO = J; KKO = K+1 ! Indexes of other Cartesian face @ K+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 2, JECE) ! Y Face in high K
                      ENDIF
                      IF(KKO<1 .OR. KKO>M%KBAR) CYCLE JAXIS_IEC_LOOP
                   ENDIF
                   OICF     = OFACE(2); OJCF     = OFACE(3)
                   IF(OICF<1 .OR. OJCF<1) CYCLE JAXIS_IEC_LOOP
                   OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                   OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                   IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                       IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                          CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                          CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                          ! IF(I==35 .AND. J== 1 .AND. K==15) &
                          ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS1 J= 1=',I,J,K,IIO,JJO,KKO
                          ! IF(I==35 .AND. J==30 .AND. K==15) &
                          ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS1 J=30=',I,J,K,IIO,JJO,KKO
                          CYCLE JAXIS_JCF_LOOP
                       ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                          M%NUNK_F = M%NUNK_F + 1
                          CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                          CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                          ! IF(I==35 .AND. J== 1 .AND. K==15) &
                          ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS2 J= 1=',I,J,K,IIO,JJO,KKO
                          ! IF(I==35 .AND. J==30 .AND. K==15) &
                          ! WRITE(LU_ERR,*) X3AXIS,'ETYPE_RGGAS, CFGAS2 J=30=',I,J,K,IIO,JJO,KKO
                          CYCLE JAXIS_JCF_LOOP
                       ENDIF
                   ENDIF
                END SELECT JAXIS_IEC_SELECT
             ENDDO JAXIS_IEC_LOOP
          ENDIF
       ENDDO JAXIS_JCF_LOOP

   CASE(KAXIS)
       ACRT = M%DY(J)*M%DX(I)
       KAXIS_JCF_LOOP : DO JCF=1,CF%NFACE
          IF(CF%UNKF(JCF)<1) THEN ! .AND. CF%AREA(JCF) < CCVOL_LINK*ACRT+TWENTY_EPSILON_EB) THEN
             LO_UNKZ = CF%UNKZ(LOW_IND, JCF)
             HI_UNKZ = CF%UNKZ(HIGH_IND,JCF)
             ! Loop edges to find next face in X2AXIS, X3AXIS plane:
             KAXIS_IEC_LOOP : DO IEC=2,CF%CEDGES(1,JCF)+1
                IEDGE = CF%CEDGES(IEC,JCF)
                KAXIS_IEC_SELECT : SELECT CASE(CF%EDGE_LIST(1,IEDGE))
                CASE(CC_ETYPE_RGGAS) KAXIS_IEC_SELECT ! Cut-faces Regular Gas Edge.
                   LOHI   = CF%EDGE_LIST(2,IEDGE)
                   X2AXIS = CF%EDGE_LIST(3,IEDGE)
                   IF (X2AXIS==IAXIS) THEN ! Edge pointed in X3AXIS=JAXIS direction.
                      X3AXIS = JAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) II=I-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE KAXIS_IEC_LOOP ! Index in CC_RCEDGE
                      !WRITE(LU_ERR,*) 'EDGE=',IERC,II,JJ,KK,X3AXIS,', CF=',I,J,K,X1AXIS
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ I-1
                         IIO=I-1; JJO=J; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-2) ! Z Face in low I
                      ELSE ! Indexes of other Cartesian face @ I+1
                         IIO=I+1; JJO=J; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 2) ! Z Face in high I
                      ENDIF
                      IF(IIO<1 .OR. IIO>M%IBAR) CYCLE KAXIS_IEC_LOOP
                   ELSEIF(X2AXIS==JAXIS) THEN ! Edge pointed in X3AXIS=IAXIS direction.
                      X3AXIS = IAXIS; II=I; JJ=J; KK=K; IF(LOHI==LOW_IND) JJ=J-1
                      IERC   = M%ECVAR(II,JJ,KK,CC_IDCE,X3AXIS); IF(IERC==0) CYCLE KAXIS_IEC_LOOP
                      IF(LOHI==LOW_IND) THEN ! Indexes of other Cartesian face @ J-1
                         IIO=I; JJO=J-1; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3,-1) ! Z Face in low J
                      ELSE ! Indexes of other Cartesian face @ J+1
                         IIO=I; JJO=J+1; KKO=K; OFACE(1:3) = M%CC_RCEDGE(IERC)%FACE_LIST(1:3, 1) ! Z Face in high J
                      ENDIF
                      IF(JJO<1 .OR. JJO>M%JBAR) CYCLE KAXIS_IEC_LOOP
                   ENDIF
                   ! Now look for other potential face and link CF%UNKF(JCF):
                   KAXIS_OFACE_SELECT : SELECT CASE (OFACE(1))
                   CASE(CC_FTYPE_RGGAS) KAXIS_OFACE_SELECT
                      OLO_UNKZ = M%CCVAR(IIO  ,JJO  ,KKO  ,CC_UNKZ)
                      OHI_UNKZ = M%CCVAR(IIO  ,JJO  ,KKO+1,CC_UNKZ)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to regular face:
                         IF(M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS) = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)
                         CF%LINK_LEV(JCF) = -1
                         CYCLE KAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_RCGAS) KAXIS_OFACE_SELECT
                      OICF     = OFACE(2)
                      OLO_UNKZ = M%RC_FACE(OICF)%UNKZ(LOW_IND)
                      OHI_UNKZ = M%RC_FACE(OICF)%UNKZ(HIGH_IND)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to RC face:
                         IF(M%RC_FACE(OICF)%UNKF<1) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             M%RC_FACE(OICF)%UNKF = M%NUNK_F
                         ENDIF
                         CF%UNKF(JCF)     = M%RC_FACE(OICF)%UNKF
                         CF%LINK_LEV(JCF) = -1
                         CYCLE KAXIS_JCF_LOOP
                      ENDIF
                   CASE(CC_FTYPE_CFGAS) KAXIS_OFACE_SELECT
                      OICF     = OFACE(2); OJCF     = OFACE(3)
                      OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                      OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                      IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                          IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                             CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                             CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                             CYCLE KAXIS_JCF_LOOP
                          ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                             M%NUNK_F = M%NUNK_F + 1
                             CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                             CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                             CYCLE KAXIS_JCF_LOOP
                          ENDIF
                      ENDIF
                   END SELECT KAXIS_OFACE_SELECT

                CASE(CC_ETYPE_CFGAS) KAXIS_IEC_SELECT ! Gas cut-edge.
                   IECE = CF%EDGE_LIST(2,IEDGE);  JECE = CF%EDGE_LIST(3,IEDGE)
                   II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                   X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                   IF (X3AXIS==JAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = IAXIS
                      IF(II==I-1) THEN ! LOWER cut-face in the I direction.
                         IIO = I-1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-2, JECE) ! Z Face in low I
                      ELSEIF(II==I) THEN ! UPPER cut-face in the I direction.
                         IIO = I+1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 2, JECE) ! Z Face in high I
                      ENDIF
                      IF(IIO<1 .OR. IIO>M%IBAR) CYCLE KAXIS_IEC_LOOP
                   ELSEIF (X3AXIS==IAXIS) THEN ! X2AXIS = JAXIS
                      IF(JJ==J-1) THEN ! LOWER cut-face in the J direction.
                         IIO = I; JJO = J-1; KKO = K ! Indexes of other Cartesian face @ J-1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2,-1, JECE) ! Z Face in low J
                      ELSEIF(JJ==J) THEN ! UPPER cut-face in the J direction.
                         IIO = I; JJO = J+1; KKO = K ! Indexes of other Cartesian face @ J+1
                         OFACE(2:3) = M%CUT_EDGE(IECE)%FACE_LIST(1:2, 1, JECE) ! Z Face in high J
                      ENDIF
                      IF(JJO<1 .OR. JJO>M%JBAR) CYCLE KAXIS_IEC_LOOP
                   ENDIF
                   OICF     = OFACE(2); OJCF     = OFACE(3)
                   IF(OICF<1 .OR. OJCF<1) CYCLE KAXIS_IEC_LOOP
                   OLO_UNKZ = M%CUT_FACE(OICF)%UNKZ(LOW_IND,OJCF)
                   OHI_UNKZ = M%CUT_FACE(OICF)%UNKZ(HIGH_IND,OJCF)
                   IF(LO_UNKZ==OLO_UNKZ .OR. HI_UNKZ==OHI_UNKZ .OR. ALL_FLG) THEN ! Link cut-face to other cut-face:
                       IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
                          CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
                          CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
                          CYCLE KAXIS_JCF_LOOP
                       ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
                          M%NUNK_F = M%NUNK_F + 1
                          CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
                          CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
                          CYCLE KAXIS_JCF_LOOP
                       ENDIF
                   ENDIF
                END SELECT KAXIS_IEC_SELECT
             ENDDO KAXIS_IEC_LOOP
          ENDIF
       ENDDO KAXIS_JCF_LOOP
   END SELECT

ENDDO ICF_LOOP_2

ELSE

ICF_LOOP_3 : DO ICF=1,M%N_CUTFACE_MESH
   CF => M%CUT_FACE(ICF); IF(CF%STATUS/=CC_GASPHASE) CYCLE ICF_LOOP_3
   I  = CF%IJK(IAXIS); J = CF%IJK(JAXIS); K = CF%IJK(KAXIS); X1AXIS = CF%IJK(KAXIS+1)
   SELECT CASE(X1AXIS)
   CASE(IAXIS)

      ACRT = M%DY(J)*M%DZ(K)
      IAXIS_JCF_LOOP_2 : DO JCF=1,CF%NFACE
         IF(CF%UNKF(JCF)<1) THEN ! Loop edges to find next face in X2AXIS, X3AXIS plane:
            COUNT=0; IIO=0; JJO=0; KKO=0
            IAXIS_IEC_LOOP_2 : DO IEC=2,CF%CEDGES(1,JCF)+1
               IEDGE = CF%CEDGES(IEC,JCF)
               SELECT CASE(CF%EDGE_LIST(1,IEDGE))
               CASE(CC_ETYPE_RGGAS) ! Cut-faces Regular Gas Edge.
                  LOHI   = CF%EDGE_LIST(2,IEDGE)
                  X2AXIS = CF%EDGE_LIST(3,IEDGE)
                  IF (X2AXIS==JAXIS) THEN ! Edge pointed in X3AXIS=KAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I; JJO=J-1; KKO=K ! Indexes of other Cartesian face @ J-1
                     ELSE;                   IIO=I; JJO=J+1; KKO=K ! Indexes of other Cartesian face @ J+1
                     ENDIF; IF(JJO<1 .OR. JJO>M%JBAR) CYCLE IAXIS_IEC_LOOP_2
                  ELSEIF(X2AXIS==KAXIS) THEN ! Edge pointed in X3AXIS=JAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I; JJO=J; KKO=K-1 ! Indexes of other Cartesian face @ K-1
                     ELSE;                   IIO=I; JJO=J; KKO=K+1 ! Indexes of other Cartesian face @ K+1
                     ENDIF; IF(KKO<1 .OR. KKO>M%KBAR) CYCLE IAXIS_IEC_LOOP_2
                  ENDIF
               CASE(CC_ETYPE_CFGAS) ! Gas cut-edge.
                  IECE = CF%EDGE_LIST(2,IEDGE)
                  II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                  X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                  IF (X3AXIS==KAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = JAXIS
                     IF(JJ==J-1) THEN;   IIO = I; JJO = J-1; KKO = K ! X Face in low J
                     ELSEIF(JJ==J) THEN; IIO = I; JJO = J+1; KKO = K ! X Face in high J
                     ENDIF; IF(JJO<1 .OR. JJO>M%JBAR) CYCLE IAXIS_IEC_LOOP_2
                  ELSEIF (X3AXIS==JAXIS) THEN ! X2AXIS = KAXIS
                     IF(KK==K-1) THEN;   IIO = I; JJO = J; KKO = K-1 ! X Face in low K
                     ELSEIF(KK==K) THEN; IIO = I; JJO = J; KKO = K+1 ! X Face in high K
                     ENDIF; IF(KKO<1 .OR. KKO>M%KBAR) CYCLE IAXIS_IEC_LOOP_2
                  ENDIF
               CASE DEFAULT; CYCLE IAXIS_IEC_LOOP_2
               END SELECT
               IF(COUNT+1>SIZE_FACE) THEN
                  ALLOCATE(FACELAUX(4,SIZE_FACE+20),FACEARAUX(SIZE_FACE+20))
                  FACELAUX(1:4,1:SIZE_FACE) = FACE_LIST(1:4,1:SIZE_FACE)
                  FACEARAUX(1:SIZE_FACE)    = FACE_AREA(1:SIZE_FACE)
                  CALL MOVE_ALLOC(FROM=FACELAUX,TO=FACE_LIST)
                  CALL MOVE_ALLOC(FROM=FACEARAUX,TO=FACE_AREA)
                  SIZE_FACE = SIZE_FACE+20
               ENDIF
               IF(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS)>0) THEN ! Underlying cut-face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_CFGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = M%CUT_FACE(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS))%ALPHA_CF*ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_IDRC,X1AXIS)>0) THEN ! RC face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RCGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_CGSC,X1AXIS)==CC_GASPHASE) THEN ! Regular face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RGGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ENDIF
            ENDDO  IAXIS_IEC_LOOP_2
            IF(COUNT==0) CYCLE IAXIS_JCF_LOOP_2
            CALL SET_UNKF_CF
         ENDIF
      ENDDO IAXIS_JCF_LOOP_2

   CASE(JAXIS)

      IF (TWO_D) CYCLE ICF_LOOP_3
      ACRT = M%DX(I)*M%DZ(K)
      JAXIS_JCF_LOOP_2 : DO JCF=1,CF%NFACE
         IF(CF%UNKF(JCF)<1) THEN
            COUNT=0; IIO=0; JJO=0; KKO=0
            JAXIS_IEC_LOOP_2 : DO IEC=2,CF%CEDGES(1,JCF)+1
               IEDGE = CF%CEDGES(IEC,JCF)
               SELECT CASE(CF%EDGE_LIST(1,IEDGE))
               CASE(CC_ETYPE_RGGAS) ! Cut-faces Regular Gas Edge.
                  LOHI   = CF%EDGE_LIST(2,IEDGE)
                  X2AXIS = CF%EDGE_LIST(3,IEDGE)
                  IF (X2AXIS==IAXIS) THEN ! Edge pointed in X3AXIS=KAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I-1; JJO=J; KKO=K ! Indexes of other Cartesian face @ I-1
                     ELSE;                   IIO=I+1; JJO=J; KKO=K ! Indexes of other Cartesian face @ I+1
                     ENDIF; IF(IIO<1 .OR. IIO>M%IBAR) CYCLE JAXIS_IEC_LOOP_2
                  ELSEIF(X2AXIS==KAXIS) THEN ! Edge pointed in X3AXIS=IAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I; JJO=J; KKO=K-1 ! Indexes of other Cartesian face @ K-1
                     ELSE;                   IIO=I; JJO=J; KKO=K+1 ! Indexes of other Cartesian face @ K+1
                     ENDIF; IF(KKO<1 .OR. KKO>M%KBAR) CYCLE JAXIS_IEC_LOOP_2
                  ENDIF
               CASE(CC_ETYPE_CFGAS) ! Gas cut-edge.
                  IECE = CF%EDGE_LIST(2,IEDGE)
                  II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                  X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                  IF (X3AXIS==KAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = IAXIS
                     IF(II==I-1) THEN;   IIO = I-1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I-1
                     ELSEIF(II==I) THEN; IIO = I+1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I+1
                     ENDIF; IF(IIO<1 .OR. IIO>M%IBAR) CYCLE JAXIS_IEC_LOOP_2
                  ELSEIF (X3AXIS==IAXIS) THEN ! X2AXIS = KAXIS
                     IF(KK==K-1) THEN;   IIO = I; JJO = J; KKO = K-1 ! Indexes of other Cartesian face @ K-1
                     ELSEIF(KK==K) THEN; IIO = I; JJO = J; KKO = K+1 ! Indexes of other Cartesian face @ K+1
                     ENDIF; IF(KKO<1 .OR. KKO>M%KBAR) CYCLE JAXIS_IEC_LOOP_2
                  ENDIF
               CASE DEFAULT; CYCLE JAXIS_IEC_LOOP_2
               END SELECT
               IF(COUNT+1>SIZE_FACE) THEN
                  ALLOCATE(FACELAUX(4,SIZE_FACE+20),FACEARAUX(SIZE_FACE+20))
                  FACELAUX(1:4,1:SIZE_FACE) = FACE_LIST(1:4,1:SIZE_FACE)
                  FACEARAUX(1:SIZE_FACE)    = FACE_AREA(1:SIZE_FACE)
                  CALL MOVE_ALLOC(FROM=FACELAUX,TO=FACE_LIST)
                  CALL MOVE_ALLOC(FROM=FACEARAUX,TO=FACE_AREA)
                  SIZE_FACE = SIZE_FACE+20
               ENDIF
               IF(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS)>0) THEN ! Underlying cut-face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_CFGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = M%CUT_FACE(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS))%ALPHA_CF*ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_IDRC,X1AXIS)>0) THEN ! RC face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RCGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_CGSC,X1AXIS)==CC_GASPHASE) THEN ! Regular face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RGGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ENDIF
            ENDDO  JAXIS_IEC_LOOP_2
            IF(COUNT==0) CYCLE JAXIS_JCF_LOOP_2
            CALL SET_UNKF_CF
         ENDIF
      ENDDO JAXIS_JCF_LOOP_2

   CASE(KAXIS)

      ACRT = M%DY(J)*M%DX(I)
      KAXIS_JCF_LOOP_2 : DO JCF=1,CF%NFACE
         IF(CF%UNKF(JCF)<1) THEN
            COUNT=0; IIO=0; JJO=0; KKO=0
            KAXIS_IEC_LOOP_2 : DO IEC=2,CF%CEDGES(1,JCF)+1
               IEDGE = CF%CEDGES(IEC,JCF)
               SELECT CASE(CF%EDGE_LIST(1,IEDGE))
               CASE(CC_ETYPE_RGGAS) ! Cut-faces Regular Gas Edge.
                  LOHI   = CF%EDGE_LIST(2,IEDGE)
                  X2AXIS = CF%EDGE_LIST(3,IEDGE)
                  IF (X2AXIS==IAXIS) THEN ! Edge pointed in X3AXIS=JAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I-1; JJO=J; KKO=K ! Indexes of other Cartesian face @ I-1
                     ELSE;                   IIO=I+1; JJO=J; KKO=K ! Indexes of other Cartesian face @ I+1
                     ENDIF; IF(IIO<1 .OR. IIO>M%IBAR) CYCLE KAXIS_IEC_LOOP_2
                  ELSEIF(X2AXIS==JAXIS) THEN ! Edge pointed in X3AXIS=IAXIS direction.
                     IF(LOHI==LOW_IND) THEN; IIO=I; JJO=J-1; KKO=K ! Indexes of other Cartesian face @ J-1
                     ELSE;                   IIO=I; JJO=J+1; KKO=K ! Indexes of other Cartesian face @ J+1
                     ENDIF; IF(JJO<1 .OR. JJO>M%JBAR) CYCLE KAXIS_IEC_LOOP_2
                  ENDIF
               CASE(CC_ETYPE_CFGAS) ! Gas cut-edge.
                  IECE = CF%EDGE_LIST(2,IEDGE)
                  II=M%CUT_EDGE(IECE)%IJK(IAXIS); JJ=M%CUT_EDGE(IECE)%IJK(JAXIS); KK=M%CUT_EDGE(IECE)%IJK(KAXIS);
                  X3AXIS = M%CUT_EDGE(IECE)%IJK(KAXIS+1);
                  IF (X3AXIS==JAXIS) THEN ! X3AXIS==Edge Axis, X2AXIS = IAXIS
                     IF(II==I-1) THEN;   IIO = I-1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I-1
                     ELSEIF(II==I) THEN; IIO = I+1; JJO = J; KKO = K ! Indexes of other Cartesian face @ I+1
                     ENDIF; IF(IIO<1 .OR. IIO>M%IBAR) CYCLE KAXIS_IEC_LOOP_2
                  ELSEIF (X3AXIS==IAXIS) THEN ! X2AXIS = JAXIS
                     IF(JJ==J-1) THEN;   IIO = I; JJO = J-1; KKO = K ! Indexes of other Cartesian face @ J-1
                     ELSEIF(JJ==J) THEN; IIO = I; JJO = J+1; KKO = K ! Indexes of other Cartesian face @ J+1
                     ENDIF; IF(JJO<1 .OR. JJO>M%JBAR) CYCLE KAXIS_IEC_LOOP_2
                  ENDIF
               CASE DEFAULT; CYCLE KAXIS_IEC_LOOP_2
               END SELECT
               IF(COUNT+1>SIZE_FACE) THEN
                  ALLOCATE(FACELAUX(4,SIZE_FACE+20),FACEARAUX(SIZE_FACE+20))
                  FACELAUX(1:4,1:SIZE_FACE) = FACE_LIST(1:4,1:SIZE_FACE)
                  FACEARAUX(1:SIZE_FACE)    = FACE_AREA(1:SIZE_FACE)
                  CALL MOVE_ALLOC(FROM=FACELAUX,TO=FACE_LIST)
                  CALL MOVE_ALLOC(FROM=FACEARAUX,TO=FACE_AREA)
                  SIZE_FACE = SIZE_FACE+20
               ENDIF
               IF(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS)>0) THEN ! Underlying cut-face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_CFGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = M%CUT_FACE(M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS))%ALPHA_CF*ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_IDRC,X1AXIS)>0) THEN ! RC face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RCGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ELSEIF(M%FCVAR(IIO,JJO,KKO,CC_CGSC,X1AXIS)==CC_GASPHASE) THEN ! Regular face.
                  COUNT=COUNT+1
                  FACE_LIST(1:4,COUNT) = (/ CC_FTYPE_RGGAS,IIO,JJO,KKO/)
                  FACE_AREA(COUNT)     = ACRT
               ENDIF
            ENDDO  KAXIS_IEC_LOOP_2
            IF(COUNT==0) CYCLE KAXIS_JCF_LOOP_2
            CALL SET_UNKF_CF
         ENDIF
      ENDDO  KAXIS_JCF_LOOP_2

   END SELECT
ENDDO ICF_LOOP_3

ENDIF

END SUBROUTINE LINK_FACES_TO_NEIGHBORS

SUBROUTINE SET_UNKF_CF

! Now define face to link to:
ILOC = MAXLOC(FACE_AREA(1:COUNT),DIM=1)
IIO  = FACE_LIST(2,ILOC); JJO = FACE_LIST(3,ILOC); KKO = FACE_LIST(4,ILOC)
SELECT CASE (FACE_LIST(1,ILOC))
CASE(CC_FTYPE_RGGAS)
   IF(M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)<1) THEN
         M%NUNK_F = M%NUNK_F + 1
         M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS) = M%NUNK_F
   ENDIF
   CF%UNKF(JCF)     = M%FCVAR(IIO,JJO,KKO,CC_UNKF,X1AXIS)
   CF%LINK_LEV(JCF) = -1
CASE(CC_FTYPE_RCGAS)
   OICF = M%FCVAR(IIO,JJO,KKO,CC_IDRC,X1AXIS)
   IF(M%RC_FACE(OICF)%UNKF<1) THEN
      M%NUNK_F = M%NUNK_F + 1
      M%RC_FACE(OICF)%UNKF = M%NUNK_F
   ENDIF
   CF%UNKF(JCF)     = M%RC_FACE(OICF)%UNKF
   CF%LINK_LEV(JCF) = -1
CASE(CC_FTYPE_CFGAS)
   OICF = M%FCVAR(IIO,JJO,KKO,CC_IDCF,X1AXIS); OJCF = 1
   IF(M%CUT_FACE(OICF)%UNKF(OJCF)>0) THEN
      CF%UNKF(JCF)     = M%CUT_FACE(OICF)%UNKF(OJCF)
      CF%LINK_LEV(JCF) = M%CUT_FACE(OICF)%LINK_LEV(OJCF) - 1
   ELSEIF(CF%AREA(JCF)+M%CUT_FACE(OICF)%AREA(OJCF) > CCVOL_LINK * ACRT) THEN
      M%NUNK_F = M%NUNK_F + 1
      CF%UNKF(JCF) = M%NUNK_F; M%CUT_FACE(OICF)%UNKF(OJCF) = M%NUNK_F
      CF%LINK_LEV(JCF) = -1; M%CUT_FACE(OICF)%LINK_LEV(OJCF) = -1
   ENDIF
END SELECT
RETURN
END SUBROUTINE SET_UNKF_CF

END SUBROUTINE GET_LINKED_FACE_INDEXES_F


END MODULE CC_INIT
