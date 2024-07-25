#ifndef THIN_WALL_H
#define THIN_WALL_H

namespace mesh::wall {

struct ThinWall {
   int BC_INDEX=0;      ///< Index within the array BOUNDARY_COORD
   int OD_INDEX=0;      ///< Index within the array BOUNDARY_ONE_D
   int TD_INDEX=0;      ///< Index within the array BOUNDARY_THR_D
   int B1_INDEX=0;      ///< Index within the array BOUNDARY_PROP1
   int SURF_INDEX=0;    ///< Index of the SURFace conditions
   int BOUNDARY_TYPE=0; ///< Descriptor: SOLID, MIRROR, OPEN, INTERPOLATED, etc
   int OBST_INDEX=0;    ///< Index of the OBSTruction
   int IEC=0;           ///< Orientation index (1=constant x, 2=constant y, 3-constant z)
   int WALL_INDEX_M=0;  ///< Lower adjacent wall cell index
   int WALL_INDEX_P=0;  ///< Upper adjacent wall cell index
   int N_REALS=0;       ///< Number of reals to pack into restart or send/recv buffer
   int N_INTEGERS=0;    ///< Number of integers to pack into restart or send/recv buffer
   int N_LOGICALS=0;    ///< Number of logicals to pack into restart or send/recv buffer
};

} // end namespace mesh::wall

#endif
