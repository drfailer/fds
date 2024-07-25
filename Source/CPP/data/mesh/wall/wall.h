#ifndef MESH_WALL_WALL_H
#define MESH_WALL_WALL_H

namespace mesh::wall {

struct Wall {
   int BC_INDEX=0;       ///< Index within the array BOUNDARY_COORD
   int OD_INDEX=0;       ///< Index within the array BOUNDARY_ONE_D
   int TD_INDEX=0;       ///< Index within the array BOUNDARY_THR_D
   int B1_INDEX=0;       ///< Index within the array BOUNDARY_PROP1
   int B2_INDEX=0;       ///< Index within the array BOUNDARY_PROP2
   int BR_INDEX=0;       ///< Index within the array BOUNDARY_RADIA
   int SURF_INDEX=0;     ///< Index of the SURFace conditions
   int BOUNDARY_TYPE=0;  ///< Descriptor: SOLID, MIRROR, OPEN, INTERPOLATED, etc
   int OBST_INDEX=0;     ///< Index of the OBSTruction
   int VENT_INDEX=0;     ///< Index of the VENT containing this cell
   int JD11_INDEX=0;
   int JD12_INDEX=0;
   int JD21_INDEX=0;
   int JD22_INDEX=0;
   int CUT_FACE_INDEX=0;
   int N_REALS=0;        ///< Number of reals to pack into restart or send/recv buffer
   int N_INTEGERS=0;     ///< Number of integers to pack into restart or send/recv buffer
   int N_LOGICALS=0;     ///< Number of logicals to pack into restart or send/recv buffer
};

} // end namespace wall

#endif
