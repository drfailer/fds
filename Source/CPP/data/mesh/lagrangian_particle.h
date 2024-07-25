#ifndef MESH_LAGRANGIAN_PARTICLE_H
#define MESH_LAGRANGIAN_PARTICLE_H

namespace mesh {

struct LagrangianParticle {
  int BC_INDEX=0;             ///< Coordinate variables
  int OD_INDEX=0;             ///< Variables devoted to 1-D heat conduction in depth
  int B1_INDEX=0;             ///< Variables devoted to surface properties
  int B2_INDEX=0;             ///< Variables devoted to surface properties
  int BR_INDEX=0;             ///< Variables devoted to radiation intensities
  int TAG;                    ///< Unique integer identifier for the particle
  int CLASS_INDEX=0;          ///< LAGRANGIAN_PARTICLE_CLASS of particle
  int INITIALIZATION_INDEX=0; ///< Index for INIT that placed the particle
  int ORIENTATION_INDEX=0;    ///< Index in the array of all ORIENTATIONs
  int WALL_INDEX=0;           ///< If liquid droplet has stuck to a wall, this is the WALL cell index
  int DUCT_INDEX=0;           ///< Index of duct
  int INIT_INDEX=0;           ///< Index of INIT line
  int DUCT_CELL_INDEX=0;      ///< Index of duct cell
  int CFACE_INDEX=0;          ///< Index of immersed boundary CFACE that the droplet has attached to
  int PROP_INDEX=0;           ///< Index of the PROPERTY_TYPE assigned to the particle
  int N_REALS=0;              ///< Number of reals to pack into restart or send/recv buffer
  int N_INTEGERS=0;           ///< Number of integers to pack into restart or send/recv buffer
  int N_LOGICALS=0;           ///< Number of logicals to pack into restart or send/recv buffer

  bool SHOW=false;          ///< Show the particle in Smokeview
  bool SPLAT=false;         ///< The liquid droplet has hit a solid
  bool EMBER=false;         ///< The particle can break away and become a burning ember
  bool PATH_PARTICLE=false; ///< Indicator of a particle with a specified path

  double U=0;        ///< \f$ x \f$ velocity component of particle (m/s)
  double V=0;        ///< \f$ y \f$ velocity component of particle (m/s)
  double W=0;        ///< \f$ z \f$ velocity component of particle (m/s)
  double PWT=1;      ///< Weight factor of particle; i.e. the number of real particles it represents
  double RE=0;       ///< Reynolds number based on particle diameter
  double MASS=0;     ///< Particle mass (kg)
  double T_INSERT=0; ///< Time when particle was inserted (s)
  double DX=1;       ///< Length factor used in POROUS_DRAG calculation (m)
  double DY=1;       ///< Length factor used in POROUS_DRAG calculation (m)
  double DZ=1;       ///< Length factor used in POROUS_DRAG calculation (m)
  double LENGTH=-1;  ///< Cylinder or plate length used for POROUS_DRAG or SCREEN_DRAG (m)
  double C_DRAG=0;   ///< Drag coefficient
  double RADIUS=0;   ///< Radius (m)
  double ACCEL_X=0;  ///< Acceleration in x direction (m/s2)
  double ACCEL_Y=0;  ///< Acceleration in y direction (m/s2)
  double ACCEL_Z=0;  ///< Acceleration in z direction (m/s2)
  double RVC=-1;     ///< Reciprocal of cell volume containing particle (1/m3)
};

} // end namespace mesh

#endif
