cmake_minimum_required(VERSION 3.16)
project(loop-fds LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

set(Hedgehog_INCLUDE_DIR Lib/hedgehog/)

# Need to defined variable Hedgehog_INCLUDE_DIR to the hedgehog path hedgehog
find_package(Hedgehog REQUIRED)
find_package(MPI REQUIRED)

# OpenBLAS
# find_package(OpenBLAS QUIET)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${Hedgehog_CXX_FLAGS})
set(CMAKE_Fortran_MODULE_DIRECTORY modules)
set(CMAKE_Fortran_COMPILER mpifort)
set(CMAKE_Fortran_FLAGS
  "-m64 -O2 -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics -fPIC -fopenmp -cpp"
)

if (NOT MSVC)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -Wall -Werror -pedantic -g -fPIC -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -g -fPIC -fopenmp")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -O3 -fPIC -fopenmp")
endif (NOT MSVC)

# files
set(cpp_files
        Source/CPP/main.cpp
)

set(fortran_files
  Source/Fortran/prec.f90
  Source/Fortran/cons.f90
  Source/Fortran/chem.f90
  Source/Fortran/prop.f90
  Source/Fortran/devc.f90
  Source/Fortran/type.f90
  Source/Fortran/data.f90
  Source/Fortran/mesh.f90
  Source/Fortran/func.f90
  Source/Fortran/gsmv.f90
  Source/Fortran/smvv.f90
  Source/Fortran/rcal.f90
  Source/Fortran/turb.f90
  Source/Fortran/soot.f90
  Source/Fortran/pois.f90
  Source/Fortran/geom.f90
  Source/Fortran/ccib.f90
  Source/Fortran/radi.f90
  Source/Fortran/part.f90
  Source/Fortran/vege.f90
  Source/Fortran/ctrl.f90
  Source/Fortran/hvac.f90
  Source/Fortran/mass.f90
  Source/Fortran/wall.f90
  Source/Fortran/fire.f90
  Source/Fortran/velo.f90
  Source/Fortran/pres.f90
  Source/Fortran/init.f90
  Source/Fortran/dump.f90
  Source/Fortran/read.f90
  Source/Fortran/divg.f90
  Source/Fortran/fds.f90
)

# executable
add_executable(fds ${cpp_files} ${fortran_files})
target_link_libraries(fds ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
