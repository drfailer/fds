# CMakeLists.txt for the fds_hh target (Hedgehog dataflow integration)
#
# Build with:
#   cmake -DUSE_HEDGEHOG=ON -DHEDGEHOG_INCLUDE_DIR=/path/to/hedgehog ..
#   cmake --build . --target fds_hh

# Enable C++ for this subdirectory
enable_language(CXX)

# Require C++20 for Hedgehog
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hedgehog is header-only; user must provide the include path
if(NOT DEFINED HEDGEHOG_INCLUDE_DIR)
    # Try to find it as a sibling directory
    if(EXISTS "${CMAKE_SOURCE_DIR}/../hedgehog/hedgehog")
        set(HEDGEHOG_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../hedgehog")
        message(STATUS "Found Hedgehog at ${HEDGEHOG_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR
            "HEDGEHOG_INCLUDE_DIR not set. Please provide the path to the Hedgehog library, e.g.:\n"
            "  cmake -DHEDGEHOG_INCLUDE_DIR=/path/to/hedgehog ..")
    endif()
endif()

# Hedgehog requires threading and TBB
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# Create the fds_hh executable
add_executable(fds_hh
    # C++ sources
    main_hh.cpp
    # Fortran sources: the driver module and C interface
    fds_driver.f90
    fds_c_interface.f90
    # All FDS Fortran sources (minus main.f90)
    ${FDS_FORTRAN_SOURCES}
)

# Include directories
target_include_directories(fds_hh PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}     # For fds_fortran_interface.h and subdirs
    ${HEDGEHOG_INCLUDE_DIR}         # For hedgehog/hedgehog.h
    ${CMAKE_SOURCE_DIR}             # For the top-level directory
)

# Link threading and TBB (required by Hedgehog)
target_link_libraries(fds_hh PRIVATE Threads::Threads)
target_link_libraries(fds_hh PRIVATE TBB::tbb)

# Link MPI (mandatory for FDS)
target_link_libraries(fds_hh PRIVATE MPI::MPI_Fortran)

# Link MKL if available
if(MKL_FOUND)
    target_compile_definitions(fds_hh PRIVATE WITH_MKL)
    target_link_libraries(fds_hh PRIVATE MKL::MKL)
endif()

# Link hypre if available
if(USE_HYPRE AND HYPRE_FOUND)
    target_link_libraries(fds_hh PRIVATE HYPRE)
    target_compile_definitions(fds_hh PRIVATE WITH_HYPRE)
    target_compile_definitions(fds_hh PRIVATE HYPRE_PP="${HYPRE_VERSION}")
endif()

# Link sundials if available
if(USE_SUNDIALS AND SUNDIALS_FOUND)
    target_link_libraries(fds_hh PRIVATE SUNDIALS::fcvode_mod)
    target_link_libraries(fds_hh PRIVATE SUNDIALS::fnvecserial_mod)
    target_compile_definitions(fds_hh PRIVATE WITH_SUNDIALS)
    target_compile_definitions(fds_hh PRIVATE SUNDIALS_PP="${SUNDIALS_VERSION}")
endif()

# Propagate compile definitions from the main fds target
target_compile_definitions(fds_hh PRIVATE BUILDDATE_PP="${BUILD_DATE}")
target_compile_definitions(fds_hh PRIVATE GITHASH_PP="${GIT_HASH}${GIT_STAT}-${GIT_BRANCH}")
target_compile_definitions(fds_hh PRIVATE GITDATE_PP="${GIT_DATE}")

# Fortran compiler flags (same as the main fds target)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fds_hh PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-cpp -std=f2018 -ffpe-summary=none -fall-intrinsics>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(fds_hh PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fpp>)
    if(LINUX)
        target_compile_definitions(fds_hh PRIVATE USE_IFPORT)
    endif()
    if(NOT WIN32)
        target_compile_options(fds_hh PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-no-wrap-margin>)
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(fds_hh PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:10448>)
endif()

install(TARGETS fds_hh)
