cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

project(
    fds
    VERSION 6.9.1
    LANGUAGES Fortran C
)
enable_language(Fortran)

# Enable C++ if building with Hedgehog
if(USE_HEDGEHOG)
    enable_language(CXX)
endif()

option(USE_HYPRE           "Use the hypre library"                    ON)
option(USE_SYSTEM_HYPRE    "Use the hypre library from the system"    OFF)

option(USE_SUNDIALS        "Use the sundials library"                 ON)
option(USE_SYSTEM_SUNDIALS "Use the sundials library from the system" OFF)

option(USE_OPENMP          "Use OpenMP"                               ON)
option(USE_HEDGEHOG        "Build fds_hh target with Hedgehog dataflow" OFF)

# The existing FDS makefile links MKL statically, but this is optional and can
# be changed to dynamic
set(MKL_LINK static CACHE STRING "Linking method for MKL (static or dynamic)")
set_property(CACHE MKL_LINK PROPERTY STRINGS static dynamic)

# Fortran source files shared between fds and fds_hh targets (everything except main.f90)
set(FDS_FORTRAN_SOURCES
    ${CMAKE_SOURCE_DIR}/Source/prec.f90
    ${CMAKE_SOURCE_DIR}/Source/cons.f90
    ${CMAKE_SOURCE_DIR}/Source/chem.f90
    ${CMAKE_SOURCE_DIR}/Source/prop.f90
    ${CMAKE_SOURCE_DIR}/Source/devc.f90
    ${CMAKE_SOURCE_DIR}/Source/type.f90
    ${CMAKE_SOURCE_DIR}/Source/data.f90
    ${CMAKE_SOURCE_DIR}/Source/mesh.f90
    ${CMAKE_SOURCE_DIR}/Source/func.f90
    ${CMAKE_SOURCE_DIR}/Source/gsmv.f90
    ${CMAKE_SOURCE_DIR}/Source/smvv.f90
    ${CMAKE_SOURCE_DIR}/Source/rcal.f90
    ${CMAKE_SOURCE_DIR}/Source/turb.f90
    ${CMAKE_SOURCE_DIR}/Source/soot.f90
    ${CMAKE_SOURCE_DIR}/Source/pois.f90
    ${CMAKE_SOURCE_DIR}/Source/geom.f90
    ${CMAKE_SOURCE_DIR}/Source/ccib.f90
    ${CMAKE_SOURCE_DIR}/Source/radi.f90
    ${CMAKE_SOURCE_DIR}/Source/part.f90
    ${CMAKE_SOURCE_DIR}/Source/vege.f90
    ${CMAKE_SOURCE_DIR}/Source/ctrl.f90
    ${CMAKE_SOURCE_DIR}/Source/hvac.f90
    ${CMAKE_SOURCE_DIR}/Source/mass.f90
    ${CMAKE_SOURCE_DIR}/Source/imkl.f90
    ${CMAKE_SOURCE_DIR}/Source/wall.f90
    ${CMAKE_SOURCE_DIR}/Source/fire.f90
    ${CMAKE_SOURCE_DIR}/Source/velo_kernels.f90
    ${CMAKE_SOURCE_DIR}/Source/velo.f90
    ${CMAKE_SOURCE_DIR}/Source/pres.f90
    ${CMAKE_SOURCE_DIR}/Source/init.f90
    ${CMAKE_SOURCE_DIR}/Source/dump.f90
    ${CMAKE_SOURCE_DIR}/Source/read.f90
    ${CMAKE_SOURCE_DIR}/Source/divg.f90
)

add_executable(fds
    ${CMAKE_SOURCE_DIR}/Source/main.f90
    ${FDS_FORTRAN_SOURCES}
)
target_include_directories(fds PRIVATE .)

# Get various properties about the time and git revision. These can be
# overridden which is important for building in situations where we don't have
# the git repo.
if (NOT(WIN32) AND NOT(BUILD_DATE_XLF))
    string(TIMESTAMP BUILD_DATE_XLF "%b %d\\, %Y  %H:%M:%S")
endif()
if (NOT(BUILD_DATE))
    string(TIMESTAMP BUILD_DATE "%b %d, %Y  %H:%M:%S ")
endif()
if (NOT(GIT_DATE))
    execute_process(COMMAND git log -1 --format=%cd          OUTPUT_VARIABLE GIT_DATE       OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if (NOT(GIT_BRANCH))
    execute_process(COMMAND git rev-parse --abbrev-ref HEAD  OUTPUT_VARIABLE GIT_BRANCH     OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if (NOT(GIT_HASH))
    execute_process(COMMAND git describe --long --abbrev=7   OUTPUT_VARIABLE GIT_HASH       OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if (NOT(GIT_DIRTY))
    execute_process(COMMAND git diff --shortstat Source/*.f90 OUTPUT_VARIABLE GIT_DIRTY)
endif()
if (GIT_DIRTY STREQUAL "")
    set(GIT_STAT "")
else()
    set(GIT_STAT "-dirty")
endif()

target_compile_definitions(fds PRIVATE BUILDDATE_PP="${BUILD_DATE}")
target_compile_definitions(fds PRIVATE GITHASH_PP="${GIT_HASH}${GIT_STAT}-${GIT_BRANCH}")
target_compile_definitions(fds PRIVATE GITDATE_PP="${GIT_DATE}")

# Use MPI (mandatory)
# When building with Hedgehog, we enable C++ for Hedgehog itself, not for MPI
# However, HYPRE and SUNDIALS detect C++ is enabled and try to find MPI::MPI_CXX
# We find it here to satisfy their build system, but FDS only uses MPI::MPI_Fortran
if(USE_HEDGEHOG)
    find_package(MPI REQUIRED COMPONENTS Fortran C CXX)
else()
    find_package(MPI REQUIRED COMPONENTS Fortran)
endif()
target_link_libraries(fds PRIVATE MPI::MPI_Fortran)

# FDS use lp64 (4-byte integer), MKL often defaults to ilp64 (8-byte integer)
set(MKL_INTERFACE lp64)
# Turn on BLACS
set(ENABLE_BLACS ON)
# We don't need ScaLAPACK
set(ENABLE_SCALAPACK OFF)
# Use MKL if found
find_package(MKL CONFIG)
if (MKL_FOUND)
    message(STATUS "${MKL_IMPORTED_TARGETS}")
    get_target_property(MKLS MKL::MKL INTERFACE_COMPILE_OPTIONS)
    message(STATUS "MKL Compile Options: ${MKLS}")
    get_target_property(MKLS MKL::MKL INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "MKL Include Directories: ${MKLS}")
    target_compile_definitions(fds PRIVATE WITH_MKL)
    target_link_libraries(fds PRIVATE MKL::MKL)
endif()

if (USE_OPENMP)
    # Use OpenMP
    find_package(OpenMP REQUIRED)
    target_link_libraries(fds PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Set compiler flags for various compilers
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fds PRIVATE -cpp -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(fds PRIVATE -fpp )
    # When compiling with Intel on Linux, enable the IFPORT library
    if (LINUX)
        target_compile_definitions(fds PRIVATE USE_IFPORT)
    endif()
    # When compiling with Intel and not on Windows apply -no-wrap-margin.
    if (NOT(WIN32))
        target_compile_options(fds PRIVATE -no-wrap-margin)
    endif()
endif()

# If we are using the old Intel Fortran compiler (ifort) suppress the warning
# that this is an old compiler.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    target_compile_options(fds PRIVATE /Qdiag-disable:10448)
endif()

if(WIN32)
    set(BUILD_SHARED_LIBS OFF)
endif()

# Find or build hypre if the option is set to use it
if(USE_HYPRE)
    # Unless forced to use the system version, download and build hypre
    if(NOT(USE_SYSTEM_HYPRE))
        include(FetchContent)
        # As we are not using the system hypre, we need to choose the version we
        # want
        set(HYPRE_GIT_VERSION "3.0.0" )
        FetchContent_Declare(
            HYPRE
            GIT_REPOSITORY https://github.com/hypre-space/hypre.git
            # Currently we need to refer directly to a commit which contains a
            # patch to handle HYPRE_FMANGLE
            GIT_TAG        6d01417697aa84d118eb1e1eb91315df42575be5 # v${HYPRE_GIT_VERSION}
            SOURCE_SUBDIR  src
            OVERRIDE_FIND_PACKAGE
        )
        # Set some options for building hypre
        if (WIN32)
            set(HYPRE_ENABLE_FMANGLE CAPS)
        endif()
        FetchContent_MakeAvailable(HYPRE)
    endif()

    # Use hypre if found
    find_package(HYPRE)
    if(HYPRE_FOUND)
        if(NOT(HYPRE_VERSION))
            # We aren't using the system version so we'll revert to the git
            # version
            set(HYPRE_VERSION "${HYPRE_GIT_VERSION}")
        endif()
        target_link_libraries(fds PRIVATE HYPRE)
        target_compile_definitions(fds PRIVATE WITH_HYPRE)
        target_compile_definitions(fds PRIVATE HYPRE_PP="${HYPRE_VERSION}")
    endif()
endif()

# Find or build sundials if the option is set to use it
if(USE_SUNDIALS)
    # Unless forced to use the system version, download and build sundials
    if(NOT(USE_SYSTEM_SUNDIALS))
        include(FetchContent)
        # As we are not using the system sundials, we need to choose the version
        # we want
        set(SUNDIALS_GIT_VERSION "7.5.0")
        FetchContent_Declare(
            SUNDIALS
            GIT_REPOSITORY https://github.com/LLNL/sundials.git
            GIT_TAG        v${SUNDIALS_GIT_VERSION}
            OVERRIDE_FIND_PACKAGE
        )
        # Set some options for building sundials
        set(ENABLE_MPI ON CACHE BOOL "" FORCE)
        set(BUILD_FORTRAN_MODULE_INTERFACE ON CACHE BOOL "" FORCE)
        set(EXAMPLES_ENABLE_C OFF CACHE BOOL "" FORCE)
        set(EXAMPLES_ENABLE_CXX OFF CACHE BOOL "" FORCE)
        set(EXAMPLES_ENABLE_F2003 OFF CACHE BOOL "" FORCE)
        set(EXAMPLES_INSTALL OFF CACHE BOOL "" FORCE)
        set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
        set(BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
        # if(HYPRE_FOUND)
        #     set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
        #     FetchContent_MakeAvailable(HYPRE SUNDIALS)
        # else()
            FetchContent_MakeAvailable(SUNDIALS)
        # endif()
    endif()

    # Use sundials if found
    find_package(SUNDIALS)
    if(SUNDIALS_FOUND)
        if(NOT(SUNDIALS_VERSION))
            # We aren't using the system version so we'll revert to the git
            # version
            set(SUNDIALS_VERSION "${SUNDIALS_GIT_VERSION}")
        endif()
        target_link_libraries(fds PRIVATE SUNDIALS::fcvode_mod)
        target_link_libraries(fds PRIVATE SUNDIALS::fnvecserial_mod)
        target_compile_definitions(fds PRIVATE WITH_SUNDIALS)
        target_compile_definitions(fds PRIVATE SUNDIALS_PP="${SUNDIALS_VERSION}")
    endif()
endif()

install(TARGETS fds)

# Hedgehog dataflow integration
if(USE_HEDGEHOG)
    add_subdirectory(Source/hedgehog)
endif()

include(CTest)
enable_testing()
add_test(NAME "FDS Executes"
    COMMAND fds)
