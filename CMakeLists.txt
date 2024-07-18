cmake_minimum_required(VERSION 3.16)
project(loop-fds LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

set(Hedgehog_INCLUDE_DIR Lib/hedgehog/)

# Need to defined variable Hedgehog_INCLUDE_DIR to the hedgehog path hedgehog
find_package(Hedgehog REQUIRED)

# OpenBLAS
# find_package(OpenBLAS QUIET)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${Hedgehog_CXX_FLAGS})
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")

if (NOT MSVC)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -Wextra -Werror -Wuninitialized -pedantic -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Werror -Wuninitialized -pedantic -g")
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Werror -Wuninitialized -pedantic -O3")
endif (NOT MSVC)

# files
set(cpp_files
        Source/CPP/main.cpp
)

set(fortran_files
        # Source/Fortran/main.f90
)

# executable
add_executable(fds ${cpp_files} ${fortran_files})