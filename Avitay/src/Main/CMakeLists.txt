#add_executable(main main.f90)
#target_link_libraries(main Input_file General Diagnostics Calculations Problem Mesh BC_Vertex BC_Cell Equation_of_state Quantities Q_Vertex Q_Cell)
#target_link_libraries(main Input_file Calculations Problem Mesh BC_Vertex BC_Cell Equation_of_state Q_Vertex Q_Cell Quantities General)
#target_link_libraries(main General Mesh )
#include_directories($ENV{HDF5_PATH}/include $ENV{MPI_PATH}/include $ENV{JSON_FORTRAN_PATH} $ENV{SILO_PATH}/include $ENV{SCR_PATH}/scr/include)
message($ENV{JSONFORTRAN})
include_directories($ENV{MPI_PATH}/include $ENV{JSONFORTRAN} $ENV{SILO_PATH}/include )
#link_directories($ENV{HDF5_PATH}/lib $ENV{MPI_PATH}/lib $ENV{JSON_FORTRAN_PATH}/lib $ENV{SILO_PATH}/lib $ENV{SCR_PATH}/scr/lib64)
link_directories($ENV{MPI_PATH}/lib $ENV{JSONFORTRAN}/lib $ENV{SILO_PATH}/lib )

# Find CUDA for GPU builds (needed for cudaProfilerStart/Stop)
if(BUILD_TARGET STREQUAL "GPU")
    # For NVHPC compiler, try to find CUDA toolkit
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
        message(STATUS "Found CUDA Toolkit: ${CUDAToolkit_VERSION}")
    else()
        # If CUDAToolkit not found, try to locate it manually for NVHPC
        # NVHPC typically has CUDA in its installation directory
        execute_process(
            COMMAND bash -c "dirname $(which mpif90) | xargs dirname"
            OUTPUT_VARIABLE NVHPC_ROOT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(NVHPC_ROOT)
            set(CUDA_LIB_PATHS
                "${NVHPC_ROOT}/cuda/lib64"
                "${NVHPC_ROOT}/cuda/lib"
                "${NVHPC_ROOT}/compilers/lib"
            )
            foreach(path ${CUDA_LIB_PATHS})
                if(EXISTS "${path}")
                    link_directories(${path})
                    message(STATUS "Added CUDA library path: ${path}")
                endif()
            endforeach()
        endif()
    endif()
endif()

add_library(Main problem.f90)
target_link_libraries(Main mpi z )#silo)
add_executable(main main.f90)
#target_link_libraries(main Diagnostics General)
#target_link_libraries(main Main Diagnostics Mesh Time_step Rezone_and_Advect Input Parallel General BC Quantities BC_Cell BC_Vertex Q_Vertex Q_Cell Material EOS)
target_link_libraries(main Main Diagnostics Mesh Time_step Rezone_and_Advect Input Parallel Quantities Q_Vertex Q_Cell BC BC_Cell BC_Vertex Material EOS General)

# Link CUDA runtime library for GPU builds (needed for cudaProfilerStart/Stop)
if(BUILD_TARGET STREQUAL "GPU")
    # For NVHPC compiler, the CUDA libraries are in the compiler's library path
    # Use -Mcudalib=cudart to link with CUDA runtime via the compiler
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -cuda")
    message(STATUS "Added -cuda linker flag for CUDA runtime library")
endif()

# Set executable name with suffix based on build target (set in parent CMakeLists.txt)
if(DEFINED EXEC_SUFFIX)
    set_target_properties(main PROPERTIES OUTPUT_NAME "main${EXEC_SUFFIX}")
    message(STATUS "Executable will be named: main${EXEC_SUFFIX}")
endif()

#target_link_libraries(main Main Diagnostics Input General Parallel BC BC_Cell BC_Vertex Material Time_step Rezone_and_Advect Mesh Quantities BC_Cell Q_Vertex BC_Vertex Q_Cell Material EOS)

#include_directories($ENV{HDF5_PATH}/include)
#link_directories($ENV{HDF5_PATH}/lib)
#add_library(Main problem.f90)
#target_link_libraries(Main hdf5 z)
