cmake_minimum_required (VERSION 3.5.2)

project(BACKUS)
enable_language (Fortran)

#set(CMAKE_Fortran_COMPILER mpiifx)
#set(FC "mpiifx")
#set(CMAKE_Fortran_FLAGS "-fPIC -O0 -cpp -g -fpp -fpe0 -traceback -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug -check bounds -check format -check output_conversion -check pointers  -check uninit -qopenmp")
#set(CMAKE_Fortran_FLAGS "-fPIC -O2 -cpp -g -fpp -fpe0 -traceback -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug -check bounds -check format -check output_conversion -check pointers  -check uninit -qopenmp")
#set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O2 -qopenmp")
add_definitions($ENV{DEFS})

# Parse build target from environment variable (set by build scripts)
set(BUILD_TARGET $ENV{BUILD_TARGET})
set(SCALSALE_CHECK_MODE $ENV{SCALSALE_CHECK_MODE})
set(SCALSALE_EXTRA_CPU_FLAGS $ENV{SCALSALE_EXTRA_CPU_FLAGS})
set(SCALSALE_EXTRA_GPU_FLAGS $ENV{SCALSALE_EXTRA_GPU_FLAGS})

# Set executable suffix based on build target
if(BUILD_TARGET STREQUAL "GPU")
    set(EXEC_SUFFIX "_gpu")
    add_definitions(-DGPU_BUILD)
    message(STATUS "========================================")
    message(STATUS "Building GPU version (main_gpu)")
    message(STATUS "========================================")
elseif(BUILD_TARGET STREQUAL "CPU")
    set(EXEC_SUFFIX "_cpu")
    add_definitions(-DCPU_BUILD)
    message(STATUS "========================================")
    message(STATUS "Building CPU version (main_cpu)")
    message(STATUS "========================================")
else()
    set(EXEC_SUFFIX "")
    message(WARNING "No BUILD_TARGET specified, building default GPU version")
    add_definitions(-DGPU_BUILD)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
set(ENV{JSONFORTRAN} "/home/reemh/json-fortran/build_nvhpc")

    if(BUILD_TARGET STREQUAL "GPU")
        # GPU-specific flags for NVHPC
        # IF USING NVIDIA A100 GPU:
        if(SCALSALE_CHECK_MODE STREQUAL "1")
            set(CMAKE_Fortran_FLAGS "-fPIC -gpu=cc80,ptxinfo,nofma,noflushz,nounroll -Kieee -Mnofma -Minfo=accel,mp,inline -Mbounds -traceback -mp=gpu -cpp -fopenmp -O0 -g -Minline=levels:10,size:10")
            message(STATUS "Using GPU A100 correctness flags")
        else()
            set(CMAKE_Fortran_FLAGS "-fPIC -gpu=cc80,ptxinfo -Minfo=accel,mp,inline -mp=gpu -cpp -fopenmp -O3 -g -Minline=levels:10,size:10")
            message(STATUS "Using GPU A100 performance flags")
        endif()

        # IF USING NVIDIA V100 GPU:
	#set(CMAKE_Fortran_FLAGS "-fPIC -gpu=cc70,ptxinfo -Minfo=accel,mp,inline -mp=gpu -cpp -fopenmp -O3 -g -Minline=levels:10,size:10")
	#message(STATUS "Using GPU V100")

    elseif(BUILD_TARGET STREQUAL "CPU")
        # CPU-optimized flags for NVHPC (no GPU offloading)
        if(SCALSALE_CHECK_MODE STREQUAL "1")
            set(CMAKE_Fortran_FLAGS "-fPIC -Mnofma -Kieee -mp -Minfo=mp,inline -cpp -fopenmp -O0 -g -Minline=levels:10,size:10")
            message(STATUS "Using NVHPC CPU-only correctness flags")
        else()
            set(CMAKE_Fortran_FLAGS "-fPIC -mp -Minfo=mp,inline -cpp -fopenmp -O3 -g -Minline=levels:10,size:10")
            message(STATUS "Using NVHPC CPU-only performance flags")
        endif()
    else()
        # Default NVHPC flags (GPU)
        set(CMAKE_Fortran_FLAGS "-fPIC -gpu=cc70,ptxinfo -Minfo=accel,mp,inline -mp=gpu -cpp -fopenmp -O3 -g -Minline=levels:10,size:10")
        message(STATUS "Using default NVHPC GPU flags (V100)")
    endif()

#-mavx512f -Ofast
else()
set(ENV{JSONFORTRAN} "/home/u39966a550667980b66e637513213840/json-fortran/json-fortran/build")

    if(BUILD_TARGET STREQUAL "GPU")
        # Intel with GPU offloading (Intel GPUs)
        if(SCALSALE_CHECK_MODE STREQUAL "1")
            set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O0 -g -qopenmp -fopenmp-targets=spir64")
            message(STATUS "Using Intel GPU correctness flags (Ponte Vecchio)")
        else()
            set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O3 -qopenmp -fopenmp-targets=spir64")
            message(STATUS "Using Intel GPU performance flags (Ponte Vecchio)")
        endif()
    elseif(BUILD_TARGET STREQUAL "CPU")
        # Intel CPU-only build
        if(SCALSALE_CHECK_MODE STREQUAL "1")
            set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O0 -g -qopenmp")
            message(STATUS "Using Intel CPU-only correctness flags")
        else()
            set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O3 -qopenmp")
            message(STATUS "Using Intel CPU-only performance flags")
        endif()
    else()
        # Default Intel flags (CPU-only)
        set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O3 -qopenmp")
        message(STATUS "Using default Intel CPU flags")
    endif()
endif()

if(BUILD_TARGET STREQUAL "CPU" AND NOT "${SCALSALE_EXTRA_CPU_FLAGS}" STREQUAL "")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SCALSALE_EXTRA_CPU_FLAGS}")
    message(STATUS "Appending extra CPU flags: ${SCALSALE_EXTRA_CPU_FLAGS}")
elseif(BUILD_TARGET STREQUAL "GPU" AND NOT "${SCALSALE_EXTRA_GPU_FLAGS}" STREQUAL "")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SCALSALE_EXTRA_GPU_FLAGS}")
    message(STATUS "Appending extra GPU flags: ${SCALSALE_EXTRA_GPU_FLAGS}")
endif()

#set(CMAKE_Fortran_FLAGS "-g -qopenmp -fpp -fPIC -traceback -fpe0")
#set(CMAKE_Fortran_FLAGS "-fPIC -O0 -g -fpp -fpe0 -traceback -qopenmp -cpp -fp-model precise -fp-model precise -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug")
add_definitions($ENV{DEFS})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/archive)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/exec)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/Scripts)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(SRC_FILES 
              ${PROJECT_SOURCE_DIR}/Parallel
              ${PROJECT_SOURCE_DIR}/General
              ${PROJECT_SOURCE_DIR}/Input
	     	  ${PROJECT_SOURCE_DIR}/Boundary_conditions
  	     	  ${PROJECT_SOURCE_DIR}/Quantities
	     	  ${PROJECT_SOURCE_DIR}/Boundary_conditions/Vertex
              ${PROJECT_SOURCE_DIR}/Boundary_conditions/Cell
	     	  ${PROJECT_SOURCE_DIR}/Quantities/Cell 
	     	  ${PROJECT_SOURCE_DIR}/Quantities/Vertex
	          ${PROJECT_SOURCE_DIR}/Mesh     	  	     	  
	     	  ${PROJECT_SOURCE_DIR}/Material/Equation_of_state 
	     	  ${PROJECT_SOURCE_DIR}/Material 
	     	  ${PROJECT_SOURCE_DIR}/Rezone_and_Advect
  			  ${PROJECT_SOURCE_DIR}/Time_step
  			  ${PROJECT_SOURCE_DIR}/Diagnostics
#  			  ${PROJECT_SOURCE_DIR}/CR
			  ${PROJECT_SOURCE_DIR}/Main
			  )

foreach(dir ${SRC_FILES})
	add_subdirectory(${dir})
endforeach(dir)
