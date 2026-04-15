#!/bin/bash

# ScalSALE Dual CPU/GPU Build System
# Usage: ./clean.sh <NVHPC|INTEL> [gpu] [cpu]
# Examples:
#   ./clean.sh NVHPC gpu cpu   # Build both versions with NVHPC
#   ./clean.sh NVHPC gpu       # Build GPU version only
#   ./clean.sh NVHPC cpu       # Build CPU version only
#   ./clean.sh NVHPC           # Build both versions by default

# Check if compiler argument is provided
if [ $# -eq 0 ]; then
    echo "ERROR: Compiler not specified"
    echo "Usage: ./clean.sh <NVHPC|INTEL> [gpu] [cpu]"
    exit 1
fi

COMPILER=$1
shift

# Parse build targets from remaining arguments
BUILD_GPU=0
BUILD_CPU=0

# If no arguments after compiler, default to building both versions
if [ $# -eq 0 ]; then
    BUILD_GPU=1
    BUILD_CPU=1
else
    for arg in "$@"; do
        if [ "$arg" = "gpu" ]; then
            BUILD_GPU=1
        elif [ "$arg" = "cpu" ]; then
            BUILD_CPU=1
        else
            echo "Unknown argument: $arg"
            echo "Usage: ./clean.sh [NVHPC|INTEL] [gpu] [cpu]"
            exit 1
        fi
    done
fi

echo "================================================"
echo "ScalSALE Dual Build System"
echo "================================================"
echo "Compiler: $COMPILER"
echo "Build GPU: $BUILD_GPU"
echo "Build CPU: $BUILD_CPU"
echo "================================================"

# Set compiler
if [ "$COMPILER" = "NVHPC" ]; then
    FC_COMPILER=mpif90
    echo "Using NVIDIA HPC Compiler (mpif90)"
else
    FC_COMPILER=mpiifx
    echo "Using Intel Compiler (mpiifx)"
fi

# Function to build a specific target
build_target() {
    local TARGET=$1
    local BUILD_DIR=$2

    echo ""
    echo "========================================"
    echo "Building $TARGET version..."
    echo "========================================"

    # Clean and create build directory
    rm -rf $BUILD_DIR
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR

    # Configure with CMake
    FC=$FC_COMPILER BUILD_TARGET=$TARGET cmake ../src/

    # Build serially (no parallelism) to ensure Fortran module dependencies are respected
    if make; then
        echo ""
        echo "SUCCESS: $TARGET build completed"
        echo "Executable: src/exec/main_${TARGET,,}"
        echo ""
    else
        echo ""
        echo "ERROR: $TARGET build failed"
        echo ""
        cd ../src/Scripts
        exit 1
    fi

    cd ../src/Scripts
}

# Remove old build directories (only for targets being built)
if [ $BUILD_GPU -eq 1 ]; then
    rm -rf ../../build_gpu 2>/dev/null
fi
if [ $BUILD_CPU -eq 1 ]; then
    rm -rf ../../build_cpu 2>/dev/null
fi
rm -rf ../../build 2>/dev/null

# Build GPU version
if [ $BUILD_GPU -eq 1 ]; then
    build_target "GPU" "../../build_gpu"
fi

# Build CPU version
if [ $BUILD_CPU -eq 1 ]; then
    build_target "CPU" "../../build_cpu"
fi

# Clean up
rm core* 2>/dev/null
rm *.so 2>/dev/null

echo ""
echo "================================================"
echo "Build Summary"
echo "================================================"
if [ $BUILD_GPU -eq 1 ]; then
    if [ -f ../exec/main_gpu ]; then
        echo "GPU executable: src/exec/main_gpu [OK]"
    else
        echo "GPU executable: [FAILED]"
    fi
fi

if [ $BUILD_CPU -eq 1 ]; then
    if [ -f ../exec/main_cpu ]; then
        echo "CPU executable: src/exec/main_cpu [OK]"
    else
        echo "CPU executable: [FAILED]"
    fi
fi
echo "================================================"
