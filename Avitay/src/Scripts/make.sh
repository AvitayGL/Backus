#!/bin/bash

# ScalSALE Incremental Build System
# Usage: ./make.sh [gpu] [cpu]
# Examples:
#   ./make.sh gpu cpu   # Rebuild both versions
#   ./make.sh gpu       # Rebuild GPU version only
#   ./make.sh cpu       # Rebuild CPU version only
#   ./make.sh           # Rebuild GPU version by default

# Parse build targets
BUILD_GPU=0
BUILD_CPU=0

# If no arguments, default to rebuilding both versions
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
            echo "Usage: ./make.sh [gpu] [cpu]"
            exit 1
        fi
    done
fi

echo "================================================"
echo "ScalSALE Incremental Build"
echo "================================================"
echo "Build GPU: $BUILD_GPU"
echo "Build CPU: $BUILD_CPU"
echo "================================================"

# Function to rebuild a specific target
rebuild_target() {
    local TARGET=$1
    local BUILD_DIR=$2

    if [ ! -d "$BUILD_DIR" ]; then
        echo ""
        echo "ERROR: $BUILD_DIR does not exist"
        echo "Run ./clean.sh first to create initial build"
        return 1
    fi

    echo ""
    echo "========================================"
    echo "Rebuilding $TARGET version..."
    echo "========================================"

    cd $BUILD_DIR

    # Build serially (no parallelism) to ensure Fortran module dependencies are respected
    if make; then
        echo ""
        echo "SUCCESS: $TARGET rebuild completed"
        echo ""
    else
        echo ""
        echo "ERROR: $TARGET rebuild failed"
        echo ""
        cd ../src/Scripts
        exit 1
    fi

    cd ../src/Scripts
}

# Rebuild GPU version
if [ $BUILD_GPU -eq 1 ]; then
    rebuild_target "GPU" "../../build_gpu"
fi

# Rebuild CPU version
if [ $BUILD_CPU -eq 1 ]; then
    rebuild_target "CPU" "../../build_cpu"
fi

echo ""
echo "================================================"
echo "Rebuild Summary"
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
