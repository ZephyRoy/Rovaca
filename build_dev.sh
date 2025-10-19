#!/bin/bash

# Config dynamic library path
export LD_LIBRARY_PATH=./

# Define variables
BUILD_DIR="build-dev"
RELEASE_DIR="release-dev"
JOBS=32

# Exit on error
set -e

# Clean up old build directory
if [ -d "$BUILD_DIR" ]; then
    echo "Cleaning existing build directory..."
    rm -rf "$BUILD_DIR"
fi

if [ -d "$RELEASE_DIR" ]; then
    echo "Cleaning existing release directory..."
    rm -rf "$RELEASE_DIR"
fi

# Make new directories
echo "Creating build directories..."
mkdir -p "$BUILD_DIR" "$RELEASE_DIR"

# Execute build
cd "$BUILD_DIR"
echo "Configuring cmake..."
cmake -DCMAKE_INSTALL_PREFIX="../$RELEASE_DIR" \
    -DLIBDEFLATE_PATH="./" \
    -DBOOST_PATH="./" \
    ..

echo "Building and installing..."
make -j"$JOBS" && make install

echo "Build completed successfully!"