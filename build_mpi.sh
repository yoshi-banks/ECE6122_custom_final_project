#!/bin/bash

rm -rf build
mkdir build && cd build

cmake -DUSE_MPI=ON ..

# Build
make -j$(nproc)
