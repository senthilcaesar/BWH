#!/bin/bash

install_dir="${HOME}/Programme"

if [ ! -d "$install_dir" ]; then
    mkdir -p "$install_dir"
fi

cd $install_dir
rm -rf LightGBM
git clone --recursive https://github.com/microsoft/LightGBM
cd LightGBM
mkdir build
cd build
cmake -DUSE_OPENMP=OFF -DBUILD_STATIC_LIB=ON ..
make -j4
