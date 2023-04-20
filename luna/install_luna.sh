#!/bin/bash
# Please run this script as a root user
# Below scripts creates the static libraries ( libfftw3.a and lig_lightgbm.a )
# lib_lightgbm.a and libfftw3.a are baked into executable luna

home='/Users/sq566'

wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
./configure
make
make install

git clone --recursive https://github.com/microsoft/LightGBM
cd LightGBM
mkdir build
cd build
cmake -DUSE_OPENMP=OFF -DBUILD_STATIC_LIB=ON ..
make -j4

git clone https://github.com/remnrem/luna-base.git
cd luna-base
make -j6 ARCH=MAC LGBM=1 LGBM_PATH=$home/Programme/LightGBM FFTW=/usr/local
cp $home/Programme/luna-base/luna /usr/local/bin/luna
cp $home/Programme/luna-base/destrat /usr/local/bin/destrat
cp $home/Programme/luna-base/behead /usr/local/bin/behead
cp $home/Programme/luna-base/fixrows /usr/local/bin/fixrows

