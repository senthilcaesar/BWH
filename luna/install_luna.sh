#!/bin/bash

install_lunaBase=1
install_lunaR=1
fftw_dir=$HOME/Programme/fftw-3.3.8
lgbm_dir=$HOME/Programme/LightGBM
module load R/3.6.3

# Install base luna
if [ ${install_lunaBase} -eq 1 ]
then
    cd $HOME
    git clone https://github.com/remnrem/luna-base.git
    cd luna-base
    make -j4 FFTW=${fftw_dir} LGBM=1 LGBM_PATH=${lgbm_dir}
    cp luna destrat behead fixrows /usr/local/bin
else
    echo "Skipping luna-base installation"
fi

# Install luna R
lib_install_loc=$HOME/R/x86_64-pc-linux-gnu-library/3.6
if [ ${install_lunaR} -eq 1 ]
then
    cd $HOME
    git clone https://github.com/remnrem/luna.git
    FFTW=${fftw_dir} LGBM=1 LGBM_PATH=${lgbm_dir} R CMD INSTALL -l ${lib_install_loc} luna
else
    echo "Skipping luna-R installation"
fi
