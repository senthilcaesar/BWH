#!/bin/bash
install_lunaBase=1
install_lunaR=1
fftw_dir=/opt/homebrew/Cellar/fftw/3.3.10_1
lgbm_dir=/Users/sq566/Downloads/tmp/LightGBM

# Install base luna
if [ ${install_lunaBase} -eq 1 ]
then
    cd $HOME/Programme
    git clone https://github.com/remnrem/luna-base.git
    cd luna-base
    make -j4 ARCH=MAC FFTW=${fftw_dir} LGBM=1 LGBM_PATH=${lgbm_dir}
    cp $HOME/Programme/luna-base/luna /usr/local/bin/luna
    cp $HOME/Programme/luna-base/destrat /usr/local/bin/detstrat
    cp $HOME/Programme/luna-base/behead /usr/local/bin/behead
    cp $HOME/Programme/luna-base/fixrows /usr/local/bin/fixrows

    echo "Luna-base installation successfull"
else
    echo "Skipping luna-base installation"
fi

echo "------------------------------------------------------------------------------"
echo "------------------------------------------------------------------------------"

# Install luna R
if [ ${install_lunaR} -eq 1 ]
then
    cd $HOME/Programme
    git clone https://github.com/remnrem/luna.git
    #FFTW=${fftw_dir} R CMD INSTALL luna
    FFTW=${fftw_dir} LGBM=1 LGBM_PATH=${lgbm_dir} R CMD INSTALL luna
    echo "Luna R installation successfull"
else
    echo "Skipping luna-R installation"
fi
