#!/bin/bash
install_lunaBase=0
install_lunaR=0
fftw_dir='/opt/homebrew/Cellar/fftw/3.3.10_1'
lgbm_dir='/Users/sq566/Downloads/tmp/LightGBM'

# Install base luna
if [ ${install_lunaBase} -eq 1 ]
then
    cd $HOME
    git clone https://github.com/remnrem/luna-base.git
    cd luna-base
    make ARCH=MAC FFTW=${fftw_dir} -j4 LGBM=1 LGBM_PATH=${lgbm_dir}
    echo alias luna=$HOME/luna-base/luna >> $HOME/.bash_profile
    echo alias luna=$HOME/luna-base/destrat >> $HOME/.bash_profile
    echo alias luna=$HOME/luna-base/behead >> $HOME/.bash_profile
    source ~/.bash_profile
else
    echo "Skipping luna-base installation"
fi

# Install luna R
if [ ${install_lunaR} -eq 1 ]
then
    cd $HOME
    git clone https://github.com/remnrem/luna.git
    FFTW=${fftw_dir} R CMD INSTALL luna
else
    echo "Skipping luna-R installation"
fi
