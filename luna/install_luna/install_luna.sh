#!/bin/bash
# Please run this script as a root user

rm -rf luna-base
git clone https://github.com/remnrem/luna-base.git
cd luna-base
make -j6 ARCH=MAC LGBM=1 LGBM_PATH=/Users/sq566/Programme/LightGBM FFTW=/Users/sq566/Programme/fftw-3.3.10
cp $HOME/Programme/luna-base/luna /usr/local/bin/luna
cp $HOME/Programme/luna-base/destrat /usr/local/bin/destrat
cp $HOME/Programme/luna-base/behead /usr/local/bin/behead
cp $HOME/Programme/luna-base/fixrows /usr/local/bin/fixrows

