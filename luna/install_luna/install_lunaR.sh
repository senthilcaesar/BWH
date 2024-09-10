rm -rf luna
git clone https://github.com/remnrem/luna.git

FFTW="/Users/sq566/Programme/fftw-3.3.10/" \
 LGBM=1 \
 LGBM_PATH="/Users/sq566/Programme/LightGBM/" \
 EXTRA_PKG_LIBS="-L/Users/sq566/Programme/LightGBM/ -L/Users/sq566/Programme/fftw-3.3.10/ -l_lightgbm" \
 R CMD INSTALL luna

