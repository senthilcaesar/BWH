# Install FFTW
install_dir="${HOME}/Programme"
mkdir -p "$install_dir"
cd $install_dir

wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
./configure --prefix=$install_dir/fftw-3.3.10/ CFLAGS="-arch arm64"
make
make install

# Install LGBM
cd $install_dir
rm -rf LightGBM
git clone --recursive https://github.com/microsoft/LightGBM
cd LightGBM
mkdir build
cd build
cmake -DUSE_OPENMP=OFF -DBUILD_STATIC_LIB=ON ..
make -j4

# Install luna-base
# Please run this script as a root user
cd $install_dir
rm -rf luna-base
git clone https://github.com/remnrem/luna-base.git
cd luna-base
make -j6 ARCH=MAC LGBM=1 LGBM_PATH=$install_dir/LightGBM FFTW=$install_dir/fftw-3.3.10
cp $HOME/Programme/luna-base/luna /usr/local/bin/luna
cp $HOME/Programme/luna-base/destrat /usr/local/bin/destrat
cp $HOME/Programme/luna-base/behead /usr/local/bin/behead
cp $HOME/Programme/luna-base/fixrows /usr/local/bin/fixrows

# Install luna R
cd $install_dir
rm -rf luna
git clone https://github.com/remnrem/luna.git

FFTW="${install_dir}/fftw-3.3.10/" \
 LGBM=1 \
 LGBM_PATH="${install_dir}/LightGBM/" \
 EXTRA_PKG_LIBS="-L${install_dir}/LightGBM/ -L${install_dir}/fftw-3.3.10/ -l_lightgbm" \
 R CMD INSTALL luna

