wget https://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvzf fftw-3.3.10.tar.gz
cd fftw-3.3.10
install_dir="${HOME}/Programme"

if [ ! -d "$install_dir" ]; then
    mkdir -p "$install_dir"
fi

./configure --prefix=$install_dir/fftw-3.3.10/ CFLAGS="-arch arm64"
make
make install
cp $install_dir/fftw-3.3.10/lib/libfftw3.a /usr/local/lib
