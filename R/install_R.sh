#!/bin/bash

R_VERSION=3.6.3
curl -O https://cran.rstudio.com/src/base/R-3/R-${R_VERSION}.tar.gz
tar -xzvf R-${R_VERSION}.tar.gz
cd R-${R_VERSION}

export PATH=$HOME/Programme/bzip2-1.0.8/bin:$PATH
export PATH=$HOME/Programme/xz-5.2.8/bin:$PATH
export PATH=$HOME/Programme/zlib-1.2.13/bin:$PATH
export PATH=$HOME/Programme/pcre-8.40/bin:$PATH
export PATH=$HOME/Programme/curl-7.86.0/bin:$PATH
export PATH=$HOME/Programme/openssl/bin:$PATH  

export JAVA_HOME=$HOME/java/jdk/jdk-11.0.12+7
export PATH="$JAVA_HOME/bin:$PATH"

export LD_LIBRARY_PATH=$HOME/Programme/bzip2-1.0.8/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/Programme/xz-5.2.8/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/Programme/zlib-1.2.13/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/Programme/pcre-8.40/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/Programme/curl-7.86.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/Programme/openssl/lib:$LD_LIBRARY_PATH    

export CFLAGS="-I$HOME/Programme/bzip2-1.0.8/include -I$HOME/Programme/xz-5.2.8/include -I$HOME/Programme/zlib-1.2.13/include -I$HOME/Programme/pcre-8.40/include -I$HOME/Programme/curl-7.86.0/include -I$HOME/Programme/openssl/include"
export LDFLAGS="-L$HOME/Programme/bzip2-1.0.8/lib -L$HOME/Programme/xz-5.2.8/lib -L$HOME/Programme/zlib-1.2.13/lib -L$HOME/Programme/pcre-8.40/lib -L$HOME/Programme/curl-7.86.0/lib -L$HOME/Programme/openssl/lib"

./configure \
--prefix=/PHShome/sq566/R/R-${R_VERSION} \
--enable-memory-profiling \
--enable-R-shlib \
--with-blas \
--with-lapack \
--with-readline=no --with-x=no

make
make install

FFTW=/PHShome/sq566/Programme/fftw-3.3.8 R CMD INSTALL -l /PHShome/sq566/R/x86_64-pc-linux-gnu-library/3.6 luna
