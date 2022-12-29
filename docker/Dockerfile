FROM rocker/shiny

WORKDIR $HOME/Programme

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get install -y git g++ emacs vim zlib1g-dev libgit2-dev libssl-dev libssh2-1-dev libxml2-dev libcurl4-openssl-dev libfftw3-dev

# Install cmake
ADD https://cmake.org/files/v3.22/cmake-3.22.2-linux-x86_64.sh /cmake-3.22.2-linux-x86_64.sh
RUN mkdir /opt/cmake \
 && sh /cmake-3.22.2-linux-x86_64.sh --prefix=/opt/cmake --skip-license \
 && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake \
 && ln -s /opt/cmake/bin/cmake /usr/bin/cmake \
 && cmake --version

# Install LightGBM
RUN mkdir $HOME/Programme && cd $HOME/Programme \
&& git clone --recursive https://github.com/microsoft/LightGBM \
&& cd LightGBM \
&& mkdir build \
&& cd build \
&& cmake .. \
&& make -j4 

# Install luna-base
RUN cd $HOME/Programme \
 && git clone https://github.com/remnrem/luna-base.git \
 && cd luna-base \
 && make -j 4 FFTW=/usr LGBM=1 LGBM_PATH=$HOME/Programme/LightGBM \
 && ln -s $HOME/Programme/luna-base/luna /usr/local/bin/luna \
 && ln -s $HOME/Programme/luna-base/destrat /usr/local/bin/destrat \
 && ln -s $HOME/Programme/luna-base/behead /usr/local/bin/behead \
 && ln -s $HOME/Programme/luna-base/fixrows /usr/local/bin/fixrows

# Install luna-R dependencies
RUN R -e "install.packages('git2r', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('plotrix', 'geosphere', 'viridis', 'data.table', 'xtable', 'DT', 'shinyFiles', 'shinydashboard', 'lubridate', 'wkb'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('aws.s3')"
RUN R -e "install.packages('shinythemes')"

# Install luna-R
RUN cd $HOME/Programme \
&& git clone https://github.com/remnrem/luna.git \
&& FFTW=/usr R CMD INSTALL luna

RUN echo 'options(defaultPackages=c(getOption("defaultPackages"),"luna" ) )' > ~/.Rprofile

# expose port
EXPOSE 3838

# copy the app directory into the image
RUN mkdir /srv/shiny-server/moonlight
COPY * /srv/shiny-server/moonlight

# run app on container start
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/moonlight', host = '0.0.0.0', port = 3838)"]
