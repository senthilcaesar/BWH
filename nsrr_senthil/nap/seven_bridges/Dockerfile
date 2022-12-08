FROM remnrem/luna:latest

WORKDIR /build
ENV DEBIAN_FRONTEND=noninteractive

RUN R -e "install.packages('matlab',repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',repos='http://cran.rstudio.com/')"

RUN cd /build \
 && git clone https://gitlab-scm.partners.org/zzz-public/nsrr.git

CMD [ "/bin/bash" ]

