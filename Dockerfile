FROM rocker/rstudio:latest


ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir /home/your_giotto

WORKDIR /home/your_giotto

        
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3.6 python3-pip python3-setuptools python3-dev \
    && apt-get install -y r-base-dev \
    && cp /etc/apt/sources.list /etc/apt/sources.list~ \
    && sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list \
    && apt-get update \
    && apt-get -y build-dep r-base-core\
    && pip3 install pandas python-igraph networkx python-louvain leidenalg scikit-learn smfishHmrf

RUN Rscript -e "install.packages('remotes')" \
    && installGithub.r RubD/Giotto \
    && rm -rf /tmp/downloaded_packages/