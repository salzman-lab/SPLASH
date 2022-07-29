# based on existing Docker image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget && apt-get install -y --no-install-recommends build-essential r-base python3.9 python3-pip python3-setuptools python3-dev
WORKDIR /app
COPY requirements.txt /app/requirements.txt
RUN export TMPDIR='/var/tmp'
RUN pip3 install --no-cache-dir -r requirements.txt

# install dependencies and clean up apt garbage
RUN apt-get update && apt-get install --no-install-recommends -y \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libgnutls-openssl27 \
    zlib1g-dev \
    libssl-dev \
    gcc \
    wget \
    make \
    perl \
    bzip2 \
    gnuplot \
    ca-certificates \
    libopenblas-dev \
    gawk && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    autoconf \
    automake \
    g++ \
    gcc \
    gfortran \
    make \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('spgs')"
RUN Rscript -e "install.packages('stringdist')"
RUN Rscript -e "install.packages('R.utils')"
RUN Rscript -e "install.packages('reshape2')"

COPY . /app

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda
ENV PATH=$PATH:/miniconda/condabin:/miniconda/bin
RUN conda install -c bioconda blast
RUN conda install -c bioconda/label/cf201901 bedtools

# install bowtie2
WORKDIR /bin
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip
RUN unzip bowtie2-2.4.5-linux-x86_64.zip
ENV PATH="/bin/bowtie2-2.4.5-linux-x86_64:${PATH}"
WORKDIR /

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar -xjf samtools-1.15.1.tar.bz2 && \
    rm samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure && \
    make && \
    make install


