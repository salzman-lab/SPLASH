# based on existing Docker image
FROM ubuntu:20.04

# install basic software required to install other packages
RUN apt-get update -y --fix-missing
RUN apt-get install -y wget
RUN apt-get install -y unzip

# install python
RUN apt-get install -y python3.9
RUN apt-get install -y python3-pip

# install python packages
COPY requirements.txt /opt/app/requirements.txt
WORKDIR /opt/app
RUN python3 -m pip install -U -r requirements.txt 
WORKDIR /

# make sure python3 is used when 'python' command is called
ENV PATH="/bin/python3.9:${PATH}"

# install bedtools
RUN apt-get install -y bedtools

# install bowtie2
WORKDIR /bin
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip
RUN unzip bowtie2-2.4.5-linux-x86_64.zip
ENV PATH="/bin/bowtie2-2.4.5-linux-x86_64:${PATH}"
WORKDIR /

# install bowtie
WORKDIR /bin
RUN wget https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip
RUN unzip bowtie-1.3.1-linux-x86_64.zip
ENV PATH="/bin/bowtie-1.3.1-linux-x86_64:${PATH}"
WORKDIR /
