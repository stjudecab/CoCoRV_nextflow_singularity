#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM ubuntu:22.04

# Set maintainer.
LABEL  description='Dockerfile for installing python 3.10, tools and copying files and scripts for main CoCoRV package' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    dirmngr libssl-dev libcurl4-openssl-dev libxml2-dev libopenblas-dev libopenblas-base liblapack3 wget curl git \
    bcftools bedtools tabix && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y \
    openjdk-8-jre-headless \
    g++ \
    python3.10 python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y build-essential python3.10-dev libpq-dev && \
    python3.10 -m pip install --trusted-host pypi.python.org pip --upgrade && \
    python3.10 -m pip install --trusted-host pypi.python.org psycopg2-binary Cython numpy scikit-learn gnomad==0.7.1 hail==0.2.128 && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y samtools libdbi-perl && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PYSPARK_SUBMIT_ARGS="--driver-memory 40G pyspark-shell"
ENV PERL5LIB="/usr/lib/x86_64-linux-gnu/perl5/5.34/"

# Define working directory.
WORKDIR /opt

# Copy Dockerfile.
COPY Dockerfile /opt/

RUN cd /opt && \
    git config --global http.sslVerify false && \   
    git clone https://bitbucket.org/Wenan/cocorv.git

COPY gnomADPCAndAncestry_docker.sh /opt/cocorv/utilities
COPY annotate_docker.sh /opt/cocorv/utilities
COPY checkFPGenes.sh /opt/cocorv/utilities
COPY postCheckCoCoRV_docker.sh /opt/cocorv/utilities

RUN chmod -R a+x /opt/cocorv/utilities

COPY annovar /opt/annovar
RUN chmod -R a+x /opt/annovar

RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/lib/x86_64-linux-gnu/perl5/5.34 /usr/lib/x86_64-linux-gnu/perl/5.34

# by default /bin/bash is executed
CMD ["/bin/bash"]
