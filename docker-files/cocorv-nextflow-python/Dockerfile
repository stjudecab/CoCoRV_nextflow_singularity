#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM ubuntu:22.04

# Set maintainer.
LABEL  description='Dockerfile for installing python 3.10, tools and copying files and scripts for main CoCoRV package' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    dirmngr libssl-dev libcurl4-openssl-dev libxml2-dev libopenblas-dev libopenblas-base liblapack3 wget curl \
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

ENV PYSPARK_SUBMIT_ARGS="--driver-memory 40G pyspark-shell"

# Define working directory.
WORKDIR /opt

# Copy Dockerfile and Changelog.
COPY Dockerfile /opt/
#COPY Changelog /opt/

# by default /bin/bash is executed
CMD ["/bin/bash"]
