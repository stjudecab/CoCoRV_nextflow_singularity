#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM ensemblorg/ensembl-vep:release_103

# Set maintainer.
LABEL  description='Dockerfile for installing python 3.10, tools and copying files and scripts for main CoCoRV package' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

USER root

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    dirmngr libssl-dev libcurl4-openssl-dev libxml2-dev libopenblas-dev libopenblas-base liblapack3 wget curl git \
    bedtools tabix samtools && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Define working directory.
WORKDIR /opt

# Copy Dockerfile.
COPY Dockerfile /opt/

RUN cd /opt && \
    git config --global http.sslVerify false && \   
    git clone https://bitbucket.org/Wenan/cocorv.git && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    make && make install

ENV BCFTOOLS_PLUGINS="/opt/bcftools/plugins"

RUN apt-get update && apt-get install -y libdbi-perl && \
    curl -k -L https://cpanmin.us | perl - App::cpanminus && \
    cpanm DBD::SQLite && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY annotate_docker.sh /opt/cocorv/utilities
COPY annotateVEPWithOptions_docker.sh /opt/cocorv/utilities
COPY annotateVEPWithOptions_docker_no_mane_v3.sh /opt/cocorv/utilities

RUN chmod -R a+x /opt/cocorv/utilities
RUN chmod -R a+r /opt/cocorv/utilities

COPY annovar /opt/annovar
RUN chmod -R a+x /opt/annovar
RUN chmod -R a+r /opt/annovar

# by default /bin/bash is executed
CMD ["/bin/bash"]
