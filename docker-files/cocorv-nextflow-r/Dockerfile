#! Dockerfile for installing tools and copying files and scripts for main CoCoRV package

# Pull base image.
FROM rocker/tidyverse:4.3.1

# Set maintainer.
LABEL  description='Dockerfile for installing tools and copying files and scripts for main CoCoRV package' \
    maintainer="Saima Sultana Tithi, St. Jude Children's Research Hospital, <saimasultana.tithi@stjude.org>"

RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends software-properties-common \
    dirmngr libssl-dev libcurl4-openssl-dev libxml2-dev wget curl \
    bcftools bedtools tabix && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get install -y ca-certificates

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r, download.file.method = 'curl', download.file.extra = '-k -L');" > ~/.Rprofile

RUN R -e 'install.packages(c("DiscreteFDR", "Rcpp", "RcppArmadillo", "data.table", "discreteMTP", "nloptr", "igraph", "devtools", "R.utils", "argparse"), \
    dependencies=TRUE)'

RUN R -e 'install.packages("BiocManager", \
    dependencies=TRUE); \
    library(BiocManager); \
    BiocManager::install(ask = F); \
    BiocManager::install(c("BiocParallel", "IRanges", "GenomicRanges", "SeqArray", "S4Vectors", ask = F))'

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/discreteMTP/discreteMTP_0.1-2.tar.gz", repos = NULL, type="source")'

RUN apt-get update && apt-get install -y libnlopt-dev

RUN R -e 'install.packages("nloptr", \
    dependencies=TRUE)'

# Define working directory.
WORKDIR /opt

# Copy Dockerfile and Changelog.
COPY Dockerfile /opt/
#COPY Changelog /opt/

COPY utilities /opt/utilities
RUN chmod -R a+x /opt/utilities

RUN cd /opt/utilities && \
    Rscript build.R

# by default /bin/bash is executed
CMD ["/bin/bash"]
