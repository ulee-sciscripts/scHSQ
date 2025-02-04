FROM ubuntu:20.04 AS builder
RUN apt-get update
ENV TZ America/New_York
ENV DEBIAN_FRONTEND noninteractive
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update
RUN apt-get install -y gnupg ca-certificates libssl-dev curl libcurl4-openssl-dev
RUN echo 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN apt-get update
RUN apt-get install -y build-essential wget git r-base r-base-dev libcairo2-dev libxt-dev
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
RUN Rscript -e 'BiocManager::install()'
RUN Rscript -e 'BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment"))'
RUN apt-get install -y cmake libpq5 libpq-dev
RUN apt-get install -y libgdal-dev libudunits2-dev
RUN apt-get install -y libharfbuzz-dev libfribidi-dev
RUN Rscript -e 'BiocManager::install("HDF5Array", upgrade="always")'
RUN Rscript -e 'remotes::install_github("cole-trapnell-lab/monocle3")'
RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"))'
RUN Rscript -e 'remotes::install_github("ulee-sciscripts/scHSQ/scHSQ", upgrade="always")'

