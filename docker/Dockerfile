# Start from the original phip-flow Dockerfile base image
FROM quay.io/hdc-workflows/ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing -qq && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y -q \
    git \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    python3 \
    python3-pip \ 
    python3-venv \
    zip \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    gnupg \
    lsb-release && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# install phippery
RUN pip install git+https://github.com/matsengrp/phippery@1.3.1

RUN pip3 install pandas

# install pre-build binary Bowtie1.3
RUN curl -fksSL https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip \
    --output bowtie-1.3.1-linux-x86_64.zip \
    && unzip bowtie-1.3.1-linux-x86_64.zip \
    && (cd /usr/bin/ && ln -s /bowtie-1.3.1-linux-x86_64/* ./)


# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

# Install Java 17 (required for Nextflow) and system build tools
RUN apt-get update && apt-get install -y \
    openjdk-17-jdk \
    curl \
    unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set environment variables for Java
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# Install R and some common system dependencies for R packages
RUN apt-get update && apt-get install -y software-properties-common dirmngr gnupg apt-transport-https ca-certificates && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 'E298A3A825C0D65DFD57CBB651716619E084DAB9' && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y r-base

# Install JAGS from source
WORKDIR /tmp

ADD https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.2.tar.gz/download JAGS-4.3.2.tar.gz

RUN tar -xzf JAGS-4.3.2.tar.gz && \
    cd JAGS-4.3.2 && \
    ./configure && \
    make && \
    make install && \
    ldconfig && \
    cd / && rm -rf /tmp/*

# Install R packages: rjags, edgeR, BEER (plus BiocManager)
RUN Rscript -e 'install.packages(c("BiocManager", "rjags"), repos="https://cloud.r-project.org")' && \
    Rscript -e 'BiocManager::install("edgeR")' && \
    Rscript -e 'BiocManager::install("beer")'

# Set working directory
WORKDIR /phipflow

RUN git clone https://github.com/matsengrp/phip-flow.git /phipflow && \
    mkdir /phipflow_selected && \
    mv /phipflow/bin /phipflow_selected/ && \
    mv /phipflow/workflows /phipflow_selected/ && \
    mv /phipflow/templates /phipflow_selected/ && \
    mv /phipflow/nextflow.config /phipflow_selected/ && \
    rm -rf /phipflow && \
    mv /phipflow_selected /phipflow 

RUN mkdir /phipflow/data
RUN mkdir /phipflow/data/seq

COPY main.nf /phipflow/main.nf

# Default command
CMD ["/bin/bash"]
