Bootstrap: debootstrap
OSVersion: xenial
MirrorURL: http://us.archive.ubuntu.com/ubuntu/

###############################################################
# Command used to build:
# sudo singularity build runpipe.simg runpipeRecipe
###############################################################

%labels
    MAINTAINER Alberto
    VERSION v1.1

%post
    locale-gen en_US.UTF-8

    echo "Installing necessary packages..."
    echo "deb http://us.archive.ubuntu.com/ubuntu/ xenial universe" >> /etc/apt/sources.list
    apt-get update && apt-get install --yes git cmake gcc g++ libncurses-dev libhdf5-cpp-11 libhdf5-dev python3-docopt vim unzip openjdk-8-jdk-headless wget gfortran libbz2-dev liblzma-dev libpcre++-dev libcurl4-openssl-dev libssl-dev pandoc texlive-latex-extra libxml2-dev libmariadb-client-lgpl-dev libreadline6-dev libreadline6 libtbb-dev

    cd /opt
    wget https://cran.r-project.org/src/base/R-3/R-3.4.3.tar.gz
    tar xzf R-3.4.3.tar.gz
    cd R-3.4.3
    ./configure --with-x=no 
    make && make install
    R -e 'source("https://bioconductor.org/biocLite.R");biocLite(ask=FALSE);biocLite(c("devtools","DESeq2","edgeR","ReportingTools","regionReport","pachterlab/sleuth","ballgown","DT","pheatmap"));devtools::install_github(c("docopt/docopt.R","alyssafrazee/RSkittleBrewer"))'
    
    echo "Cloning repository..."
    git clone --recursive https://Alberto024:testpassword1234@github.com/UNR-CFB/rna-seq.git /rna-seq

    echo "Installing pipeline..."
    cd /rna-seq
    echo 'export RNASEQDIR=/rna-seq/build' >> $SINGULARITY_ENVIRONMENT
    /rna-seq/scripts/Pipeline/autoSetup.sh

%environment
    RNASEQDIR=/rna-seq/build
    PATH="${PATH}:/rna-seq/build:/rna-seq/build/hisat2:/rna-seq/build/ncbi-blast/bin:/rna-seq/build/samtools/bin:/rna-seq/build/subread/bin:/rna-seq/scripts/Pipeline"
    LANG=en_US.UTF-8
    LANGUAGE=en_US
    export RNASEQDIR PATH LANG LANGUAGE

%help
    Go to https://github.com/UNR-CFB/rna-seq for more help OR try "./runpipe.simg -h" for pipeline help

%runscript
    exec /rna-seq/scripts/Pipeline/runPipe "$@"
