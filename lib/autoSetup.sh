#!/bin/bash

################################################################
# Assuming you're in lahontan/ directory
################################################################
RnaSeq=$(pwd)
Bin="${RnaSeq}/bin"
Source="${RnaSeq}/src"
Lib="${RnaSeq}/lib"

mkdir "${Bin}"

################################################################
# Sources for tools not on Github
################################################################

BLAST_SRC="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz"
TRIM_SRC="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip"
FASTQC_SRC="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip"
FCOUNTS_SRC="https://sourceforge.net/projects/subread/files/subread-1.6.0/subread-1.6.0-Linux-x86_64.tar.gz"
SAMTOOLS_SRC="https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2"

################################################################
# Download Source Files not on Github
################################################################

cd "${Source}"

wget "${BLAST_SRC}"
mkdir "${Source}/ncbi-blast" && tar -xzf "${Source}/ncbi-blast-2.7.1+-x64-linux.tar.gz" -C "${Source}/ncbi-blast" --strip-components 1

wget "${TRIM_SRC}"
unzip "${Source}/Trimmomatic-0.36.zip" && mv "${Source}/Trimmomatic-0.36" "${Source}/trimmomatic"

wget "${FASTQC_SRC}"
unzip "${Source}/fastqc_v0.11.6.zip" && mv "${Source}/FastQC" "${Source}/fastqc"

wget "${FCOUNTS_SRC}"
mkdir "${Source}/subread" && tar -xzf "${Source}/subread-1.6.0-Linux-x86_64.tar.gz" -C "${Source}/subread" --strip-components 1

wget "${SAMTOOLS_SRC}"
mkdir "${Source}/samtools" && tar -xjf "${Source}/samtools-1.6.tar.bz2" -C "${Source}/samtools" --strip-components 1

################################################################
# Installing Tools
################################################################

# Installing gffread
cd "${Source}/GFF/gffread"
make
ln -sr "${Source}/GFF/gffread/gffread" "${Bin}/gffread"

# Installing stringtie
cd "${Source}/stringtie"
make release
ln -sr "${Source}/stringtie/stringtie" "${Bin}/stringtie"

# Installing gffcompare
cd "${Source}/GFF/gffcompare"
make release
ln -sr "${Source}/GFF/gffcompare/gffcompare" "${Bin}/gffcompare"

# Installing samtools
cd "${Source}/samtools"
./configure
make
ln -sr "${Source}/samtools/samtools" "${Bin}/samtools"

# Installing seqtk
cd "${Source}/seqtk"
make
ln -sr "${Source}/seqtk/seqtk" "${Bin}/seqtk"

# Installing stranded_classifier.py
ln -sr "${Source}/stranded_classifier/stranded_classifier.py" "${Bin}/stranded_classifier.py"

# Installing kallisto
cd "${Source}/kallisto"
mkdir build
cd build
cmake ..
make
ln -sr "${Source}/kallisto/build/src/kallisto" "${Bin}/kallisto"

# Installing bowtie2
cd "${Source}/bowtie2"
make
ln -sr "${Source}/bowtie2/bowtie2" "${Bin}/bowtie2"
ln -sr "${Source}/bowtie2/bowtie2-align-s" "${Bin}/bowtie2-align-s"
ln -sr "${Source}/bowtie2/bowtie2-align-l" "${Bin}/bowtie2-align-l"
ln -sr "${Source}/bowtie2/bowtie2-build" "${Bin}/bowtie2-build"
ln -sr "${Source}/bowtie2/bowtie2-build-s" "${Bin}/bowtie2-build-s"
ln -sr "${Source}/bowtie2/bowtie2-build-l" "${Bin}/bowtie2-build-l"
ln -sr "${Source}/bowtie2/bowtie2-inspect" "${Bin}/bowtie2-inspect"
ln -sr "${Source}/bowtie2/bowtie2-inspect-s" "${Bin}/bowtie2-inspect-s"
ln -sr "${Source}/bowtie2/bowtie2-inspect-l" "${Bin}/bowtie2-inspect-l"

# Installing fastqc
cd "${Source}/fastqc"
chmod +x "${Source}/fastqc/fastqc"
ln -sr "${Source}/fastqc/fastqc" "${Bin}/fastqc"

# Installing hisat2
cd "${Source}/hisat2"
make
ln -sr "${Source}/hisat2/hisat2" "${Bin}/hisat2"
ln -sr "${Source}/hisat2/hisat2-align-s" "${Bin}/hisat2-align-s"
ln -sr "${Source}/hisat2/hisat2-align-l" "${Bin}/hisat2-align-l"
ln -sr "${Source}/hisat2/hisat2-build" "${Bin}/hisat2-build"
ln -sr "${Source}/hisat2/hisat2-build-s" "${Bin}/hisat2-build-s"
ln -sr "${Source}/hisat2/hisat2-build-l" "${Bin}/hisat2-build-l"
ln -sr "${Source}/hisat2/hisat2-inspect" "${Bin}/hisat2-inspect"
ln -sr "${Source}/hisat2/hisat2-inspect-s" "${Bin}/hisat2-inspect-s"
ln -sr "${Source}/hisat2/hisat2-inspect-l" "${Bin}/hisat2-inspect-l"
ln -sr "${Source}/hisat2/hisat2-extract_exons.py" "${Bin}/extract_exons.py"
ln -sr "${Source}/hisat2/hisat2-extract_splice_sites.py" "${Bin}/extract_splice_sites.py"

# Installing BLAST
cd "${Source}/ncbi-blast"
ln -sr "${Source}/ncbi-blast/*" "${Bin}/"

# Installing featureCounts
cd "${Source}/subread"
ln -sr "${Source}/subread/bin/featureCounts" "${Bin}/featureCounts"

# Installing Trimmomatic
cd "${Source}/trimmomatic"
ln -sr "${Source}/trimmomatic/trimmomatic-0.36.jar" "${Bin}/trimmomatic-0.36.jar"

################################################################
# Adding to PATH
################################################################

echo "export RNASEQDIR=${Bin}" >> $HOME/.bashrc
echo 'export PATH=$PATH'":${Bin}:${Lib}" >> $HOME/.bashrc
