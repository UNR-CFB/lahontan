#!/bin/bash

# Assuming your in rna-seq/ directory

RnaSeq=$(pwd)

mkdir "${RnaSeq}/build"

cd "${RnaSeq}/src/GFF/gffread"
make
wait
ln -sr "${RnaSeq}/src/GFF/gffread/gffread" "${RnaSeq}/build/gffread"

cd "${RnaSeq}/src/stringtie"
make release
wait
ln -sr "${RnaSeq}/src/stringtie/stringtie" "${RnaSeq}/build/stringtie"

cd "${RnaSeq}/src/GFF/gffcompare"
make release
wait
ln -sr "${RnaSeq}/src/GFF/gffcompare/gffcompare" "${RnaSeq}/build/gffcompare"

cd "${RnaSeq}/src/samtools"
./configure
wait
make
wait
make prefix="${RnaSeq}/build/samtools" install
wait

cd "${RnaSeq}/src/seqtk"
make
wait
ln -sr "${RnaSeq}/src/seqtk/seqtk" "${RnaSeq}/build/seqtk"

ln -sr "${RnaSeq}/src/stranded_classifier/stranded_classifier.py" "${RnaSeq}/build/stranded_classifier.py"

cd "${RnaSeq}/src/kallisto"
mkdir build
cd build
cmake ..
wait
make
wait
ln -sr "${RnaSeq}/src/kallisto/build/src/kallisto" "${RnaSeq}/build/kallisto"

cd "${RnaSeq}/src/bowtie2"
make
wait
ln -sr "${RnaSeq}/src/bowtie2/bowtie2" "${RnaSeq}/build/bowtie2"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-align-s" "${RnaSeq}/build/bowtie2-align-s"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-align-l" "${RnaSeq}/build/bowtie2-align-l"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-build" "${RnaSeq}/build/bowtie2-build"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-build-s" "${RnaSeq}/build/bowtie2-build-s"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-build-l" "${RnaSeq}/build/bowtie2-build-l"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-inspect" "${RnaSeq}/build/bowtie2-inspect"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-inspect-s" "${RnaSeq}/build/bowtie2-inspect-s"
ln -sr "${RnaSeq}/src/bowtie2/bowtie2-inspect-l" "${RnaSeq}/build/bowtie2-inspect-l"

cd "${RnaSeq}/bin/FastQC"
chmod +x fastqc
ln -sr "${RnaSeq}/bin/FastQC/fastqc" "${RnaSeq}/build/fastqc"

cd "${RnaSeq}/bin"
cp -r hisat2 "${RnaSeq}/build"
cp -r ncbi-blast "${RnaSeq}/build"
cp -r subread "${RnaSeq}/build"
cp -r Trimmomatic "${RnaSeq}/build"

cd "${RnaSeq}"
echo "export RNASEQDIR=${RnaSeq}/build" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/build/hisat2" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/build/ncbi-blast/bin" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/build/samtools/bin" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/build/subread/bin" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/scripts/Pipeline" >> $HOME/.bashrc
echo 'export PATH=$PATH:'"${RnaSeq}/build" >> $HOME/.bashrc
