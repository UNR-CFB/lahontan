#!/bin/bash

PROCS=$(nproc)

# Record hostname
hostname > host.txt

# FastQC run 1
mkdir fastqc1
{ time fastqc -t $PROCS -o fastqc1 *.fq.gz; } > fastqc_1.txt 2>&1

# Trimmomatic
{ time java -jar $RNASEQDIR/Trimmomatic/trimmomatic-0.35.jar \
	PE -threads $PROCS -phred64 \
	read_1.fq.gz read_2.fq.gz \
	read1.P.trim.fq.gz read1.U.trim.fq.gz \
	read2.P.trim.fq.gz read2.U.trim.fq.gz \
	ILLUMINACLIP:$RNASEQDIR/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25; } > trimmomatic.txt 2>&1

# FastQC run 2
mkdir fastqc2
{ time fastqc -t $PROCS -o fastqc2 *.fq.gz; } > fastqc_2.txt 2>&1

# SeqTK sub sample
{ time seqtk sample -s100 read1.P.trim.fq.gz 10000 | seqtk seq -A - > sampled.read1.fa; } > seqtk_1.txt 2>&1
{ time seqtk sample -s100 read2.P.trim.fq.gz 10000 | seqtk seq -A - > sampled.read2.fa; } > seqtk_2.txt 2>&1

# BLAST the sub samples
{ time blastn -query sampled.read1.fa \
	-db ../../../../refs/mus_musculus_20160216/Mus_musculus.GRCm38.cdna.all \
	-out sampled.read1_vscdna.out \
	-task blastn-short -outfmt '6 std sstrand' \
	-max_target_seqs 1 -num_threads $PROCS; } > blastn_1.txt 2>&1
{ time blastn -query sampled.read2.fa \
	-db ../../../../refs/mus_musculus_20160216/Mus_musculus.GRCm38.cdna.all \
	-out sampled.read2_vscdna.out \
	-task blastn-short -outfmt '6 std sstrand' \
	-max_target_seqs 1 -num_threads $PROCS; } > blastn_2.txt 2>&1

# Stranded classifier.  Should only continue if True.  NOT EVALUATED.
{ time stranded_classifier.py -1 sampled.read1_vscdna.out -2 sampled.read2_vscdna.out; } > stranded_classifier.txt 2>&1

# Hisat2
{ time hisat2 -k 5 -p $PROCS --dta --phred64 \
--known-splicesite-infile ../../../../refs/mus_musculus_20160216/splice_sites.txt \
-x ../../../../refs/mus_musculus_20160216/mus_musculus.dna.transformed \
-1 read1.P.trim.fq.gz -2 read2.P.trim.fq.gz \
-S aligned.sam; } > hisat2.txt 2>&1

# SAM -> BAM
{ time samtools view -bT ../../../../refs/mus_musculus_20160216/mus_musculus.dna.transformed.fa \
-@$PROCS aligned.sam -o aligned.bam; } > samtools.txt 2>&1
rm aligned.sam

# To counts
{ time featureCounts -T $PROCS -p -C --primary --ignoreDup -t exon -g gene_id \
-a ../../../../refs/mus_musculus_20160216/Mus_musculus.GRCm38.83.chr.gtf -o aligned.counts aligned.bam; } > featurecounts.txt 2>&1
