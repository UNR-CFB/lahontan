#!/bin/bash
### PROJECT-WIDE VARIABLES SECTION
### these variables were used for a particular organization scheme for the
### beetle projects. Structure is probably a bit too over-wrought
# work_root=/home/rltillett/data/3beetles/basespace
# sequence_dir=$work_root/glued
# hisat_dir=$work_root/dendro_hisat
# hisat_ref=$hisat_dir/dendro_hisat
# known_genes=$hisat_dir/dendro_known.txt
# out_dir=$work_root/hisat_mas
# scripts_dir=$work_root/scripts
# prefix_file=$scripts_dir/dendro_prefix

### EXTRA FANCY LOOP SECTION
### Auto-generation of the sample input prefixes for reads 1 and 2
### empty version of the loop syntax for loaded-by-file is also appended to 
### this comment
# prefix_file=$scripts_dir/dendro_prefix
# cd $sequence_dir
# ls -1 *.gz | cut -d '.' -f1 | sort | uniq > $prefix_file
# cd $work_root
# while read -r prefix <&9; do
# <SOME CMDS ITERATING ON "${prefix}">
# done 9< ${prefix_file}

### EXPLICIT LOOP SECTION
### Simplest possible bash looping method below
### if input files are named like "sample1_r1.fq.gz" then a logical 
### place to break out a variable named "$prefix" is "sample1" 
for prefix in \
sample1 \
sample2 \
sample11 \
sampleq1qxxx
do
hisat -p 30 --known-splicesite-infile $known_genes -x $hisat_ref -1 "${sequence_dir}/${prefix}_r1.fq.gz" -2 "${sequence_dir}/${prefix}_r2.fq.gz" -S "${out_dir}/${prefix}.sam" &> "${out_dir}/${prefix}.log"
samtools view -bS "${out_dir}/${prefix}.sam" > "${out_dir}/${prefix}.notSRs.bam"
pigz "${out_dir}/${prefix}.notSRs.sam"
done

