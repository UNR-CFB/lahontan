#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

################################################################
# Preprocessing of reference data
################################################################

function preProcess {
	cd "${Ref}"
	makeblastdb -in "${Cdna}" -dbtype nucl -out ""$basename".cdna.all"
	extract_splice_sites.py "${Gtf}" > splice_sites.txt
	extract_exons.py "${Gtf}" > known_exons.txt
	hisat2-build "${Genome}" "$basename"
	samtools faidx "${Genome}"
}

preProcess
