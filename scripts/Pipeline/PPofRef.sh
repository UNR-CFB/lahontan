#!/bin/bash

Usage="PPofRef.sh [-h]

where:
    -h          Show this screen"

while getopts ':h' option; do                           
    case "${option}" in                                 
        h) echo "${Usage}"; exit;;                      
        ?)                                              
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;                     
    esac                                                
done                                                    
shift $(( OPTIND -1 ));                                 

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
