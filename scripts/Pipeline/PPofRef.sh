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

################################################################
# Preprocessing of reference data
################################################################

function checkDirectory {
    local DIR=$1
    if [ -d "${DIR}" ]; then
        :
    else
        echo "${DIR} does not exist"; exit
    fi
}
function checkFile {
    local Location=$1
    local Filename=$2
    if [ -f "${Location}/${Filename}" ]; then
        :
    else
        echo "${Location}/${Filename} does not exist"; exit
    fi
}

function preProcess {
    local Ref=$1
    local Cdna=$2
    local Gtf=$3
    local Genome=$4
    local basename=$5

	cd "${Ref}"

	makeblastdb -in "${Cdna}" -dbtype nucl -out ""$basename".cdna.all"
	extract_splice_sites.py "${Gtf}" > splice_sites.txt
	extract_exons.py "${Gtf}" > known_exons.txt
	hisat2-build "${Genome}" "$basename"
	samtools faidx "${Genome}"
}

refPath=$1
cdnaName=$2
gtfName=$3
genomeName=$4
baseName=$5
checkDirectory "${refPath}"
checkFile "${refPath}" "${cdnaName}"
checkFile "${refPath}" "${gtfName}"
checkFile "${refPath}" "${genomeName}"

preProcess "${refPath}" "${cdnaName}" "${gtfName}" "${genomeName}" "${baseName}"
