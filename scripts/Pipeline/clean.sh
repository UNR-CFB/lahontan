#!/bin/bash

Usage="clean.sh [-h] -r | -d | -a | -p | -s <samplename>

where:
    -h                  Shows this message
    -r		        resets REFERENCE directory
    -d		        resets ALL SAMPLE directories in \$Data
    -a		        resets REFERENCE,POSTPROCESSING and ALL SAMPLE directories
    -p              resets POSTPROCESSING directory
    -s <samplename> 	resets a SPECIFIC SAMPLE directory"

whatToClean=false
while getopts ':hrdaps:' option; do
    case "$option" in
        h) echo "${Usage}"; exit;;
        r) whatToClean='reference'; break;;
        d) whatToClean='data'; break;;
        a) whatToClean='all'; break;;
        p) whatToClean='pp'; break;;
        s) whatToClean="${OPTARG}"; break;;
        ?) 
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;
    esac
done
shift $(( OPTIND - 1 ));

###############################################################
# Clean Up
###############################################################

function cleanRef {
	mv "${Ref}"/"${Genome}" "${Ref}"/"${Cdna}" "${Ref}"/"${Gtf}" /tmp/
	rm -rf "${Ref}"/
	mkdir "${Ref}"
	mv /tmp/"${Genome}" /tmp/"${Cdna}" /tmp/"${Gtf}" "${Ref}"/
}

function cleanSample {
	local sample=$1
    if [ -d "${Data}"/"${sample}" ]; then
        if [ -f "${Data}"/"${sample}"/*trim* ]; then
	        rm "${Data}"/"${sample}"/*trim*
        fi
        cd "${Data}/${sample}"
        rm $(ls -I "*.gz")
    else
        echo "Invalid sample name: "${sample}""
    fi
}

function cleanData {
	cd "${Data}"
	for directory in */; do
		dir=$(echo "${directory}" | cut -sf 1 -d '/')
		cleanSample "${dir}"
	done
}

function cleanPostprocessing {
    if [ -f "${Postprocessing}/Metadata.json" ]; then
        mv "${Postprocessing}/Metadata.json" /tmp/Metadata.json
        rm -rf "${Postprocessing}"
        mkdir "${Postprocessing}"
        mv /tmp/Metadata.json "${Postprocessing}/Metadata.json"
    else
        rm -rf "${Postprocessing}"
        mkdir "${Postprocessing}"
    fi
}

function cleanAll {
	cleanRef
	cleanData
    cleanPostprocessing
}

function main {
    local thingToClean=$1
    case "${thingToClean}" in
        [r]* ) cleanRef; exit;;
        [d]* ) cleanData; exit;;
        [p]* ) cleanPostprocessing; exit;;
        [a]* ) cleanAll; exit;;
        [s]* ) cleanSample "${thingToClean}"; exit;;
    esac
}

Genome=$1
Cdna=$2
Gtf=$3
Ref=$4
Data=$5
Postprocessing=$6

if [ "${whatToClean}" == false ]; then
    echo "Need valid argument"
    echo "${Usage}"; exit
else
    main "${whatToClean}"
fi
