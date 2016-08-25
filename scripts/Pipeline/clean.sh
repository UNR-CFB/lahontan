#!/bin/bash

Usage="clean.sh [-h] -r | -d | -a | -p | -s <samplename>

where:
    -h                  Shows this message
    -r		        resets REFERENCE directory
    -d		        resets ALL SAMPLE directories in \$Data
    -a		        resets REFERENCE,POSTPROCESSING and ALL SAMPLE directories
    -p              resets POSTPROCESSING directory
    -s <samplename> 	resets a SPECIFIC SAMPLE directory"

whatToClean=
while getopts ':hrdas:' option; do
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

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"; exit
fi

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
	    rm "${Data}"/"${sample}"/*trim*
	    mv "${Data}"/"${sample}"/*.gz /tmp/
	    rm -rf "${Data}"/"${sample}"/
	    mkdir "${Data}"/"${sample}"/
	    mv /tmp/*.gz "${Data}"/"${sample}"/	
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
    rm -rf "${Postprocessing}"
    mkdir "${Postprocessing}"
}

function cleanAll {
	cleanRef
	cleanData
    cleanPostprocessing
}

function main {
    local thingToClean=$1
    case "${thingToClean}" in
        [reference] ) cleanRef; exit;;
        [data] ) cleanData; exit;;
        [pp] ) cleanPostprocessing; exit;;
        [all] ) cleanAll; exit;;
        * ) cleanSample "${thingToClean}"; exit;;
    esac
}

main "${whatToClean}"
