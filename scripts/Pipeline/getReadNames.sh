#!/bin/bash

Usage="getReadNames.sh [-h] <samplename>

where:
    -h                  Show this message
    <samplename>        Mandatory argument of Sample Name"

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

###############################################################
# Getting read names from sample folders
###############################################################

sample=$1 #Name of specific sample directory

function getReadNames {
    local samplename=$1

    if [ -d "${Data}"/"${samplename}" ]; then
        > /tmp/readNames
        for file in "${Data}"/"${sample}"/*; do
        	echo "${file}" >> /tmp/readNames
        done
        
        numFiles=$(wc -l < /tmp/readNames)
        if [ "${numFiles}" -eq 2 ]; then
        	awk -F '/' '{print $NF}' /tmp/readNames
        else
        	echo "Need only specific reads in each sample directory"
            echo "There are "${numFiles}" in "${sample}"; need only two"
        	cat /tmp/readNames
        fi
        
        rm /tmp/readNames
    else
        echo "Invalid sample name: "${samplename}""
    fi
}

getReadNames "${sample}"
