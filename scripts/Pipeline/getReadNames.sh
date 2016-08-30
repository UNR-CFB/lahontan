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


###############################################################
# Getting read names from sample folders
###############################################################

function getReadNames {
    local samplename=$1
    local Data=$2

    if [ -d "${Data}"/"${samplename}" ]; then
        > /tmp/readNames."${samplename}"
        for file in "${Data}"/"${samplename}"/*; do
        	echo "${file}" >> /tmp/readNames."${samplename}"
        done
        
        numFiles=$(wc -l < /tmp/readNames."${samplename}")
        if [ "${numFiles}" -eq 2 ]; then
        	awk -F '/' '{print $NF}' /tmp/readNames."${samplename}"
        else
        	echo "Need only specific reads in each sample directory"
            echo "There are "${numFiles}" in "${samplename}"; need only two"
        	cat /tmp/readNames."${samplename}"
        fi
        
        rm /tmp/readNames."${samplename}"
    else
        echo "Invalid sample name: "${samplename}""
    fi
}

sample=$1 #Name of specific sample directory
dataPath=$2

getReadNames "${sample}" "${dataPath}"
