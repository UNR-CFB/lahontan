#!/bin/bash

Usage="runall.sh [-h] <Data>

where
    -h          Show this screen
    <Data>      Mandatory argument: path to Data directory"

while getopts ':h' option; do
    case "${option}" in
        h) echo "${Usage}"; exit;;
        ?)
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;
    esac
done
shift $(( OPTIND - 1 ));

function checkDirectory {
    local DIR=$1
    if [ -d "${DIR}" ]; then
        :
    else
        echo "${DIR} does not exist"; exit
    fi
}

Data=$1
fastq=$2
PROCS=$3
Ref=$4
Genome=$5
basename=$6
Project=$7

checkDirectory "${Data}"
checkDirectory "${Data}"/"${sample}"
checkDirectory "${Ref}"

mkdir "${Project}"/runPipeNotify

cd "${Data}"

for directory in */; do
	dir=$(echo "${directory}" | cut -sf 1 -d '/')
	timeit3.sh "${dir}" "${Data}" "${fastq}" "${PROCS}" "${Ref}" "${Genome}" "${basename}" "${Project}" &
done
