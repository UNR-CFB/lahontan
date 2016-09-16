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
    cd "${Data}/${sample}"
    if [ -d "${Data}"/"${sample}" ]; then
        numero=$(echo * | wc -w)
        if [ "${numero}" -eq 2 ]; then
            :
        else
            pattern="${Data}/${sample}/*trim*"
            files=( $pattern )
            if [ -f "${files[0]}" ]; then
	            rm "${Data}/${sample}"/*trim*
            fi
            pattern2="${Data}/${sample}/fastq*/"
            files2=( ${pattern2} )
            if [ -d "${files2[0]}" ]; then
                #counter=$((counter+1))
                counter=0
                while true; do
                    if [ -d "${files2[$counter]}" ]; then
                        rm -rf "${files2[$counter]}"
                        let counter++
                    else
                        break
                    fi
                done
            fi
            cd "${Data}/${sample}"
            if [ "${numero}" -eq 2 ]; then
                :
            else
                rm $(ls -I "*.gz")
            fi
        fi
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
