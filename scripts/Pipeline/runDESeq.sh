#!/bin/bash

Usage="runDESeq.sh [-h]

where
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

###############################################################
# Running DESeq
###############################################################

function runR {
    cd "${Postprocessing}"
    { time Rscript "makeReport.r"; } > makeReportTime.log 2>&1 &
}

runR
