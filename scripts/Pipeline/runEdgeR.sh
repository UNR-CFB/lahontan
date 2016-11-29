#!/bin/bash

Usage="runEdgeR.sh [-h]

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

###############################################################
# Running DESeq
###############################################################

function runR {
    local Postprocessing=$1

    cd "${Postprocessing}"
    { time Rscript "makeEdge.r"; } > makeEdgeTime.log 2>&1
    echo '' > "${Postprocessing}/.doneE"
}

postProcessingPath=$1
runR "${postProcessingPath}"
