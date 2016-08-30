#!/bin/bash

Usage="formatFeatures.sh [-h]

where:
    -h          Shows this screen"

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
# Combining featureCounts output for all samples
###############################################################

function makeNiceCounts {
    local Postprocessing=$1
    local Data=$2

    mkdir "${Postprocessing}"/ProjectCounts
    
    cd "${Data}"
    
    counter=1
    for directory in $(ls -d sample*/); do
    	dir=$(echo "${directory}" | cut -sf 1 -d '/')
    	if [ "${counter}" -eq 1 ]; then
    		cp "${Data}"/"${dir}"/aligned."${dir}".counts.three "${Postprocessing}"/ProjectCounts
    		let counter++
    	else
    		cp "${Data}"/"${dir}"/aligned."${dir}".counts.one "${Postprocessing}"/ProjectCounts
    		let counter++
    	fi	
    done
    
    paste "${Postprocessing}"/ProjectCounts/*.three "${Postprocessing}"/ProjectCounts/*.one > "${Postprocessing}"/NiceCounts.dat
    
    rm -rf "${Postprocessing}"/ProjectCounts/
}

postProcessingPath=$1
dataPath=$2

makeNiceCounts "${postProcessingPath}" "${dataPath}"
