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

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Combining featureCounts output for all samples
###############################################################

function makeNiceCounts {
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

makeNiceCounts
