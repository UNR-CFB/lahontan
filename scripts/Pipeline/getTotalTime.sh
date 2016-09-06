#!/bin/bash

Usage="getTotalTime.sh [-h]

where:
    -h              Shows this screen"

while getopts ':h' option; do                           
    case "${option}" in                                 
        h) echo "${Usage}"; exit;;                      
        ?)                                              
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;                     
    esac                                                
done                                                    
shift $(( OPTIND -1 ));                                 

################################################################
# Pastes together Pipeline Times
################################################################

function getTime {
    local Postprocessing=$1
    local Data=$2

	cd "${Data}"
	
	> "${Postprocessing}"/totalTime.dat
	for directory in */; do
		dir=$(echo "${directory}" | cut -sf 1 -d '/')
		echo "${dir}" >> "${Postprocessing}"/totalTime.dat
		tail -n 4 "${Data}"/"${dir}"/Runtime."${dir}".log >> "${Postprocessing}"/totalTime.dat
	done
	
	echo '' >> "${Postprocessing}"/totalTime.dat
	echo 'Total time is the biggest number out of these:' >> "${Postprocessing}"/totalTime.dat
	grep 'real' "${Postprocessing}"/totalTime.dat | awk '{print $2}' >> "${Postprocessing}"/totalTime.dat
}

function printTime {
    local Postprocessing=$1

    cd "${Postprocessing}"

    maxTime=$(sed -n -e '/biggest/,$p' totalTime.dat | tail -n +2 | sort --reverse | head -n +1)
    echo "${maxTime}"
}

postProcessingPath=$1
dataPath=$2

getTime "${postProcessingPath}" "${dataPath}"
#printTime "${postProcessingPath}"
