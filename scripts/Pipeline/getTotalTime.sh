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

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

################################################################
# Pastes together Pipeline Times
################################################################

def getTime {
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

getTime
