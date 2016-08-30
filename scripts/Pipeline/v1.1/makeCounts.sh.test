#!/bin/bash

Usage="makeCounts.sh [-h] <toName>

where:
    -h          Show this screen
    <toName>	Mandatory argument: Name of Counts file 
                to be saved [default: Counts.dat]"

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
# Making Counts.dat
###############################################################


function makeCounts {
    local toName=$1
    local Postprocessing=$2

	cd "${Postprocessing}"
	
	cut -sf 2 --complement NiceCounts.dat > tempData
	
	head -n +1 tempData | cut -sf 1 --complement > "${toName}"
	tail -n +2 tempData >> "${toName}"
	
	rm tempData
}

postProcessingPath=$1
nameforCounts=$2

makeCounts "${nameforCounts}" "${postProcessingPath}"
