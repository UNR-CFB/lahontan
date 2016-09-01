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

###############################################################
# Running DESeq
###############################################################

function runR {
    local Postprocessing=$1

    cd "${Postprocessing}"
    { time Rscript "makeReport.r"; } > makeReportTime.log 2>&1 &
}
function askx11 {
    while true; do
        read -p "Is the x server available or is X11 forwarding enabled?(y,n) " answer
        case "${answer}" in
            [Yy]* )
                break;;
            [Nn]* )
                echo "You must enable X11 forwarding or run from local host machine"
                exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done

}

postProcessingPath=$1
#askx11
runR "${postProcessingPath}"
