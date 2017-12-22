#!/bin/bash

Usage="makeSyms.sh [-h]

where:
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
# Making symlinks from original folder to sample folders
###############################################################

function checkDirectory {
    local DIR=$1
    if [ -d "${DIR}" ]; then
        :
    else
        echo "${DIR} does not exist"; exit
    fi
}

function makeOriginalSyms {
    local ogOriginal=$1
    local Original=$2
    checkDirectory "${Original}"

    for file in "${ogOriginal}"/*; do
    	name=$(echo "${file}" | awk -F '/' '{print $NF}')
    
    	ln -s "${file}" "${Original}"/"${name}"
    done
}

function makeDataSyms {
    local Original=$1
    local Data=$2
    checkDirectory "${Data}"

    samplecounter=1
    directorycounter=0
    OriginalFiles=$(echo "${Original}/*" | sort)
    for file in ${OriginalFiles}; do
    	name=$(echo "${file}" | awk -F '/' '{print $NF}')
        samplename=$(printf "sample_%02d" "${samplecounter}")
    
    	ln -sr "${file}" "${Data}"/"${samplename}"/"${name}"
    
    	let directorycounter++
    	if [ $(($directorycounter%2)) -eq 0 ]; then
    		let samplecounter++
    	fi
    done
}

function makeRefSyms {
    local ogReference=$1
    local Reference=$2
    checkDirectory "${Reference}"

    for file in "${ogReference}"/*; do
    	name=$(echo "${file}" | awk -F '/' '{print $NF}')
    
    	ln -s "${file}" "${Reference}"/"${name}"
    done

    if [ -f "${ogReference}/.init" ]; then
        ln -s "${ogReference}/.init" "${Reference}/.init"
    fi
}
function makeJsonSym {
    local JSONfile=$1
    local postprocessingPath=$2
    local name=$3

    checkDirectory "${postprocessingPath}"

    if [ "${name}" == false ]; then
        name=$(echo "${JSONfile}" | awk -F '/' '{print $NF}')
    fi

    ln -s "${JSONfile}" "${postprocessingPath}/${name}"
}
projectPath=$1
originalPath=$2
refPath=$3
jsonPath=$4
checkDirectory "${projectPath}"
checkDirectory "${originalPath}"
checkDirectory "${refPath}"

makeOriginalSyms "${originalPath}" "${projectPath}/Original"
makeDataSyms "${projectPath}/Original" "${projectPath}/Data"
makeRefSyms "${refPath}" "${projectPath}/Reference"

if [ "${jsonPath}" == 'false' ]; then
    :
else
    makeJsonSym "${jsonPath}" "${projectPath}/Postprocessing" "Metadata.json"
fi
