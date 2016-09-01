#!/bin/bash

Usage="makeStructure.sh [-h] [-f] <Projectpath> <num of samples>

where:
    -h                  Shows this screen
    -f                  Force: recursively removes \$Data Directory
                        to make new samples
                        Warning: Will lose all sample data
    <Projectpath>       Path to Project Directory
    <num of samples>	An integer: Amount of samples in experiment"

force=false
while getopts ':hf' option; do
    case "${option}" in
        h) echo "${Usage}"; exit;;
        f) force=true; break;;
        ?) 
            printf "Invalid Argument: -%s\n" "${OPTARG}"
            echo "${Usage}"; exit;;
    esac
done
shift $(( OPTIND -1 ));

###############################################################
# Structure Maker
###############################################################

function makeStructure {
    local Project=$1
    local numSamp=$2
    local Data="${Project}/Data"
    local Ref="${Project}/Reference"
    local Original="${Project}/Original"
    local Postprocessing="${Project}/Postprocessing"

    if [ "${force}" = true ]; then
        if [ ! -d "${Project}" ]; then
        	mkdir "${Project}"
        fi

        if [ ! -d "${Data}" ]; then
        	mkdir "${Data}"
        else
            rm -rf "${Data}"
            mkdir "${Data}"
        fi
        
        if [ ! -d "${Ref}" ]; then
        	mkdir "${Ref}"
        fi
        
        if [ ! -d "${Original}" ]; then
        	mkdir "${Original}"
        fi

        if [ ! -d "${Postprocessing}" ]; then
        	mkdir "${Postprocessing}"
        fi

        for number in `seq -f "%02g" 1 "${numSamp}"`; do
        	mkdir "${Data}"/sample_$number
        done
    else
        if [ ! -d "${Project}" ]; then
        	mkdir "${Project}"
        fi
        
        if [ ! -d "${Data}" ]; then
        	mkdir "${Data}"
        fi
        
        if [ ! -d "${Ref}" ]; then
        	mkdir "${Ref}"
        fi

        if [ ! -d "${Original}" ]; then
        	mkdir "${Original}"
        fi
        
        if [ ! -d "${Postprocessing}" ]; then
        	mkdir "${Postprocessing}"
        fi

        for number in `seq -f "%02g" 1 "${numSamp}"`; do
        	mkdir "${Data}"/sample_$number
        done
    fi
}

projectPath=$1
numberofSamples=$2
makeStructure "${projectPath}" "${numberofSamples}"
