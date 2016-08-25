#!/bin/bash

Usage="makeStructure.sh [-h] [-f] <num of samples>

where:
    -h                  Shows this screen
    -f                  Force: recursively removes \$Data Directory
                        to make new samples
                        Warning: Will lose all sample data
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

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Structure Maker
###############################################################

numberofSamples=$1

function makeStructure {
    local numSamp=$1

    re='^[0-9]+$' # regular expression to test argument validity
    if [[ "${numSamp}" =~ $re ]] ; then
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
            
            if [ ! -d "${Postprocessing}" ]; then
            	mkdir "${Postprocessing}"
            fi

            for number in `seq -w 1 "${numSamp}"`; do
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
            
            if [ ! -d "${Postprocessing}" ]; then
            	mkdir "${Postprocessing}"
            fi

            for number in `seq -w 1 "${numSamp}"`; do
            	mkdir "${Data}"/sample_$number
            done
        fi
    else
        echo "Need number of samples as an argument" >&2; exit 1
    fi
}

makeStructure "${numberofSamples}"
