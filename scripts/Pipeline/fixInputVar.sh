#!/bin/bash

Usage="fixInputVar.sh [-h] <path>

where:
    -h          Shows this message
    <path>      Mandatory argument: path to the INPUT file

NOTE: Need to source your .bashrc after running
ex:$ source "${HOME}"/.bashrc"

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
# Change bash environmental variable "Input_Field"
################################################################

inputLocation=$1

function changeInputLocation {
    local inputPATH=$1

    if [ -f "${inputPATH}" ]; then
        if grep -q "Project=" "${inputPATH}"; then
            if grep -q Input_Field "$HOME/.bashrc"; then
            	grep -v "Input_Field" "$HOME/.bashrc" > /tmp/.bashrc
            	mv /tmp/.bashrc "$HOME/.bashrc"
            	echo "export Input_Field="${inputPATH}"" >> $HOME/.bashrc
            else
            	echo "export Input_Field="${inputPATH}"" >> $HOME/.bashrc
            fi
        else
            echo "Need an actual INPUT file"
        fi
    else
        echo "Need path to input file as an argument"
    fi
}

changeInputLocation "${inputLocation}"
