#!/bin/bash

Usage="runJSON.sh [-h]

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

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Running JSON scripts
###############################################################

cd "${Postprocessing}"

if [ -f "Metadata.json" ]; then
    while true; do
        read -p "Metadata.json exists, would you like to overwrite Metadata.json?(y,n,createnew) " yn1
        case "${yn1}" in
            [Yy]* ) makeJSON.py;;
            [Nn]* ) break;;
            [c]* )   
                while true; do 
                    read -p "Would you like to create a different Metadata file?(y,n) " yn2
                    case "${yn2}" in
                        [Yy]* ) 
                            read -p "What will the new Metadata file be called? " newName
                            if [ "${newName}" == "Metadata.json" ]; then
                                echo "You need a different name."
                            else
                                makeJSON.py -f "${newName}"; break
                            fi;;
                        [Nn]* ) break;;
                        * ) echo "Please answer yes or no."
                    esac
                done
                break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
else
    makeJSON.py
fi

while true; do
    read -p "Would you like to review the JSON file generated?(y,n) " yn3
    case "${yn3}" in
        [Yy]* ) 
            if [ -f "${newName}" ]; then
                vim "${newName}"
            else
                vim Metadata.json
            fi;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done

if [ -f "${newName}" ]; then
    while true; do
        read -p "Would you like to use the new JSON file to make the stuff?(y,n) " yn4
        case "${yn4}" in
            [Yy]* ) 
                while true; do
                    read -p "Would you like to keep the name of Counts.dat?(y,n) " yn5
                    case "${yn5}" in
                        [Yy]* ) makeCounts.sh "Counts.dat"; break;;
                        [Nn]* )
                            read -p "What should the new Counts file be called?(y,n) " countsName
                            if [ "${countsName}" == 'Counts.dat' ]; then
                                echo 'You need a different name'
                            else
                                makeCounts.sh "${countsName}"; break
                            fi;;
                        * ) echo "Please answer yes or no.";;
                    esac
                done
                while true; do
                    read -p "Would you like to keep the name of Cols.dat?(y,n) " yn6
                    case "${yn6}" in
                        [Yy]* ) makeCols.py -f "${newName}"; break;;
                        [Nn]* )
                            read -p "What should the new Column file be called?(y,n) " colsName
                            if [ "${colsName}" == 'Cols.dat' ]; then
                                echo 'You need a different name'
                            else
                                makeCols.py -f "${newName}" -t "${colsName}"; break
                            fi;;
                        * ) echo "Please answer yes or no.";;
                    esac
                done
                while true; do
                    read -p "Would you like to keep the name of makeReport.r?(y,n) " yn7
                    case "${yn7}" in
                        [Yy]* ) makeReportr.py -f "${newName}"; break;;
                        [Nn]* )
                            read -p "What should the new makeReport.r file be called?(y,n) " mrName
                            if [ "${mrName}" == 'makeReport.r' ]; then
                                echo 'You need a different name'
                            else
                                makeCols.py -f "${newName}" -t "${mrName}"; break
                            fi;;
                        * ) echo "Please answer yes or no.";;
                    esac
                done
                break;;
            [Nn]* )
                makeCounts.sh "Counts.dat"
                makeCols.py
                makeReportr.py
                break;;
            * ) echo "Please answer yes or no.";;
        esac
    done
else
    makeCounts.sh "Counts.dat"
    makeCols.py
    makeReportr.py
fi
