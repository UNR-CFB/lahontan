#!/bin/bash

# Possible arguments: ref,data,all,<sample>
#	ref		resets reference directory
#	data		resets every sample directory in $Data
#	all		resets reference and all sample directories
#	<sample> 	resets specific sample directory

if [ -f /home/alberton/Pipeline/INPUT ]; then
	source /home/alberton/Pipeline/INPUT
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Clean Up
###############################################################

function cleanRef {
	mv "${Ref}"/"${Genome}" "${Ref}"/"${Cdna}" "${Ref}"/"${Gtf}" /tmp/
	rm -rf "${Ref}"/
	mkdir "${Ref}"
	mv /tmp/"${Genome}" /tmp/"${Cdna}" /tmp/"${Gtf}" "${Ref}"/
}

function cleanSample {
	local sample=$1
	rm "${Data}"/"${sample}"/*trim*
	mv "${Data}"/"${sample}"/*.gz /tmp/
	rm -rf "${Data}"/"${sample}"/
	mkdir "${Data}"/"${sample}"/
	mv /tmp/*.gz "${Data}"/"${sample}"/	
}

function cleanData {
	for directory in $(ls -d "${Data}"/_*/); do
		cleanSample "${directory}"	
	done
}

function cleanAll {
	cleanRef
	cleanData
}

clean=$1
if [ "${clean}" == "data" ]; then
	cleanData
elif [ "${clean}" == "ref" ]; then
	cleanRef
else
	if [ "${clean}" == "all" ]; then
		cleanAll
	else
		if [ -d "${Data}"/"${clean}" ]; then
			cleanSample "${clean}"				
		else
			echo "Need a real sample directory name"
		fi
	fi
fi
