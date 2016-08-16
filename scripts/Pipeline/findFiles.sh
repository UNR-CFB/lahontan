#!/bin/bash

if [ -f /home/alberton/Pipeline/INPUT ]; then
	source /home/alberton/Pipeline/INPUT
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Getting read names from sample folders
###############################################################

sample=$1 #Name of specific sample directory
#Data="/home/alberton/intern-demo-At/fastq/find_test"

> /tmp/readNames
for file in "${Data}"/"${sample}"/*; do
	echo "${file}" >> /tmp/readNames
done

numFiles=$(wc -l < /tmp/readNames)
if [ "${numFiles}" -eq 2 ]; then
	awk -F '/' '{print $NF}' /tmp/readNames
else
	echo "Need only specific reads in each sample directory"
	cat /tmp/readNames
fi

rm /tmp/readNames
