#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Making Counts.dat
###############################################################

cd "${Project}"

cut -sf 2 --complement KarenCounts.dat > tempData

head -n +1 tempData | cut -sf 1 --complement > Counts.dat
tail -n +2 tempData >> Counts.dat

rm tempData
