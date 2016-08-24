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

toName=$1
cd "${Project}"

cut -sf 2 --complement KarenCounts.dat > tempData

head -n +1 tempData | cut -sf 1 --complement > "${toName}"
tail -n +2 tempData >> "${toName}"

rm tempData
