#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Making symlinks from original folder to sample folders
###############################################################

samplecounter=1
directorycounter=0
for file in "${Original}"/*; do
	name=$(echo "${file}" | awk -F '/' '{print $NF}')

	ln -s "${file}" "${Data}"/sample_$samplecounter/"${name}"

	let directorycounter++
	if [ $(($directorycounter%2)) -eq 0 ]; then
		let samplecounter++
	else
		:
	fi
done
