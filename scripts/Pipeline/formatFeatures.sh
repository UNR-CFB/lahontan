#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Combining featureCounts output for all samples
###############################################################

mkdir "${Data}"/../ProjectCounts

cd "${Data}"

counter=1
for directory in $(ls -d sample*/); do
	dir=$(echo "${directory}" | cut -sf 1 -d '/')
	if [ "${counter}" -eq 1 ]; then
		cp "${Data}"/"${dir}"/aligned."${dir}".counts.three "${Data}"/../ProjectCounts
		let counter++
	else
		cp "${Data}"/"${dir}"/aligned."${dir}".counts.one "${Data}"/../ProjectCounts
		let counter++
	fi	
done

paste "${Data}"/../ProjectCounts/*.three "${Data}"/../ProjectCounts/*.one > "${Data}"/../Counts.dat

rm -rf "${Data}"/../ProjectCounts/
