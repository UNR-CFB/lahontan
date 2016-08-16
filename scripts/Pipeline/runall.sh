#!/bin/bash

if [ -f /home/alberton/Pipeline/INPUT ]; then
	source /home/alberton/Pipeline/INPUT
else
	echo "NEED INPUT FILE"
	exit
fi

cd "${Data}"

for directory in $(ls -d */); do
	dir=$(echo "${directory}" | cut -sf 1 -d '/')
	timeit3.sh "${dir}"
done
