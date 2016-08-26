#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

cd "${Data}"

for directory in */; do
	dir=$(echo "${directory}" | cut -sf 1 -d '/')
	timeit3.sh "${dir}" &
	#sleep 10
done
