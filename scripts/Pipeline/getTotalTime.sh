#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

cd "${Data}"

> "${Project}"/totalTime.dat
for directory in */; do
	dir=$(echo "${directory}" | cut -sf 1 -d '/')
	echo "${dir}" >> "${Project}"/totalTime.dat
	tail -n 4 "${Data}"/"${dir}"/Runtime."${dir}".log >> "${Project}"/totalTime.dat
done

echo '' >> "${Project}"/totalTime.dat
echo 'Total time is the biggest number out of these:' >> "${Project}"/totalTime.dat
grep 'real' "${Project}"/totalTime.dat | awk '{print $2}' >> "${Project}"/totalTime.dat
