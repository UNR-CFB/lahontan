#!/bin/bash

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Making Cols.dat
###############################################################

cd "${Project}"

echo condition > Cols.dat

for name in $(head -n +1 Counts.dat); do
	echo "${name}	blah" >> Cols.dat
done
