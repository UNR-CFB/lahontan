#!/bin/bash

# Arguments:
#	<num of samples>	makes sample directories

if [ -f "${Input_Field}" ]; then
	source "${Input_Field}"
else
	echo "NEED INPUT FILE"
	exit
fi

###############################################################
# Structure Maker
###############################################################

numSamp=$1

# Test argument validity
re='^[0-9]+$'
if ! [[ "${numSamp}" =~ $re ]] ; then
   echo "Need number of samples as an argument" >&2; exit 1
fi

if [ ! -d "${Project}" ]; then
	mkdir "${Project}"
fi

if [ ! -d "${Data}" ]; then
	mkdir "${Data}"
fi

if [ ! -d "${Ref}" ]; then
	mkdir "${Ref}"
fi

for number in `seq -w 1 "${numSamp}"`; do
	mkdir "${Data}"/sample_$number
done
