#!/bin/bash

# fixInputVar.sh <path1>
#Arguments:
#	path1:	</path/to/INPUT>

inputPATH=$1

if [ -f "${inputPATH}" ]; then
	source "${inputPATH}"
else
	echo "NEED PATH TO INPUT FILE"
	exit
fi

if grep -q Input_Field "$HOME/.bashrc"; then
	grep -v "Input_Field" "$HOME/.bashrc" > /tmp/.bashrc
	mv /tmp/.bashrc "$HOME/.bashrc"
	echo "export Input_Field="${inputPATH}"" >> $HOME/.bashrc
else
	echo "export Input_Field="${inputPATH}"" >> $HOME/.bashrc
fi
