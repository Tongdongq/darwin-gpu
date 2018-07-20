#!/bin/bash

# compile with different flags for different optimizations

rm reference_guided

set -e

if [[ -n "$options" ]]; then
	printf "Error \$options was already set\n"
	exit 1
fi

options=""
maxregcount=128

for var in "$@"
do
	case "$var" in
	"BATCH")
		options="$options -D BATCH";;
	"GPU")
		options="$options -D GPU";;
	"CMAT")
		options="$options -D COALESCE_MATRICES";;
	"CBASES")
		options="$options -D COALESCE_BASES";;
	"TIME")
		options="$options -D TIME";;
	"STREAM")
		options="$options -D STREAM";;
	"CDIR")
		options="$options -D COMPRESS_DIR";;
	"GASAL")
		options="$options -D GASAL";;
	"STABLE")
		options="$options -D STABLE";;
	"32")
		maxregcount=32;;
	"64")
		maxregcount=64;;
	"128")
		maxregcount=128;;
	"256")
		maxregcount=256;;
	*)
		printf "Error unkown option '$var'\n"
		exit 1
	esac
done

if [ "$TACC" -eq "1" ]; then
	options="$options -D TACC=1"
	printf "Running on TACC\n"
else
	printf "Running on ce-cuda\n"
fi

old_options=$options
options="$options -D Z_COMPILE_USED"
gpu_options="$options --maxrregcount=$maxregcount"

# make 'options' visible to the makefile, note that 'echo $options' after the script is empty, because the script cannot set a variable in the parent environment
export options=$options
export gpu_options=$gpu_options

# compile
make clean
make

# let user know what options he used
printf "\nOptions used:\n$old_options\n"
head -n14 Makefile | grep "NCFLAGS" | grep -v "#"


