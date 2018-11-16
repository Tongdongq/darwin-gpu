#!/bin/bash

# executes gact

#batch size
T=${1:-1}
BLOCKS=${2:-1}
TPB=${3:-3}

#time ./reference_guided sacCer3.fa reads.fa $T $BLOCKS $TPB
#time ./reference_guided reference.fasta reads.fasta $T $BLOCKS $TPB
time ./reference_guided reference.fasta reference.fasta $T $BLOCKS $TPB

last="$(head -n20 Makefile | grep "NCFLAGS" | grep -v "#")"
if [[ $last == *"-O0"* ]]; then
	printf "\nNote: probably compiled with -O0\n\n"
fi
if [[ $last == *"-O3"* ]]; then
	printf "\nNote: probably compiled with -O3\n\n"
fi
if [[ $last == *"-lineinfo"* ]]; then
	printf "\nNote: probably compiled with -lineinfo\n\n"
fi
last="$(head -n20 Makefile | grep "CFLAGS" | grep -v "#" | tail -n 1)"
if [[ $last == *"-O0"* ]]; then
	printf "\nNote: probably compiled with -O0\n\n"
fi
if [[ $last == *"-O3"* ]]; then
	printf "\nNote: probably compiled with -O3\n\n"
fi
if [[ $last == *"-lineinfo"* ]]; then
	printf "\nNote: probably compiled with -lineinfo\n\n"
fi

