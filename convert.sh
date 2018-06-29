#!/bin/bash

# copy fasta files

# possible values: 1, 2, 6_8MB, 50MB, 130MB, 225MB, s1

if [ "$TACC" -eq "1" ]; then
	src="$DATA/DAZZ_DB-master/human"
	printf "Running on TACC\n"
else
	src="../DAZZ_DB-master/human"
	printf "Running on ce-cuda\n"
fi

cp $src/s1.1.subreads.fasta reference.fasta
cp $src/s1.2.subreads.fasta reads.fasta

