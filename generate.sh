#!/bin/bash

TYPE=${1:-1}
LENGTH=${2:-1100}
N=${3:-10000}
G_SIZE=${4:-1}
CROSSTALK=${5:-88}
THRES=${6:-5}

if [ "$TACC" -eq "1" ]; then
	src="$DATA/DAZZ_DB-master/human"
	printf "Running on TACC\n"
else
	src="../DAZZ_DB-master/human"
	printf "Running on ce-cuda\n"
fi

cd $src
./generate $TYPE $LENGTH $N $G_SIZE $CROSSTALK $THRES

