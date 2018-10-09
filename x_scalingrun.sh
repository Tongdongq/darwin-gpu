#!/bin/bash

function bench_and_compare {
	sleep 5h
	date
	./generate.sh 5 10000 10000 20000 1
	./convert.sh s1
	./z_compile.sh
	date
}

function compare {
	./z_compile.sh GPU GASAL 64 STREAM CMAT
	./run.sh 8 32 64
	cat darwin.*.out | sort > tss
	diff tss tg_50MB | head
}

function debug {
	./z_compile.sh GPU CMAT CBASES 64 STREAM
	./run.sh 1 1 32 > tg
	./z_compile.sh GPU CMAT CBASES 64 STREAM
	./run.sh 1 1 32 > ts
}

function profile {
	sleep 4h
	./profile.sh f 1 256 64
}

#compare







