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
	./generate.sh 4 2000 5000 4000 2
	./convert.sh s1
	./z_compile.sh
	./run.sh | grep "ref_id" | sort > tg
	./z_compile.sh GPU CMAT CBASES 64
	./run.sh 64 32 | grep "ref_id" | sort > ts
	diff tg ts
}

function debug {
	./z_compile.sh
	./run.sh > tg
	./z_compile.sh BATCH
	./run.sh 1 1 > tb
}

debug







