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
	sleep 4h
	./z_compile.sh
	./run.sh 8 | grep "ref_id" | sort > tg
	./z_compile.sh GPU CMAT CBASES 64 STREAM
	./run.sh 8 48 64 | grep "ref_id" | sort > ts
	diff tg ts
	./run.sh 1 256 64 | grep "ref_id" | sort > ts
	diff tg ts
}

function debug {
	./z_compile.sh GPU CMAT CBASES 64 STREAM STABLE
	./run.sh 1 1 32 > tg
	./z_compile.sh GPU CMAT CBASES 64 STREAM
	./run.sh 1 1 32 > ts
}

debug







