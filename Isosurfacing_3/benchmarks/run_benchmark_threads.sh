#!/bin/bash -l

if [ $# -lt 2 ] ; then
    echo "Arguments 'scenario' and 'max_threads' required!"
    exit 1
fi

N=1000
SCENARIO=$1
MAX_THREADS=$2

CELLS=$(echo "$N^3" | bc)

echo "cells;threads;time"

for (( i=1; i<=$MAX_THREADS; i++ ))
do
    echo -n "$CELLS;$i;"
    export OMP_NUM_THREADS=$i
    ../../build/Benchmark $SCENARIO $N
done
