#!/bin/bash -l

if [ $# -lt 3 ] ; then
    echo "Arguments 'scenario', 'num_threads' and 'max_grid_size' required!"
    exit 1
fi

SCENARIO=$1
THREADS=$2
MAX_N=$3

MAX_CELLS=$(echo "$MAX_N^3" | bc)
i=1000

export OMP_NUM_THREADS=$THREADS

echo "cells;threads;time"

while [[ $i -le $MAX_CELLS ]]
do
    N=$(echo "$i" | awk '{ print int($1^(1.0/3.0)); }')
    echo -n "$i;$THREADS;"
    ../../build/Benchmark $SCENARIO $N
    i=$(echo "($i*1.2)/1" | bc)
done
