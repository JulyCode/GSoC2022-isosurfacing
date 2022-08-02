#!/bin/bash -l

./run_benchmark_threads.sh implicit_sphere 6 > threads_implicit.csv
./run_benchmark_threads.sh grid_sphere 6 > threads_grid.csv

./run_benchmark_size.sh implicit_sphere 1 500 > size_implicit.csv
./run_benchmark_size.sh grid_sphere 1 500 > size_grid.csv