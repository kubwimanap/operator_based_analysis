#!/bin/bash
# Reproduces every table and figure of the revised manuscript.
cd "$(dirname "$0")"
( export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  for e in sharp cert sens floor computed bench sens_alg trace; do python3 exps.py $e --seeds 30 > results/log_$e.txt 2>&1; done ) &
( export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  for e in baselines shift accuracy; do python3 exps2.py $e --seeds 30 > results/log2_$e.txt 2>&1; done ) &
wait
python3 exps2.py time2 --seeds 10 > results/log2_time2.txt 2>&1   # timings: run alone, all threads
python3 figures.py
echo DONE > results/DONE
