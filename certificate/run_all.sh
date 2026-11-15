#!/bin/bash
cd "$(dirname "$0")"
( export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  for e in sharp cert sens floor computed; do python3 exps.py $e --seeds 30 > results/log_$e.txt 2>&1; done ) &
( export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  python3 exps.py lazy_counts --seeds 30 > results/log_lazy_counts.txt 2>&1 ) &
wait
python3 exps.py lazy_time --seeds 5 > results/log_lazy_time.txt 2>&1
echo DONE > results/DONE
