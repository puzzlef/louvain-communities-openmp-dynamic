#!/usr/bin/env bash
src="louvain-communities-openmp-dynamic"
ulimit -s unlimited

# Download program
if [[ "$DOWNLOAD" != "0" ]]; then
  rm -rf $src
  git clone https://github.com/puzzlef/$src
  cd $src
fi

# Don't need to download program again.
export DOWNLOAD="0"

# 1. Static vs Dynamic Louvain
export MAX_THREADS="64"
./main.sh

# For scaling experiments
export NUM_THREADS_BEGIN="1"
export NUM_THREADS_END="$MAX_THREADS"
export NUM_THREADS_STEP="*=2"

# 2. With strong scaling (fixed batch size)
export BATCH_DELETIONS_BEGIN="0.0005"
export BATCH_DELETIONS_END="0.0005"
export BATCH_INSERTIONS_BEGIN="0.0005"
export BATCH_INSERTIONS_END="0.0005"
# ./main.sh "--strong-scaling"

# 3. With weak scaling
export BATCH_DELETIONS_BEGIN="0.00005"
export BATCH_DELETIONS_END="0.0032"
export BATCH_DELETIONS_STEP="*=2"
export BATCH_INSERTIONS_BEGIN="0.00005"
export BATCH_INSERTIONS_END="0.0032"
export BATCH_INSERTIONS_STEP="*=2"
export NUM_THREADS_MODE="with-batch"
# ./main.sh "--weak-scaling"
