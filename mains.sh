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
export NUM_THREADS_END="128"
export NUM_THREADS_STEP="*=2"

# 2. With strong scaling (fixed batch size)
export BATCH_DELETIONS_BEGIN="0.001"
export BATCH_DELETIONS_END="0.001"
# ./main.sh "--strong-scaling"

# 3. With weak scaling
export BATCH_DELETIONS_BEGIN="0.0001"
export BATCH_DELETIONS_END="0.0128"
export BATCH_DELETIONS_STEP="*=2"
export NUM_THREADS_MODE="with-batch"
# ./main.sh "--weak-scaling"

# Signal completion
curl -X POST "https://maker.ifttt.com/trigger/puzzlef/with/key/${IFTTT_KEY}?value1=$src"
