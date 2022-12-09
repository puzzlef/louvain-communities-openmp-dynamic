#!/usr/bin/env bash
src="louvain-openmp-static-vs-dynamic"
out="/home/resources/Documents/subhajit/$src.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
rm -rf $src
git clone https://github.com/puzzlef/$src
cd $src

# Run
g++ -std=c++17 -O3 -fopenmp main.cxx
stdbuf --output=L ./a.out ~/Data/GAP-road.mtx       1 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-twitter.mtx    0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-web.mtx        0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-kron.mtx       1 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-urand.mtx      1 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Orkut.mtx      1 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Friendster.mtx 1 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/arabic-2005.mtx    0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2005.mtx        0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/webbase-2001.mtx   0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/it-2004.mtx        0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/sk-2005.mtx        0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_V1r.mtx       1 2>&1 | tee -a "$out"
