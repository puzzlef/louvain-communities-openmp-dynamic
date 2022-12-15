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
stdbuf --output=L ./a.out ~/Data/web-Stanford.mtx      0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-BerkStan.mtx      0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-Google.mtx        0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/web-NotreDame.mtx     0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Slashdot0811.mtx  0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Slashdot0902.mtx  0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-Epinions1.mtx     0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coAuthorsDBLP.mtx     1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coAuthorsCiteseer.mtx 1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/soc-LiveJournal1.mtx  0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coPapersCiteseer.mtx  1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/coPapersDBLP.mtx      1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/indochina-2004.mtx    0 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/italy_osm.mtx         1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/great-britain_osm.mtx 1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/germany_osm.mtx       1 0 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/asia_osm.mtx          1 0 2>&1 | tee -a "$out"

# stdbuf --output=L ./a.out ~/Data/GAP-road.mtx       1 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/GAP-twitter.mtx    0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/GAP-web.mtx        0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/GAP-kron.mtx       1 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/GAP-urand.mtx      1 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/com-Orkut.mtx      1 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/com-Friendster.mtx 1 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/arabic-2005.mtx    0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/uk-2005.mtx        0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/webbase-2001.mtx   0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/it-2004.mtx        0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/sk-2005.mtx        0 0 2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out ~/Data/kmer_V1r.mtx       1 0 2>&1 | tee -a "$out"
