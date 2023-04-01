#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "src/main.hxx"

using namespace std;




// Fixed config
#ifndef TYPE
#define TYPE float
#endif
#ifndef MAX_THREADS
#define MAX_THREADS 64
#endif
#ifndef REPEAT_BATCH
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
#define REPEAT_METHOD 1
#endif




// HELPERS
// -------

template <class G, class K>
inline double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityByOmp(x, fc, M, 1.0);
}




// PERFORM EXPERIMENT
// ------------------

template <class G>
void runExperiment(const G& x) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  auto fm = [](auto v)         { return v; };
  auto fn = [](auto k, auto v) { return v; };
  vector<K> *init = nullptr;
  // Get community memberships on original graph (static).
  auto b0   = louvainStaticOmp(x, init);
  auto M    = edgeWeightOmp(x)/2;
  auto Q    = getModularity(x, b0, M);
  auto csiz = countAs(x, b0.membership, fm);
  auto sizs = valuesOf(csiz, fn);
  auto lc   = lorenzCurve(sizs);
  auto gc   = giniCoefficient(lc);
  LOG(
    "{%09.1fms, %04d iters, %03d passes, %01.9f modularity, %01.9f gini, %zu communities} %s\n",
    b0.time, b0.iterations, b0.passes, Q, gc, csiz.size(), "louvainStaticOmp"
  );
}


int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  install_sigsegv();
  char *file     = argv[1];
  bool symmetric = argc>2? stoi(argv[2]) : false;
  bool weighted  = argc>3? stoi(argv[3]) : false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  OutDiGraph<K, None, V> x;
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { x = symmetricizeOmp(x); LOG(""); print(x); printf(" (symmetricize)\n"); }
  runExperiment(x);
  printf("\n");
  return 0;
}
