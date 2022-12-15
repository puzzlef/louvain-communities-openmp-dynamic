#include <utility>
#include <random>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class K>
double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityByOmp(x, fc, M, 1);
}


template <class G, class R, class V>
auto addRandomEdges(G& a, R& rnd, size_t span, V w, int batchSize) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K, V>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v, w);
    a.addEdge(v, u, w);
    insertions.push_back(make_tuple(u, v, w));
    insertions.push_back(make_tuple(v, u, w));
  };
  for (int i=0; i<batchSize; ++i)
    retry([&]() { return addRandomEdge(a, rnd, span, w, fe); }, retries);
  updateOmpU(a);
  return insertions;
}


template <class G, class R>
auto removeRandomEdges(G& a, R& rnd, int batchSize) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    a.removeEdge(v, u);
    deletions.push_back(make_tuple(u, v));
    deletions.push_back(make_tuple(v, u));
  };
  for (int i=0; i<batchSize; ++i)
    retry([&]() { return removeRandomEdge(a, rnd, fe); }, retries);
  updateOmpU(a);
  return deletions;
}




template <class G>
void runExperiment(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int retries = 5;
  double M = edgeWeightOmp(x)/2;
  double Q = modularityOmp(x, M, 1);
  LOG("[%01.6f modularity] noop\n", Q);

  // Get community memberships on original graph (static).
  auto ak = louvainSeqStatic(x, init);
  // Batch of additions only (dynamic).
  // Remove sequential algorithms to reduce duration of experiment.
  for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
    for (int batchCount=1; batchCount<=5; ++batchCount) {
      auto   y = duplicate(x);
      auto insertions = addRandomEdges(y, rnd, x.span(), V(1), batchSize); vector<tuple<K, K>> deletions;
      double M = edgeWeightOmp(y)/2;
      // Find static Louvain (sequential).
      auto al = louvainSeqStatic(y, init, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqStatic\n",                double(batchSize), al.time, al.iterations, al.passes, getModularity(y, al, M));
      // Find static Louvain.
      auto am = louvainOmpStatic(y, init, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpStatic\n",                double(batchSize), am.time, am.iterations, am.passes, getModularity(y, am, M));
      // Find naive-dynamic Louvain.
      auto an = louvainOmpStatic(y, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpNaiveDynamic\n",          double(batchSize), an.time, an.iterations, an.passes, getModularity(y, an, M));
      // Find delta-screening based dynamic Louvain.
      auto ao = louvainOmpDynamicDeltaScreening(y, deletions, insertions, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicDeltaScreening\n", double(batchSize), ao.time, ao.iterations, ao.passes, getModularity(y, ao, M));
      // Find frontier based dynamic Louvain.
      auto ap = louvainOmpDynamicFrontier(y, deletions, insertions, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicFrontier\n",       double(batchSize), ap.time, ap.iterations, ap.passes, getModularity(y, ap, M));
    }
  }
  // Batch of deletions only (dynamic).
  // for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
  //   for (int batchCount=1; batchCount<=5; ++batchCount) {
  //     auto   y = duplicate(x);
  //     auto deletions = removeRandomEdges(y, rnd, batchSize); vector<tuple<K, K, V>> insertions;
  //     double M = edgeWeight(y)/2;
  //     // Find static Louvain (sequential).
  //     auto al = louvainSeqStatic(y, init, {repeat});
  //     printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqStatic\n",                double(-batchSize), al.time, al.iterations, al.passes, getModularity(y, al, M));
  //     // Find static Louvain.
  //     auto am = louvainOmpStatic(y, init, {repeat});
  //     printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpStatic\n",                double(-batchSize), am.time, am.iterations, am.passes, getModularity(y, am, M));
  //     // Find naive-dynamic Louvain.
  //     auto an = louvainOmpStatic(y, &ak.membership, {repeat});
  //     printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpNaiveDynamic\n",          double(-batchSize), an.time, an.iterations, an.passes, getModularity(y, an, M));
  //     // Find delta-screening based dynamic Louvain.
  //     auto ao = louvainOmpDynamicDeltaScreening(y, deletions, insertions, &ak.membership, {repeat});
  //     printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicDeltaScreening\n", double(-batchSize), ao.time, ao.iterations, ao.passes, getModularity(y, ao, M));
  //     // Find frontier based dynamic Louvain.
  //     auto ap = louvainOmpDynamicFrontier(y, deletions, insertions, &ak.membership, {repeat});
  //     printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicFrontier\n",       double(-batchSize), ap.time, ap.iterations, ap.passes, getModularity(y, ap, M));
  //   }
  // }
}


int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  install_sigsegv();
  char *file     = argv[1];
  bool symmetric = argc>2? stoi(argv[2]) : false;
  bool weighted  = argc>3? stoi(argv[3]) : false;
  int  repeat    = argc>4? stoi(argv[4]) : 5;
  OutDiGraph<K, None, V> x;  // V w = 1;
  LOG("Loading graph %s ...\n", file);
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { x = symmetricizeOmp(x); LOG(""); print(x); printf(" (symmetricize)\n"); }
  // auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runExperiment(x, repeat);
  printf("\n");
  return 0;
}
