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




template <class G, class K, class V>
double getModularity(const G& x, const LouvainResult<K>& a, V M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, V(1));
}


template <class G, class R, class K, class V>
auto addRandomEdges(G& a, R& rnd, K span, V w, int batchSize) {
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
  a.correct();
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
  a.correct();
  return deletions;
}




template <class G>
void runLouvain(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int retries  = 5;
  auto M = edgeWeight(x)/2;
  auto Q = modularity(x, M, 1.0f);
  printf("[%01.6f modularity] noop\n", Q);

  // Get community memberships on original graph (static).
  auto ak = louvainSeqStatic(x, init, {repeat});
  printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqStatic\n", 0.0, ak.time, ak.iterations, ak.passes, getModularity(x, ak, M));
  // Batch of additions only (dynamic).
  for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
    for (int batchCount=1; batchCount<=5; ++batchCount) {
      auto y = duplicate(x);
      auto insertions = addRandomEdges(y, rnd, x.span(), V(1), batchSize); vector<tuple<K, K>> deletions;
      auto M = edgeWeight(y)/2;
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
  for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
    for (int batchCount=1; batchCount<=5; ++batchCount) {
      auto y = duplicate(x);
      auto deletions = removeRandomEdges(y, rnd, batchSize); vector<tuple<K, K, V>> insertions;
      auto M = edgeWeight(y)/2;
      // Find static Louvain (sequential).
      auto al = louvainSeqStatic(y, init, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqStatic\n",                double(-batchSize), al.time, al.iterations, al.passes, getModularity(y, al, M));
      // Find static Louvain.
      auto am = louvainOmpStatic(y, init, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpStatic\n",                double(-batchSize), am.time, am.iterations, am.passes, getModularity(y, am, M));
      // Find naive-dynamic Louvain.
      auto an = louvainOmpStatic(y, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpNaiveDynamic\n",          double(-batchSize), an.time, an.iterations, an.passes, getModularity(y, an, M));
      // Find delta-screening based dynamic Louvain.
      auto ao = louvainOmpDynamicDeltaScreening(y, deletions, insertions, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicDeltaScreening\n", double(-batchSize), ao.time, ao.iterations, ao.passes, getModularity(y, ao, M));
      // Find frontier based dynamic Louvain.
      auto ap = louvainOmpDynamicFrontier(y, deletions, insertions, &ak.membership, {repeat});
      printf("[%1.0e batch_size; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainOmpDynamicFrontier\n",       double(-batchSize), ap.time, ap.iterations, ap.passes, getModularity(y, ap, M));
    }
  }
}


int main(int argc, char **argv) {
  using K = int;
  using V = TYPE;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x; V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW(x, file); println(x);
  auto y = symmetricize(x); print(y); printf(" (symmetricize)\n");
  // auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runLouvain(y, repeat);
  printf("\n");
  return 0;
}
