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
inline double matchRatio(const G& x, const vector<K>& ans, const vector<K>& ref) {
  size_t n = 0;
  x.forEachVertexKey([&](auto u) {
    if (ans[u]==ref[u]) ++n;
  });
  return double(n) / x.order();
}


template <class G, class K>
inline double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, 1.0);
}




// GENERATE BATCH
// --------------

template <class G, class R>
inline auto addRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  int retries = 5;
  vector<tuple<K, K, V>> insertions;
  auto fe = [&](auto u, auto v, auto w) {
    a.addEdge(u, v, w);
    a.addEdge(v, u, w);
    insertions.push_back(make_tuple(u, v, w));
    insertions.push_back(make_tuple(v, u, w));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return addRandomEdge(a, rnd, i, n, V(1), fe); }, retries);
  updateOmpU(a);
  return insertions;
}


template <class G, class R>
auto removeRandomEdges(G& a, R& rnd, size_t batchSize, size_t i, size_t n) {
  using K = typename G::key_type;
  int retries = 5;
  vector<tuple<K, K>> deletions;
  auto fe = [&](auto u, auto v) {
    a.removeEdge(u, v);
    a.removeEdge(v, u);
    deletions.push_back(make_tuple(u, v));
    deletions.push_back(make_tuple(v, u));
    return true;
  };
  for (size_t l=0; l<batchSize; ++l)
    retry([&]() { return removeRandomEdge(a, rnd, i, n, fe); }, retries);
  updateOmpU(a);
  return deletions;
}


template <class K, class V>
inline auto withoutWeights(const vector<tuple<K, K, V>>& edges) {
  vector<tuple<K, K>> a;
  for (auto [u, v, w] : edges)
    a.push_back(make_tuple(u, v));
  return a;
}


template <class K, class V>
inline auto withWeights(const vector<tuple<K, K>>& edges, V w) {
  vector<tuple<K, K, V>> a;
  for (auto [u, v] : edges)
    a.push_back(make_tuple(u, v, w));
  return a;
}




// PERFORM EXPERIMENT
// ------------------

template <class G, class R, class F>
inline void runAbsoluteBatches(const G& x, R& rnd, F fn) {
  size_t d = BATCH_DELETIONS_BEGIN;
  size_t i = BATCH_INSERTIONS_BEGIN;
  for (int epoch=0;; ++epoch) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      for (int sequence=0; sequence<BATCH_LENGTH; ++sequence) {
      auto deletions  = removeRandomEdges(y, rnd, d, 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, i, 1, x.span()-1);
        fn(y, d, deletions, i, insertions, sequence, epoch);
      }
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, size_t(BATCH_DELETIONS_END));
    i = min(i, size_t(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runRelativeBatches(const G& x, R& rnd, F fn) {
  double d = BATCH_DELETIONS_BEGIN;
  double i = BATCH_INSERTIONS_BEGIN;
  for (int epoch=0;; ++epoch) {
    for (int r=0; r<REPEAT_BATCH; ++r) {
      auto y  = duplicate(x);
      for (int sequence=0; sequence<BATCH_LENGTH; ++sequence) {
      auto deletions  = removeRandomEdges(y, rnd, size_t(d * x.size()/2), 1, x.span()-1);
      auto insertions = addRandomEdges   (y, rnd, size_t(i * x.size()/2), 1, x.span()-1);
        fn(y, d, deletions, i, insertions, sequence, epoch);
      }
    }
    if (d>=BATCH_DELETIONS_END && i>=BATCH_INSERTIONS_END) break;
    d BATCH_DELETIONS_STEP;
    i BATCH_INSERTIONS_STEP;
    d = min(d, double(BATCH_DELETIONS_END));
    i = min(i, double(BATCH_INSERTIONS_END));
  }
}


template <class G, class R, class F>
inline void runBatches(const G& x, R& rnd, F fn) {
  if (BATCH_UNIT=="%") runRelativeBatches(x, rnd, fn);
  else runAbsoluteBatches(x, rnd, fn);
}


template <class F>
inline void runThreadsWithBatch(int epoch, F fn) {
  int t = NUM_THREADS_BEGIN;
  for (int l=0; l<epoch && t<=NUM_THREADS_END; ++l)
    t NUM_THREADS_STEP;
  omp_set_num_threads(t);
  fn(t);
  omp_set_num_threads(MAX_THREADS);
}


template <class F>
inline void runThreadsAll(F fn) {
  for (int t=NUM_THREADS_BEGIN; t<=NUM_THREADS_END; t NUM_THREADS_STEP) {
    omp_set_num_threads(t);
    fn(t);
    omp_set_num_threads(MAX_THREADS);
  }
}


template <class F>
inline void runThreads(int epoch, F fn) {
  if (NUM_THREADS_MODE=="with-batch") runThreadsWithBatch(epoch, fn);
  else runThreadsAll(fn);
}


template <class G>
void runExperiment(const G& x) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat  = REPEAT_METHOD;
  int retries = 5;
  vector<K> *init = nullptr;
  double M = edgeWeightOmp(x)/2;
  // Follow a specific result logging format, which can be easily parsed later.
  auto glog = [&](const auto& ans, const auto& ref, const char *technique, int numThreads, const auto& y, auto M, auto deletionsf, auto insertionsf) {
    printf(
      "{-%.3e/+%.3e batchf, %03d threads} -> "
      "{%09.1fms, %09.1fms preproc, %09.1fms firstpass, %09.1fms locmove, %09.1fms aggr, %.3e affected, %04d iters, %03d passes, %01.9f modularity, %01.9f match} %s\n",
      double(deletionsf), double(insertionsf), numThreads,
      ans.time, ans.preprocessingTime, ans.firstPassTime, ans.localMoveTime, ans.aggregationTime,
      double(ans.affectedVertices), ans.iterations, ans.passes, getModularity(y, ans, M),
      matchRatio(x, ans.membership, ref.membership), technique
    );
  };
  // Get community memberships on original graph (static).
  auto b0 = louvainStaticOmp(x, init, {5});
  glog(b0, b0, "louvainStaticOmpOriginal", MAX_THREADS, x, M, 0.0, 0.0);
  #if BATCH_LENGTH>1
  vector<K> B2, B3, B4;
  #else
  const auto& B2 = b0.membership;
  const auto& B3 = b0.membership;
  const auto& B4 = b0.membership;
  #endif
  // Get community memberships on updated graph (dynamic).
  runBatches(x, rnd, [&](const auto& y, auto deletionsf, const auto& deletions, auto insertionsf, const auto& insertions, int sequence, int epoch) {
    double M = edgeWeightOmp(y)/2;
    #if BATCH_LENGTH>1
    if (sequence==0) {
      B2 = b0.membership;
      B3 = b0.membership;
      B4 = b0.membership;
    }
    #endif
    // Adjust number of threads.
    runThreads(epoch, [&](int numThreads) {
      auto flog = [&](const auto& ans, const char *technique) {
        glog(ans, b0, technique, numThreads, y, M, deletionsf, insertionsf);
      };
      auto rdeletions  = withoutWeights(insertions);
      auto rinsertions = withWeights(deletions, V(1));
      // Find static Louvain.
      auto b1d = louvainStaticOmp(y, init, {repeat});
      auto b1i = louvainStaticOmp(x, init, {repeat});
      flog(b1i, "louvainStaticOmp");
      // Find naive-dynamic Louvain.
      auto b2d = louvainStaticOmp(y, &b0 .membership, {repeat});
      auto b2i = louvainStaticOmp(x, &b2d.membership, {repeat});
      flog(b2i, "louvainNaiveDynamicOmp");
      // Find frontier based dynamic Louvain.
      auto b4d = louvainDynamicFrontierOmp(y, deletions,  insertions,  &b0 .membership, {repeat});
      auto b4i = louvainDynamicFrontierOmp(x, rdeletions, rinsertions, &b4d.membership, {repeat});
      flog(b4i, "louvainDynamicFrontierOmp");
      // Find delta-screening based dynamic Louvain.
      auto b3d = louvainDynamicDeltaScreeningOmp(y, deletions,  insertions,  &b0 .membership, {repeat});
      auto b3i = louvainDynamicDeltaScreeningOmp(x, rdeletions, rinsertions, &b3d.membership, {repeat});
      flog(b3i, "louvainDynamicDeltaScreeningOmp");
      #if BATCH_LENGTH>1
      B2 = b2.membership;
      B3 = b3.membership;
      B4 = b4.membership;
      #endif
    });
  });
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
  DiGraph<K, None, V> x;
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { x = symmetricizeOmp(x); LOG(""); print(x); printf(" (symmetricize)\n"); }
  runExperiment(x);
  printf("\n");
  return 0;
}
