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
  auto glog = [&](const auto& ans, const char *technique, int numThreads, const auto& y, auto M, auto deletionsf, auto insertionsf) {
    printf(
      "{-%.3e/+%.3e batchf, %03d threads} -> "
      "{%09.1fms, %09.1fms preproc, %09.1fms firstpass, %09.1fms locmove, %09.1fms aggr, %zu/%zu affected, %04d iters, %03d passes, %01.9f modularity} %s\n",
      double(deletionsf), double(insertionsf), numThreads,
      ans.time, ans.preprocessingTime, ans.firstPassTime, ans.localMoveTime, ans.aggregationTime,
      ans.affectedVertices, y.order(), ans.iterations, ans.passes, getModularity(y, ans, M), technique
    );
  };
  // Get community memberships on original graph (static).
  auto b0 = louvainStaticOmp(x, init, {5});
  glog(b0, "louvainStaticOmpOriginal", MAX_THREADS, x, M, 0.0, 0.0);
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
        glog(ans, technique, numThreads, y, M, deletionsf, insertionsf);
      };
      // LouvainOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double aggregationTolerance=0.8, double toleranceDrop=10, int maxIterations=20, int maxPasses=10) :
      // Find static Louvain.
      auto b1 = louvainStaticOmp(y, init, {repeat});
      flog(b1, "louvainStaticOmp");
      // Find naive-dynamic Louvain.
      auto b20 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-2, 0.8, 10, 20, 10});
      flog(b20, "louvainNaiveDynamicOmp1e_2");
      auto b21 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-4, 0.8, 10, 20, 10});
      flog(b21, "louvainNaiveDynamicOmp1e_4");
      auto b22 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-6, 0.8, 10, 20, 10});
      flog(b22, "louvainNaiveDynamicOmp1e_6");
      auto b23 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-8, 0.8, 10, 20, 10});
      flog(b23, "louvainNaiveDynamicOmp1e_8");
      auto b24 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-10, 0.8, 10, 20, 10});
      flog(b24, "louvainNaiveDynamicOmp1e_10");
      auto b25 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-12, 0.8, 10, 20, 10});
      flog(b25, "louvainNaiveDynamicOmp1e_12");
      auto b26 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-14, 0.8, 10, 20, 10});
      flog(b26, "louvainNaiveDynamicOmp1e_14");
      auto b27 = louvainStaticOmp(y, &B2, {repeat, 1.0, 1e-16, 0.8, 10, 20, 10});
      flog(b27, "louvainNaiveDynamicOmp1e_16");
      // Find frontier based dynamic Louvain.
      auto b40 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-2, 0.8, 10, 20, 10});
      flog(b40, "louvainDynamicFrontierOmp1e_2");
      auto b41 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-4, 0.8, 10, 20, 10});
      flog(b41, "louvainDynamicFrontierOmp1e_4");
      auto b42 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-6, 0.8, 10, 20, 10});
      flog(b42, "louvainDynamicFrontierOmp1e_6");
      auto b43 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-8, 0.8, 10, 20, 10});
      flog(b43, "louvainDynamicFrontierOmp1e_8");
      auto b44 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-10, 0.8, 10, 20, 10});
      flog(b44, "louvainDynamicFrontierOmp1e_10");
      auto b45 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-12, 0.8, 10, 20, 10});
      flog(b45, "louvainDynamicFrontierOmp1e_12");
      auto b46 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-14, 0.8, 10, 20, 10});
      flog(b46, "louvainDynamicFrontierOmp1e_14");
      auto b47 = louvainDynamicFrontierOmp(y, deletions, insertions, &B4, {repeat, 1.0, 1e-16, 0.8, 10, 20, 10});
      flog(b47, "louvainDynamicFrontierOmp1e_16");
      // Find delta-screening based dynamic Louvain.
      auto b30 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-2, 0.8, 10, 20, 10});
      flog(b30, "louvainDynamicDeltaScreeningOmp1e_2");
      auto b31 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-4, 0.8, 10, 20, 10});
      flog(b31, "louvainDynamicDeltaScreeningOmp1e_4");
      auto b32 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-6, 0.8, 10, 20, 10});
      flog(b32, "louvainDynamicDeltaScreeningOmp1e_6");
      auto b33 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-8, 0.8, 10, 20, 10});
      flog(b33, "louvainDynamicDeltaScreeningOmp1e_8");
      auto b34 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-10, 0.8, 10, 20, 10});
      flog(b34, "louvainDynamicDeltaScreeningOmp1e_10");
      auto b35 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-12, 0.8, 10, 20, 10});
      flog(b35, "louvainDynamicDeltaScreeningOmp1e_12");
      auto b36 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-14, 0.8, 10, 20, 10});
      flog(b36, "louvainDynamicDeltaScreeningOmp1e_14");
      auto b37 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, &B3, {repeat, 1.0, 1e-16, 0.8, 10, 20, 10});
      flog(b37, "louvainDynamicDeltaScreeningOmp1e_16");
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
