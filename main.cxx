#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "inc/main.hxx"

using namespace std;




#pragma region CONFIGURATION
#ifndef TYPE
/** Type of edge weights. */
#define TYPE float
#endif
#ifndef MAX_THREADS
/** Maximum number of threads to use. */
#define MAX_THREADS 64
#endif
#ifndef REPEAT_BATCH
/** Number of times to repeat each batch. */
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
/** Number of times to repeat each method. */
#define REPEAT_METHOD 5
#endif
#pragma endregion




#pragma region METHODS
#pragma region HELPERS
/**
 * Obtain the modularity of community structure on a graph.
 * @param x original graph
 * @param a rak result
 * @param M sum of edge weights
 * @returns modularity
 */
template <class G, class K>
inline double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, 1.0);
}
#pragma endregion




#pragma region EXPERIMENTAL SETUP
/**
 * Run a function on each batch update, with a specified range of batch sizes.
 * @param x original graph
 * @param rnd random number generator
 * @param fn function to run on each batch update
 */
template <class G, class R, class F>
inline void runBatches(const G& x, R& rnd, F fn) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  random_device dev;
  default_random_engine rng(dev());
  uniform_int_distribution<long long int> dist(1, 100000000);
  // Batch update with insertions to merge groups 1 and 2 with cascading updates.
  V w = V(1000000);
  int count = 1000;
  {
    auto y = duplicate(x);
    vector<tuple<K, K, V>> deletions;
    vector<tuple<K, K, V>> insertions;
    for (int i=0; i<count; ++i) {
      K u = K(sqrt(dist(rng)));
      K v = 10000 + K(sqrt(dist(rng)));
      insertions.push_back(make_tuple(u, v, w));
      insertions.push_back(make_tuple(v, u, w));
    }
    tidyBatchUpdateU(deletions, insertions, y);
    applyBatchUpdateOmpU(y, deletions, insertions);
    fn(y, 0.0, deletions, insertions.size()/double(x.size()), insertions, 0, 0);
  }
  // Batch update with to split group 3 with cascading updates.
  {
    auto y = duplicate(x);
    vector<tuple<K, K, V>> deletions;
    vector<tuple<K, K, V>> insertions;
    for (int i=0; i<count; ++i) {
      K u1 = K(sqrt(dist(rng)));
      K v1 = 20000 + K(sqrt(dist(rng))) / 2;
      K u2 = 10000 + K(sqrt(dist(rng)));
      K v2 = 25000 + K(sqrt(dist(rng))) / 2;
      insertions.push_back(make_tuple(u1, v1, w));
      insertions.push_back(make_tuple(v1, u1, w));
      insertions.push_back(make_tuple(u2, v2, w));
      insertions.push_back(make_tuple(v2, u2, w));
    }
    tidyBatchUpdateU(deletions, insertions, y);
    applyBatchUpdateOmpU(y, deletions, insertions);
    fn(y, 0.0, deletions, insertions.size()/double(x.size()), insertions, 0, 0);
  }
}


/**
 * Run a function on each number of threads, for a specific epoch.
 * @param epoch epoch number
 * @param fn function to run on each number of threads
 */
template <class F>
inline void runThreadsWithBatch(int epoch, F fn) {
  int t = NUM_THREADS_BEGIN;
  for (int l=0; l<epoch && t<=NUM_THREADS_END; ++l)
    t NUM_THREADS_STEP;
  omp_set_num_threads(t);
  fn(t);
  omp_set_num_threads(MAX_THREADS);
}


/**
 * Run a function on each number of threads, with a specified range of thread counts.
 * @param fn function to run on each number of threads
 */
template <class F>
inline void runThreadsAll(F fn) {
  for (int t=NUM_THREADS_BEGIN; t<=NUM_THREADS_END; t NUM_THREADS_STEP) {
    omp_set_num_threads(t);
    fn(t);
    omp_set_num_threads(MAX_THREADS);
  }
}


/**
 * Run a function on each number of threads, with a specified range of thread counts or for a specific epoch (depending on NUM_THREADS_MODE).
 * @param epoch epoch number
 * @param fn function to run on each number of threads
 */
template <class F>
inline void runThreads(int epoch, F fn) {
  if (NUM_THREADS_MODE=="with-batch") runThreadsWithBatch(epoch, fn);
  else runThreadsAll(fn);
}
#pragma endregion




#pragma region PERFORM EXPERIMENT
/**
 * Perform the experiment.
 * @param x original graph
 */
template <class G>
void runExperiment(const G& x) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  using W = LOUVAIN_WEIGHT_TYPE;
  random_device dev;
  default_random_engine rnd(dev());
  int repeat  = REPEAT_METHOD;
  double M = edgeWeightOmp(x)/2;
  // Follow a specific result logging format, which can be easily parsed later.
  auto glog = [&](const auto& ans, const char *technique, int numThreads, const auto& y, auto M, auto deletionsf, auto insertionsf) {
    printf(
      "{-%.3e/+%.3e batchf, %03d threads} -> "
      "{%09.1fms, %09.1fms mark, %09.1fms init, %09.1fms firstpass, %09.1fms locmove, %09.1fms aggr, %.3e aff, %04d iters, %03d passes, %01.9f modularity, %zu/%zu disconnected} %s\n",
      double(deletionsf), double(insertionsf), numThreads,
      ans.time, ans.markingTime, ans.initializationTime, ans.firstPassTime, ans.localMoveTime, ans.aggregationTime,
      double(ans.affectedVertices), ans.iterations, ans.passes, getModularity(y, ans, M),
      countValue(communitiesDisconnectedOmp(x, ans.membership), char(1)),
      communities(x, ans.membership).size(), technique
    );
  };
  // Get community memberships on original graph (static).
  auto b0 = louvainStaticOmp(x, {5});
  glog(b0, "louvainStaticOmpOriginal", MAX_THREADS, x, M, 0.0, 0.0);
  printf("\n");
  #if BATCH_LENGTH>1
  vector<K> B2, B3, B4;
  vector<W> VW, CW;
  #else
  const auto& B2 = b0.membership;
  const auto& B3 = b0.membership;
  const auto& B4 = b0.membership;
  const auto& VW = b0.vertexWeight;
  const auto& CW = b0.communityWeight;
  #endif
  // Get community memberships on updated graph (dynamic).
  runBatches(x, rnd, [&](const auto& y, auto deletionsf, const auto& deletions, auto insertionsf, const auto& insertions, int sequence, int epoch) {
    double M = edgeWeightOmp(y)/2;
    #if BATCH_LENGTH>1
    if (sequence==0) {
      B2 = b0.membership;
      B3 = b0.membership;
      B4 = b0.membership;
      VW = b0.vertexWeight;
      CW = b0.communityWeight;
    }
    #endif
    // Adjust number of threads.
    runThreads(epoch, [&](int numThreads) {
      auto flog = [&](const auto& ans, const char *technique) {
        glog(ans, technique, numThreads, y, M, deletionsf, insertionsf);
      };
      double resolution = 1;
      double dynamicTolerance = 1e-2 * min(deletionsf + insertionsf, 1.0);
      // Find static Louvain.
      auto b1 = louvainStaticOmp(y, {repeat});
      flog(b1, "louvainStaticOmp");
      // Find naive-dynamic Louvain.
      auto b21 = louvainNaiveDynamicOmp(y, deletions, insertions, B2, VW, CW, {repeat});
      flog(b21, "louvainNaiveDynamicOmp");
      auto b22 = louvainNaiveDynamicOmp(y, deletions, insertions, B2, VW, CW, {repeat, resolution, dynamicTolerance});
      flog(b22, "louvainNaiveDynamicOmpAdjusted");
      // Find delta-screening based dynamic Louvain.
      auto b31 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, B3, VW, CW, {repeat});
      flog(b31, "louvainDynamicDeltaScreeningOmp");
      auto b32 = louvainDynamicDeltaScreeningOmp(y, deletions, insertions, B3, VW, CW, {repeat, resolution, dynamicTolerance});
      flog(b32, "louvainDynamicDeltaScreeningOmpAdjusted");
      // Find frontier based dynamic Louvain.
      auto b41 = louvainDynamicFrontierOmp(y, deletions, insertions, B4, VW, CW, {repeat});
      flog(b41, "louvainDynamicFrontierOmp");
      auto b42 = louvainDynamicFrontierOmp(y, deletions, insertions, B4, VW, CW, {repeat, resolution, dynamicTolerance});
      flog(b42, "louvainDynamicFrontierOmpAdjusted");
      #if BATCH_LENGTH>1
      B2 = b2.membership;
      B3 = b3.membership;
      B4 = b4.membership;
      VW = b1.vertexWeight;
      CW = b1.communityWeight;
      #endif
    });
    printf("\n");
  });
}


/**
 * Create a sample graph with three groups of 100 vertices each (1 weakly connected).
 * @param degree average degree of vertices
 * @returns sample graph
 */
template <class K, class V>
DiGraph<K, None, V> createSampleGraph(int degree) {
  random_device dev;
  default_random_engine rng(dev());
  uniform_int_distribution<long long int> dist(1, 100000000);
  DiGraph<K, None, V> x;
  // Add all vertices.
  for (K u=1; u<=30000; ++u)
    x.addVertex(u);
  // Add edges within group 1.
  for (K u=1; u<=10000; ++u) {
    for (int i=0; i<degree; ++i) {
      K v = K(sqrt(dist(rng)));
      x.addEdge(u, v, 1);
    }
  }
  // Add edges within group 2.
  for (K u=10001; u<=20000; ++u) {
    for (int i=0; i<degree; ++i) {
      K v = 10000 + K(sqrt(dist(rng)));
      x.addEdge(u, v, 1);
    }
  }
  // Add edges within one half of group 3.
  for (K u=20001; u<=25000; ++u) {
    for (int i=0; i<degree; ++i) {
      K v = 20000 + K(sqrt(dist(rng))) / 2;
      x.addEdge(u, v, 1);
    }
  }
  // Add edges within the other half of group 3.
  for (K u=25001; u<=30000; ++u) {
    for (int i=0; i<degree; ++i) {
      K v = 25000 + K(sqrt(dist(rng))) / 2;
      x.addEdge(u, v, 1);
    }
  }
  // Add edges between the two halves of group 3.
  for (K i=1; i<4000*degree; ++i) {
    K u = 20000 + K(sqrt(dist(rng))) / 2;
    K v = 25000 + K(sqrt(dist(rng))) / 2;
    x.addEdge(u, v, 1);
  }
  // Update and symmetrize.
  updateOmpU(x);
  x = symmetrizeOmp(x);
  return x;
}


/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  install_sigsegv();
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Creating sample graph ...\n");
  DiGraph<K, None, V> x = createSampleGraph<K, V>(10);
  LOG(""); println(x);
  runExperiment(x);
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
