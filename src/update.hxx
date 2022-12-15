#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <omp.h>

using std::pair;
using std::vector;




// ADD EDGE
// --------
// Add an edge (in parallel).

template <class G, class K, class E>
inline void addEdgeU(G& a, K u, K v, E w=E()) {
  a.addEdge(u, v, w);
}

template <class G, class K, class E>
inline void addEdgeOmpU(G& a, K u, K v, E w=E()) {
  if (belongsOmp(u)) a.addEdge(u, v, w);
}




// UPDATE
// ------
// Update changes made to a graph.

template <class G>
inline void updateU(G& a) {
  a.update();
}

template <class G>
inline void updateOmpU(G& a) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t S = a.span();
  // Create per-thread buffers for update operation.
  int THREADS = omp_get_max_threads();
  vector<vector<pair<K, E>>*> bufs(THREADS);
  for (int i=0; i<THREADS; ++i)
    bufs[i] = new vector<pair<K, E>>();
  // Update edges of each vertex individually.
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    a.updateEdges(u, bufs[t]);
  }
  // Update the entire graph, find total vertices and edges.
  a.update();
  // Clean up.
  for (int i=0; i<THREADS; ++i)
    delete bufs[i];
}
