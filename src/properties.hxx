#pragma once
#include <vector>
#include <unordered_map>

using std::vector;
using std::unordered_map;




// WEIGHT
// ------

/**
 * Find the total outgoing edge weight of a vertex.
 * @param x original graph
 * @param u given vertex
 * @returns total outgoing weight of a vertex
 */
template <class G, class K>
inline double edgeWeight(const G& x, K u) {
  double a = 0;
  x.forEachEdgeValue(u, [&](auto w) { a += w; });
  return a;
}


/**
 * Find the total edge weight of a graph.
 * @param x original graph
 * @returns total edge weight (undirected graph => each edge considered twice)
 */
template <class G>
inline double edgeWeight(const G& x) {
  double a = 0;
  x.forEachVertexKey([&](auto u) { a += edgeWeight(x, u); });
  return a;
}

template <class G>
inline double edgeWeightOmp(const G& x) {
  using K = typename G::key_type;
  double a = 0;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a += edgeWeight(x, u);
  }
  return a;
}




// COUNT AS
// --------

/**
 * Count the number of vertices with each given value.
 * @param x original graph
 * @param vdata vertex values
 * @param fm the thing to be counted (value)
 */
template <class G, class V, class FM>
inline auto countAs(const G& x, const vector<V>& vdata, FM fm) {
  using T = decltype(fm(V()));
  unordered_map<T, size_t> a;
  x.forEachVertexKey([&](auto u) {
    auto d = fm(vdata[u]);
    a[d] += 1;
  });
  return a;
}
