#pragma once
#include <tuple>
#include "vertices.hxx"

using std::make_tuple;




// DEGREES (ALL, MIN, MAX, AVG)
// ----------------------------

template <class G, class F>
auto degreesDo(const G& x, F fn) {
  using K = typename G::key_type;
  auto  a = createContainer(x, K());
  x.forEachVertexKey([&](auto u) { a[u] = x.degree(u); fn(u, x.degree(u)); });
  return a;
}
template <class G>
inline auto degrees(const G& x) {
  auto fn = [](auto u, auto d) {};
  return degreesDo(x, fn);
}


template <class G>
auto minDegree(const G& x) {
  using K = typename G::key_type;
  K min = K(x.order());
  x.forEachVertexKey([&](auto u) {
    auto d = x.degree(u);
    if (d<min) min = d;
  });
  return min;
}

template <class G>
auto maxDegree(const G& x) {
  using K = typename G::key_type;
  K max = 0;
  x.forEachVertexKey([&](auto u) {
    auto d = x.degree(u);
    if (d>max) max = d;
  });
  return max;
}

template <class G>
inline double avgDegree(const G& x) {
  size_t N = x.order();
  return N==0? 0 : x.size()/double(N);
}


template <class G>
auto minMaxAvgDegree(const G& x) {
  using K = typename G::key_type;
  K min = K(x.order()), max = 0;
  x.forEachVertexKey([&](auto u) {
    auto d = x.degree(u);
    if (d<min) min = d;
    if (d>max) max = d;
  });
  return make_tuple(min, max, avgDegree(x));
}




// DENSITY
// -------
// Fully connectedness fraction.

template <class G>
inline double density(const G& x) {
  double N = x.order();
  return N>0? x.size()/(N*N) : 0;
}




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
