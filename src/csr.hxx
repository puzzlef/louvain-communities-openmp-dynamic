#pragma once
#include <numeric>
#include <algorithm>
#include <vector>
#include "_main.hxx"

using std::vector;
using std::iota;
using std::equal;
using std::transform;




// SOURCE-OFFSETS
// --------------

template <class G, class J, class T>
auto sourceOffsetsAs(const G& x, const J& ks, T _) {
  vector<T> a; T i = 0;
  a.reserve(x.order()+1);
  for (auto u : ks) {
    a.push_back(i);
    i += x.degree(u);
  }
  a.push_back(i);
  return a;
}
template <class G, class T>
inline auto sourceOffsetsAs(const G& x, T _) {
  return sourceOffsetsAs(x, x.vertexKeys(), _);
}

template <class G, class J>
inline auto sourceOffsets(const G& x, const J& ks) {
  return sourceOffsetsAs(x, ks, size_t());
}
template <class G>
inline auto sourceOffsets(const G& x) {
  return sourceOffsetsAs(x, size_t());
}




// DESTINATION-INDICES
// -------------------

template <class G, class J, class T>
auto destinationIndicesAs(const G& x, const J& ks, T _) {
  auto ids = valueIndicesUnorderedMap(ks); vector<T> a;
  for (auto u : ks)
    x.forEachEdgeKey(u, [&](auto v) { a.push_back(T(ids[v])); });
  return a;
}
template <class G, class T>
auto destinationIndicesAs(const G& x, T _) {
  return destinationIndicesAs(x, x.vertexKeys(), _);
}

template <class G, class J>
auto destinationIndices(const G& x, const J& ks) {
  using K = typename G::key_type;
  return destinationIndicesAs(x, ks, K());
}
template <class G>
inline auto destinationIndices(const G& x) {
  using K = typename G::key_type;
  return destinationIndicesAs(x, K());
}




// CSR-EQUAL
// ---------

template <class K, class V>
int csrCompare(const K *xv, const K *xd, const K *xe, const V *xw, const K *yv, const K *yd, const K *ye, const V *yw, K N) {
  vector<K> ix, iy;
  for (K u=0; u<N; ++u) {
    K XOFF = xv[u];
    K YOFF = yv[u];
    K XDEG = xd? xd[u] : xv[u+1] - xv[u];
    K YDEG = yd? yd[u] : yv[u+1] - yv[u];
    if (XDEG!=YDEG) return sgn(XDEG - YDEG);
    if (XDEG<=K())  continue;
    ix.resize(XDEG);
    iy.resize(XDEG);
    iota(ix.begin(), ix.end(), XOFF);
    iota(iy.begin(), iy.end(), YOFF);
    sortValues(ix, [&](auto i, auto j) { return xe[i] < xe[j]; });
    sortValues(iy, [&](auto i, auto j) { return ye[i] < ye[j]; });
    for (K j=K(); j<XDEG; ++j) {
      if (xe[ix[j]] != ye[iy[j]]) return sgn(xe[ix[j]] - ye[iy[j]]);
      if (!xw || !yw) continue;
      if (xw[ix[j]] != yw[iy[j]]) return sgn(xe[ix[j]] - ye[iy[j]]);
    }
  }
  return 0;
}
template <class K, class V>
inline int csrCompare(const vector<K>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw, const vector<K>& yv, const vector<K>& yd, const vector<K>& ye, const vector<V>& yw) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  const K *_yd = yd.empty()? nullptr : yd.data();
  const V *_xw = xw.empty()? nullptr : xw.data();
  const V *_yw = yw.empty()? nullptr : yw.data();
  if (xv.size() != yv.size()) return sgn(K(xv.size()) - K(yv.size()));
  return csrCompare(xv.data(), _xd, xe.data(), _xw, yv.data(), _yd, ye.data(), _yw, K(xv.size()-1));
}
template <class K>
inline int csrCompare(const vector<K>& xv, const vector<K>& xe, const vector<K>& yv, const vector<K>& ye) {
  vector<K> _;
  return csrCompare(xv, _, xe, _, yv, _, ye, _);
}

template <class K, class V>
inline bool csrEqual(const K *xv, const K *xd, const K *xe, const V *xw, const K *yv, const K *yd, const K *ye, const V *yw, size_t N) {
  return csrCompare(xv, xd, xe, xw, yv, yd, ye, yw, N)==0;
}
template <class K, class V>
inline bool csrEqual(const vector<K>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw, const vector<K>& yv, const vector<K>& yd, const vector<K>& ye, const vector<V>& yw) {
  return csrCompare(xv, xd, xe, xw, yv, yd, ye, yw)==0;
}
template <class K>
inline bool csrEqual(const vector<K>& xv, const vector<K>& xe, const vector<K>& yv, const vector<K>& ye) {
  return csrCompare(xv, xe, yv, ye)==0;
}




// CSR-GRAPH
// ---------

template <class G, class K, class V>
void csrGraphW(G& a, const K *xv, const K *xd, const K *xe, const V *xw, size_t N) {
  for (size_t u=0; u<N; ++u)
    a.addVertex(u);
  for (size_t u=0; u<N; ++u) {
    K OFF = xv[u];
    K DEG = xd? xd[u] : xv[u+1] - xv[u];
    for (K j=K(); j<DEG; ++j) {
      K v = xe[OFF+j];
      V w = xw? xw[OFF+j] : V(1);
      a.addEdge(u, v, w);
    }
  }
  a.correct();
}
template <class G, class K, class V>
void csrGraphW(G& a, const vector<K>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  const V *_xw = xw.empty()? nullptr : xw.data();
  csrGraphW(a, xv.data(), _xd, xe.data(), _xw, xv.size()-1);
}
template <class G, class K>
void csrGraphW(G& a, const vector<K>& xv, const vector<K>& xe) {
  vector<K> _;
  csrGraphW(a, xv.data(), _, xe.data(), _);
}

template <class G, class K, class V>
auto csrGraph(const K *xv, const K *xd, const K *xe, const V *xw, size_t N) {
  OutDiGraph<K, None, V> a; csrGraphW(a, xv, xd, xe, xw, N);
  return a;
}
template <class K, class V>
auto csrGraph(const vector<K>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw) {
  OutDiGraph<K, None, V> a; csrGraphW(a, xv, xd, xe, xw);
  return a;
}
template <class K>
auto csrGraph(const vector<K>& xv, const vector<K>& xe) {
  OutDiGraph<K> a; csrGraphW(a, xv, xe);
  return a;
}




// CSR-SUM-EDGE-VALUES
// -------------------

template <class K, class V>
V csrSumEdgeValues(const K *xv, const K *xd, const V *xw, K N, V a=V()) {
  for (K u=0; u<N; ++u) {
    K OFF = xv[u];
    K DEG = xd? xd[u] : xv[u+1] - xv[u];
    for (K j=K(); j<DEG; ++j)
      a += xw[OFF+j];
  }
  return a;
}
template <class K, class V>
inline V csrSumEdgeValues(const vector<K>& xv, const vector<K>& xd, const vector<V>& xw) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  return csrSumEdgeValues(xv.data(), _xd, xw.data(), K(xv.size()-1));
}
