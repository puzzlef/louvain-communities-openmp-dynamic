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

template <class G, class KS>
auto sourceOffsets(const G& x, const KS& ks) {
  vector<size_t> a; size_t i = 0;
  a.reserve(x.order()+1);
  for (auto u : ks) {
    a.push_back(i);
    i += x.degree(u);
  }
  a.push_back(i);
  return a;
}
template <class G>
inline auto sourceOffsets(const G& x) {
  return sourceOffsets(x, x.vertexKeys());
}




// DESTINATION-INDICES
// -------------------

template <class G, class KS>
auto destinationIndices(const G& x, const KS& ks) {
  using K = typename G::key_type;
  auto ids = valueIndicesUnorderedMap(ks); vector<K> a;
  for (auto u : ks)
    x.forEachEdgeKey(u, [&](auto v) { a.push_back(K(ids[v])); });
  return a;
}
template <class G>
auto destinationIndices(const G& x) {
  return destinationIndices(x, x.vertexKeys());
}




// CSR-EQUAL
// ---------

template <class K, class V>
int csrCompare(const size_t *xv, const K *xd, const K *xe, const V *xw, const size_t *yv, const K *yd, const K *ye, const V *yw, size_t N) {
  vector<size_t> ix, iy;
  for (size_t u=0; u<N; ++u) {
    size_t XOFF = xv[u];
    size_t YOFF = yv[u];
    size_t XDEG = xd? xd[u] : xv[u+1] - xv[u];
    size_t YDEG = yd? yd[u] : yv[u+1] - yv[u];
    if (XDEG!=YDEG) return sgn(XDEG - YDEG);
    if (XDEG<=0)  continue;
    ix.resize(XDEG);
    iy.resize(XDEG);
    iota(ix.begin(), ix.end(), XOFF);
    iota(iy.begin(), iy.end(), YOFF);
    sortValues(ix, [&](auto i, auto j) { return xe[i] < xe[j]; });
    sortValues(iy, [&](auto i, auto j) { return ye[i] < ye[j]; });
    for (size_t j=0; j<XDEG; ++j) {
      if (xe[ix[j]] != ye[iy[j]]) return sgn(xe[ix[j]] - ye[iy[j]]);
      if (!xw || !yw) continue;
      if (xw[ix[j]] != yw[iy[j]]) return sgn(xe[ix[j]] - ye[iy[j]]);
    }
  }
  return 0;
}
template <class K, class V>
inline int csrCompare(const vector<size_t>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw, const vector<size_t>& yv, const vector<K>& yd, const vector<K>& ye, const vector<V>& yw) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  const K *_yd = yd.empty()? nullptr : yd.data();
  const V *_xw = xw.empty()? nullptr : xw.data();
  const V *_yw = yw.empty()? nullptr : yw.data();
  if (xv.size() != yv.size()) return sgn(K(xv.size()) - K(yv.size()));
  return csrCompare(xv.data(), _xd, xe.data(), _xw, yv.data(), _yd, ye.data(), _yw, xv.size()-1);
}
template <class K>
inline int csrCompare(const vector<size_t>& xv, const vector<K>& xe, const vector<size_t>& yv, const vector<K>& ye) {
  vector<K> _;
  return csrCompare(xv, _, xe, _, yv, _, ye, _);
}

template <class K, class V>
inline bool csrEqual(const size_t *xv, const K *xd, const K *xe, const V *xw, const size_t *yv, const K *yd, const K *ye, const V *yw, size_t N) {
  return csrCompare(xv, xd, xe, xw, yv, yd, ye, yw, N)==0;
}
template <class K, class V>
inline bool csrEqual(const vector<size_t>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw, const vector<size_t>& yv, const vector<K>& yd, const vector<K>& ye, const vector<V>& yw) {
  return csrCompare(xv, xd, xe, xw, yv, yd, ye, yw)==0;
}
template <class K>
inline bool csrEqual(const vector<size_t>& xv, const vector<K>& xe, const vector<size_t>& yv, const vector<K>& ye) {
  return csrCompare(xv, xe, yv, ye)==0;
}




// CSR-GRAPH
// ---------

template <class G, class K, class V>
void csrGraphW(G& a, const size_t *xv, const K *xd, const K *xe, const V *xw, size_t N) {
  for (size_t u=0; u<N; ++u)
    a.addVertex(u);
  for (size_t u=0; u<N; ++u) {
    size_t OFF = xv[u];
    size_t DEG = xd? xd[u] : xv[u+1] - xv[u];
    for (size_t j=0; j<DEG; ++j) {
      K v = xe[OFF+j];
      V w = xw? xw[OFF+j] : V(1);
      a.addEdge(u, v, w);
    }
  }
  a.correct();
}
template <class G, class K, class V>
void csrGraphW(G& a, const vector<size_t>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  const V *_xw = xw.empty()? nullptr : xw.data();
  csrGraphW(a, xv.data(), _xd, xe.data(), _xw, xv.size()-1);
}
template <class G, class K>
void csrGraphW(G& a, const vector<size_t>& xv, const vector<K>& xe) {
  vector<K> _;
  csrGraphW(a, xv.data(), _, xe.data(), _);
}

template <class G, class K, class V>
auto csrGraph(const size_t *xv, const K *xd, const K *xe, const V *xw, size_t N) {
  OutDiGraph<K, None, V> a; csrGraphW(a, xv, xd, xe, xw, N);
  return a;
}
template <class K, class V>
auto csrGraph(const vector<size_t>& xv, const vector<K>& xd, const vector<K>& xe, const vector<V>& xw) {
  OutDiGraph<K, None, V> a; csrGraphW(a, xv, xd, xe, xw);
  return a;
}
template <class K>
auto csrGraph(const vector<size_t>& xv, const vector<K>& xe) {
  OutDiGraph<K> a; csrGraphW(a, xv, xe);
  return a;
}




// CSR-SUM-EDGE-VALUES
// -------------------

template <class K, class V, class A=double>
A csrSumEdgeValues(const size_t *xv, const K *xd, const V *xw, size_t N, A a=A()) {
  for (size_t u=0; u<N; ++u) {
    size_t OFF = xv[u];
    size_t DEG = xd? xd[u] : xv[u+1] - xv[u];
    for (size_t j=0; j<DEG; ++j)
      a += xw[OFF+j];
  }
  return a;
}
template <class K, class V, class A=double>
inline V csrSumEdgeValues(const vector<size_t>& xv, const vector<K>& xd, const vector<V>& xw, A a=A()) {
  const K *_xd = xd.empty()? nullptr : xd.data();
  return csrSumEdgeValues(xv.data(), _xd, xw.data(), xv.size()-1, a);
}
