#include <random>

using std::uniform_real_distribution;




// ADD RANDOM EDGE
// ---------------

template <class G, class R, class V, class FE>
inline bool addRandomEdge(const G& x, R& rnd, size_t span, V w, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(dis(rnd) * span);
  K v = K(dis(rnd) * span);
  fe(u, v, w);
  return true;
}

template <class G, class R, class V>
inline bool addRandomEdge(G& a, R& rnd, size_t span, V w) {
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(u, v, w); };
  return addRandomEdge(a, rnd, span, w, fe);
}


template <class G, class R, class V, class FE>
inline bool addRandomEdgeByDegree(const G& x, R& rnd, size_t span, V w, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  double deg = x.size() / x.span();
  size_t un = size_t(dis(rnd) * deg * span);
  size_t vn = size_t(dis(rnd) * deg * span);
  K u = 0, v = 0; size_t n = 0;
  x.forEachVertexKey([&](auto t) {
    if (un<0 && un > n+x.degree(t)) u = t;
    if (vn<0 && vn > n+x.degree(t)) v = t;
    if (un>0 && vn>=0) return;
    n += x.degree(t);
  });
  if (!u) u = K(un/deg);
  if (!v) v = K(vn/deg);
  fe(u, v, w);
  return true;
}

template <class G, class R, class V>
inline bool addRandomEdgeByDegree(G& a, R& rnd, size_t span, V w) {
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(u, v, w); };
  return addRandomEdgeByDegree(a, rnd, span, w, fe);
}




// REMOVE RANDOM EDGE
// ------------------

template <class G, class R, class K, class FE>
inline bool removeRandomEdgeFrom(const G& x, R& rnd, K u, FE fe) {
  uniform_real_distribution<> dis(0.0, 1.0);
  if (x.degree(u) == 0) return false;
  K vi = K(dis(rnd) * x.degree(u)), i = 0;
  bool removed = false;
  x.forEachEdgeKey(u, [&](auto v) {
    if (removed) return;
    if (i++ == vi) { fe(u, v); removed = true; }
  });
  return removed;
}

template <class G, class R, class K>
inline bool removeRandomEdgeFrom(G& a, R& rnd, K u) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); };
  return removeRandomEdgeFrom(a, rnd, u, fe);
}


template <class G, class R, class FE>
inline bool removeRandomEdge(const G& x, R& rnd, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  K u = K(dis(rnd) * x.span());
  return removeRandomEdgeFrom(x, rnd, u, fe);
}

template <class G, class R>
inline bool removeRandomEdge(G& a, R& rnd) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); };
  return removeRandomEdge(a, rnd, fe);
}


template <class G, class R, class FE>
inline bool removeRandomEdgeByDegree(const G& x, R& rnd, FE fe) {
  using K = typename G::key_type;
  uniform_real_distribution<> dis(0.0, 1.0);
  size_t v = size_t(dis(rnd) * x.size()), n = 0;
  bool attempted = false, removed = false;
  x.forEachVertexKey([&](auto u) {
    if (attempted) return;
    if (v > n+x.degree(u)) n += x.degree(u);
    else { removed = removeRandomEdge(x, rnd, u, fe); attempted = true; }
  });
  return removed;
}

template <class G, class R>
inline bool removeRandomEdgeByDegree(G& a, R& rnd) {
  auto fe = [&](auto u, auto v) { a.removeEdge(u, v); };
  return removeRandomEdgeByDegree(a, rnd, fe);
}
