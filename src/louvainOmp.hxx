#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "_main.hxx"
#include "properties.hxx"
#include "duplicate.hxx"
#include "modularity.hxx"
#include "louvain.hxx"

using std::tuple;
using std::vector;
using std::min;




// LOUVAIN-CHANGE-COMMUNITY
// ------------------------

/**
 * Move vertex to another community C.
 * @param vcom community each vertex belongs to (updated)
 * @param ctot total edge weight of each community (updated)
 * @param x original graph
 * @param u given vertex
 * @param c community to move to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
void louvainChangeCommunityOmp(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  #pragma omp atomic
  ctot[d] -= vtot[u];
  #pragma omp atomic
  ctot[c] += vtot[u];
  vcom[u] = c;
}




// LOUVAIN-MOVE
// ------------

/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param E tolerance
 * @param L max iterations
 * @param fa is a vertex affected?
 * @param fp process vertices whose communities have changed
 * @returns iterations performed
 */
template <class G, class K, class W, class FA, class FP>
int louvainMoveOmp(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L, FA fa, FP fp) {
  size_t S = x.span();
  int l = 0;
  for (; l<L;) {
    W el = W();
    #pragma omp parallel for schedule(auto) reduction(+:el)
    for (K u=K(); u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!fa(u)) continue;
      louvainClearScan(*vcs[t], *vcout[t]);
      louvainScanCommunities(*vcs[t], *vcout[t], x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (c)      { louvainChangeCommunityOmp(vcom, ctot, x, u, c, vtot); fp(u); }
      el += e;  // l1-norm
    } ++l;
    if (el<=E) break;
  }
  return l;
}
template <class G, class K, class W, class FA>
inline int louvainMoveOmp(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L, FA fa) {
  auto fp = [](auto u) {};
  return louvainMoveOmp(vcom, ctot, vcs, vcout, x, vtot, M, R, E, L, fa, fp);
}
template <class G, class K, class W>
inline int louvainMoveOmp(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L) {
  auto fa = [](auto u) { return true; };
  return louvainMoveOmp(vcom, ctot, vcs, vcout, x, vtot, M, R, E, L, fa);
}




// LOUVAIN-AGGREGATE
// -----------------

/**
 * Louvain algorithm's community aggregation phase.
 * @param a output graph
 * @param vcs communities vertex u is linked to (temporary buffer, updated)
 * @param vcout total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K, class W>
void louvainAggregateOmp(G& a, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  auto comv = louvainCommunityVertices(x, vcom);
  for (K c=0; c<comv.size(); ++c) {
    if (comv[c].empty()) continue;
    a.addVertex(c);
  }
  #pragma omp parallel for schedule(auto)
  for (K c=0; c<comv.size(); ++c) {
    int t = omp_get_thread_num();
    if (comv[c].empty()) continue;
    louvainClearScan(*vcs[t], *vcout[t]);
    for (K u : comv[c])
      louvainScanCommunities<true>(*vcs[t], *vcout[t], x, u, vcom);
    for (auto d : *vcs[t])
      a.addEdge(c, d, (*vcout[t])[d]);
  }
  a.correct();
}
template <class G, class K, class W>
inline auto louvainAggregateOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom) {
  G a; louvainAggregateOmp(a, vcs, vcout, x, vcom);
  return a;
}




// LOUVAIN-OMP
// -----------

template <class G, class K, class FA, class FP>
auto louvainOmp(const G& x, const vector<K>* q, const LouvainOptions& o, FA fa, FP fp) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  double R = o.resolution;
  double D = o.passTolerance;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  int    T = omp_get_max_threads();
  vector<K> vcom(S), a(S);
  vector<W> vtot(S), ctot(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  for (int t=0; t<T; ++t) {
    vcs[t]   = new vector<K>();
    vcout[t] = new vector<W>(S);
  }
  float t = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    double Q0 = modularity(x, M, R);
    G y = duplicate(x);
    fillValueU(vcom, K());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    mark([&]() {
      louvainVertexWeights(vtot, y);
      louvainInitialize(vcom, ctot, y, vtot);
      if (q) louvainInitializeFrom(vcom, ctot, x, vtot, *q);
      else   louvainInitialize(vcom, ctot, y, vtot);
      copyValues(vcom, a);
      for (l=0, p=0; M>0 && p<P;) {
        int m = 0;
        if (p==0) m = louvainMoveOmp(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L, fa, fp);
        else      m = louvainMoveOmp(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L);
        l += m; ++p;
        if (m<=1 || p>=P) { louvainLookupCommunities(a, vcom); break; }
        // K N0 = y.order();
        y = louvainAggregateOmp(vcs, vcout, y, vcom);
        // K N1 = y.order();
        // if (N1==N0) break;
        louvainLookupCommunities(a, vcom);
        PRINTFD("louvainOmp(): p=%d, l=%d, m=%d, Q=%f\n", p, l, m, modularity(y, M, R));
        double Q = D? modularity(y, M, R) : 0;
        if (D && Q-Q0<=D) break;
        fillValueU(vcom, K());
        fillValueU(vtot, W());
        fillValueU(ctot, W());
        louvainVertexWeights(vtot, y);
        louvainInitialize(vcom, ctot, y, vtot);
        E /= o.tolerenceDeclineFactor;
        Q0 = Q;
      }
    });
  }, o.repeat);
  for (int t=0; t<T; ++t) {
    delete vcs[t];
    delete vcout[t];
  }
  return LouvainResult<K>(a, l, p, t);
}
template <class G, class K, class FA>
inline auto louvainOmp(const G& x, const vector<K>* q, const LouvainOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return louvainOmp(x, q, o, fa, fp);
}
template <class G, class K>
inline auto louvainOmp(const G& x, const vector<K>* q, const LouvainOptions& o) {
  auto fa = [](auto u) { return true; };
  return louvainOmp(x, q, o, fa);
}




// LOUVAIN-OMP-STATIC
// ------------------

template <class G, class K>
inline auto louvainOmpStatic(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  return louvainOmp(x, q, o);
}




// LOUVAIN-OMP-DYNAMIC-DELTA-SCREENING
// -----------------------------------

template <class G, class K, class V>
inline auto louvainOmpDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  size_t S = x.span();
  double R = o.resolution;
  double M = edgeWeight(x)/2;
  const vector<K>& vcom = *q;
  vector<W> vtot(S), ctot(S);
  louvainVertexWeights(vtot, x);
  louvainCommunityWeights(ctot, x, vcom, vtot);
  auto vaff = louvainAffectedVerticesDeltaScreening(x, deletions, insertions, vcom, vtot, ctot, M, R);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return louvainOmp(x, q, o, fa);
}




// LOUVAIN-SEQ-DYNAMIC-FRONTIER
// ----------------------------

template <class G, class K, class V>
inline auto louvainOmpDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  size_t S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = louvainAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return louvainOmp(x, q, o, fa, fp);
}
