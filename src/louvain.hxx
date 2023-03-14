#pragma once
#include <utility>
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "Graph.hxx"
#include "duplicate.hxx"
#include "properties.hxx"
#include "modularity.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;
using std::move;
using std::get;
using std::min;




// LOUVAIN OPTIONS
// ---------------

struct LouvainOptions {
  int    repeat;
  double resolution;
  double tolerance;
  double passTolerance;
  double tolerenceDeclineFactor;
  int    maxIterations;
  int    maxPasses;

  LouvainOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double passTolerance=0, double tolerenceDeclineFactor=10, int maxIterations=20, int maxPasses=20) :
  repeat(repeat), resolution(resolution), tolerance(tolerance), passTolerance(passTolerance), tolerenceDeclineFactor(tolerenceDeclineFactor), maxIterations(maxIterations), maxPasses(maxPasses) {}
};

// Weight to be using in hashtable.
#define LOUVAIN_WEIGHT_TYPE double




// LOUVAIN RESULT
// --------------

template <class K>
struct LouvainResult {
  vector<K> membership;
  int   iterations;
  int   passes;
  float time;
  float preprocessingTime;

  LouvainResult(vector<K>&& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0) :
  membership(membership), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime) {}

  LouvainResult(vector<K>& membership, int iterations=0, int passes=0, float time=0, float preprocessingTime=0) :
  membership(move(membership)), iterations(iterations), passes(passes), time(time), preprocessingTime(preprocessingTime) {}
};




// LOUVAIN INITIALIZE
// ------------------

/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, should be initialized to 0)
 * @param x original graph
 */
template <class G, class W>
inline void louvainVertexWeights(vector<W>& vtot, const G& x) {
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      vtot[u] += w;
    });
  });
}

#ifdef OPENMP
template <class G, class W>
inline void louvainVertexWeightsOmp(vector<W>& vtot, const G& x) {
  using  K = typename G::key_type;
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    x.forEachEdge(u, [&](auto v, auto w) { vtot[u] += w; });
  }
}
#endif


/**
 * Find the total edge weight of each community.
 * @param ctot total edge weight of each community (updated, should be initialized to 0)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void louvainCommunityWeights(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ctot[c] += vtot[u];
  });
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainCommunityWeightsOmp(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u];
    #pragma omp atomic
    ctot[c] += vtot[u];
  }
}
#endif


/**
 * Initialize communities such that each vertex is its own community.
 * @param vcom community each vertex belongs to (updated, should be initialized to 0)
 * @param ctot total edge weight of each community (updated, should be initilized to 0)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 */
template <class G, class K, class W>
inline void louvainInitialize(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    vcom[u] = u;
    ctot[u] = vtot[u];
  });
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainInitializeOmp(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  size_t S = x.span();
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    vcom[u] = u;
    ctot[u] = vtot[u];
  }
}
#endif


/**
 * Initialize communities from given initial communities.
 * @param vcom community each vertex belongs to (updated, should be initialized to 0)
 * @param ctot total edge weight of each community (updated, should be initilized to 0)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param q initial community each vertex belongs to
 */
template <class G, class K, class W>
inline void louvainInitializeFrom(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  copyValuesW(vcom, q, 0, min(q.size(), vcom.size()));
  louvainCommunityWeights(ctot, x, vcom, vtot);
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainInitializeFromOmp(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  copyValuesOmpW(vcom, q, 0, min(q.size(), vcom.size()));
  louvainCommunityWeightsOmp(ctot, x, vcom, vtot);
}
#endif




// LOUVAIN COMMUNITY VERTICES
// --------------------------

template <class G, class K>
inline auto louvainCommunityVertices(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  x.forEachVertexKey([&](auto u) { a[vcom[u]].push_back(u); });
  return a;
}

#ifdef OPENMP
template <class G, class K>
inline auto louvainCommunityVerticesOmp(const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  vector2d<K> a(S);
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      if (belongsOmp(vcom[u])) a[vcom[u]].push_back(u);
    });
  }
  return a;
}
#endif




// LOUVAIN LOOKUP COMMUNITIES
// --------------------------

/**
 * Update community membership in a tree-like fashion (to handle aggregation).
 * @param a output community each vertex belongs to (updated)
 * @param vcom community each vertex belongs to (at this aggregation level)
 */
template <class K>
inline void louvainLookupCommunities(vector<K>& a, const vector<K>& vcom) {
  for (auto& v : a)
    v = vcom[v];
}

#ifdef OPENMP
template <class K>
inline void louvainLookupCommunitiesOmp(vector<K>& a, const vector<K>& vcom) {
  size_t S = a.size();
  #pragma omp parallel for schedule(auto)
  for (size_t u=0; u<S; ++u)
    a[u] = vcom[a[u]];
}
#endif




// LOUVAIN CHANGE COMMUNITY
// ------------------------

/**
 * Scan an edge community connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, class W>
inline void louvainScanCommunity(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
  if (!SELF && u==v) return;
  K c = vcom[v];
  if (!vcout[c]) vcs.push_back(c);
  vcout[c] += w;
}


/**
 * Scan communities connected to a vertex.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class G, class K, class W>
inline void louvainScanCommunities(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { louvainScanCommunity<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void louvainClearScan(vector<K>& vcs, vector<W>& vcout) {
  for (K c : vcs)
    vcout[c] = W();
  vcs.clear();
}


/**
 * Choose connected community with best delta modularity.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param vcs communities vertex u is linked to
 * @param vcout total edge weight from vertex u to community C
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns [best community, delta modularity]
 */
template <bool SELF=false, class G, class K, class W>
inline auto louvainChooseCommunity(const G& x, K u, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, const vector<K>& vcs, const vector<W>& vcout, double M, double R) {
  K cmax = K(), d = vcom[u];
  W emax = W();
  for (K c : vcs) {
    if (!SELF && c==d) continue;
    W e = deltaModularity(vcout[c], vcout[d], vtot[u], ctot[c], ctot[d], M, R);
    if (e>emax) { emax = e; cmax = c; }
  }
  return make_pair(cmax, emax);
}


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
inline void louvainChangeCommunity(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  ctot[d] -= vtot[u];
  ctot[c] += vtot[u];
  vcom[u] = c;
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainChangeCommunityOmp(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  #pragma omp atomic
  ctot[d] -= vtot[u];
  #pragma omp atomic
  ctot[c] += vtot[u];
  vcom[u] = c;
}
#endif




// LOUVAIN MOVE
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
inline int louvainMove(vector<K>& vcom, vector<W>& ctot, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L, FA fa, FP fp) {
  int l = 0;
  for (; l<L;) {
    W el = W();
    x.forEachVertexKey([&](auto u) {
      if (!fa(u)) return;
      louvainClearScan(vcs, vcout);
      louvainScanCommunities(vcs, vcout, x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
      if (c)      { louvainChangeCommunity(vcom, ctot, x, u, c, vtot); fp(u); }
      el += e;  // l1-norm
    }); ++l;
    if (el<=E) break;
  }
  return l;
}

#ifdef OPENMP
template <class G, class K, class W, class FA, class FP>
inline int louvainMoveOmp(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L, FA fa, FP fp) {
  size_t S = x.span();
  int l = 0;
  for (; l<L;) {
    W el = W();
    #pragma omp parallel for schedule(auto) reduction(+:el)
    for (K u=0; u<S; ++u) {
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
#endif


template <class G, class K, class W, class FA>
inline int louvainMove(vector<K>& vcom, vector<W>& ctot, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L, FA fa) {
  auto fp = [](auto u) {};
  return louvainMove(vcom, ctot, vcs, vcout, x, vtot, M, R, E, L, fa, fp);
}
template <class G, class K, class W>
inline int louvainMove(vector<K>& vcom, vector<W>& ctot, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, double E, int L) {
  auto fa = [](auto u) { return true; };
  return louvainMove(vcom, ctot, vcs, vcout, x, vtot, M, R, E, L, fa);
}

#ifdef OPENMP
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
#endif




// LOUVAIN AGGREGATE
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
inline void louvainAggregate(G& a, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom) {
  auto comv = louvainCommunityVertices(x, vcom);
  for (K c=0; c<comv.size(); ++c) {
    if (comv[c].empty()) continue;
    louvainClearScan(vcs, vcout);
    for (K u : comv[c])
      louvainScanCommunities<true>(vcs, vcout, x, u, vcom);
    a.addVertex(c);
    for (auto d : vcs)
      a.addEdge(c, d, vcout[d]);
  }
  a.update();
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainAggregateOmp(G& a, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom) {
  size_t  S = x.span();
  auto comv = louvainCommunityVerticesOmp(x, vcom);
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
  updateOmpU(a);
}
#endif


template <class G, class K, class W>
inline auto louvainAggregate(vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom) {
  G a; louvainAggregate(a, vcs, vcout, x, vcom);
  return a;
}

#ifdef OPENMP
template <class G, class K, class W>
inline auto louvainAggregateOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom) {
  G a; louvainAggregateOmp(a, vcs, vcout, x, vcom);
  return a;
}
#endif




// LOUVAIN AFFECTED VERTICES DELTA-SCREENING
// -----------------------------------------
// Using delta-screening approach.
// - All edge batches are undirected, and sorted by source vertex-id.
// - For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
//   `i`'s neighbors and `j*`'s community is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i`'s neighbors and `j`'s community is marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns flags for each vertex marking whether it is affected
 */
template <class G, class K, class V, class W>
inline auto louvainAffectedVerticesDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R=1) {
  size_t S = x.span();
  vector<K> vcs; vector<W> vcout(S);
  vector<bool> vertices(S), neighbors(S), communities(S);
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[vcom[v]] = true;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    louvainClearScan(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v = get<1>(insertions[i]);
      V w = get<2>(insertions[i]);
      if (vcom[u] == vcom[v]) continue;
      louvainScanCommunity(vcs, vcout, u, v, w, vcom);
    }
    auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
    if (e<=0) continue;
    vertices[u]  = true;
    neighbors[u] = true;
    communities[c] = true;
  }
  x.forEachVertexKey([&](auto u) {
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = true; });
    if (communities[vcom[u]]) vertices[u] = true;
  });
  return vertices;
}




// LOUVAIN AFFECTED VERTICES FRONTIER
// ----------------------------------
// Using frontier based approach.
// - All source and destination vertices are marked as affected for insertions and deletions.
// - For edge additions across communities with source vertex `i` and destination vertex `j`,
//   `i` is marked as affected.
// - For edge deletions within the same community `i` and `j`,
//   `i` is marked as affected.
// - Vertices whose communities change in local-moving phase have their neighbors marked as affected.

/**
 * Find the vertices which should be processed upon a batch of edge insertions and deletions.
 * @param x original graph
 * @param deletions edge deletions for this batch update (undirected, sorted by source vertex id)
 * @param insertions edge insertions for this batch update (undirected, sorted by source vertex id)
 * @param vcom community each vertex belongs to
 * @returns flags for each vertex marking whether it is affected
 */
template <class G, class K, class V>
inline auto louvainAffectedVerticesFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  size_t S = x.span();
  vector<bool> vertices(S);
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = true;
  }
  for (const auto& [u, v, w] : insertions) {
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = true;
  }
  return vertices;
}




// LOUVAIN
// -------

template <class G, class K, class FA, class FP>
auto louvainSeq(const G& x, const vector<K>* q, const LouvainOptions& o, FA fa, FP fp) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  double R = o.resolution;
  double D = o.passTolerance;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  vector<K> vcom(S), vcs, a(S);
  vector<W> vtot(S), ctot(S), vcout(S);
  float t = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    double Q0 = modularity(x, M, R);
    G y = duplicate(x);
    fillValueU(vcom, K());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    mark([&]() {
      louvainVertexWeights(vtot, y);
      if (q) louvainInitializeFrom(vcom, ctot, x, vtot, *q);
      else   louvainInitialize(vcom, ctot, y, vtot);
      copyValues(vcom, a);
      for (l=0, p=0; M>0 && p<P;) {
        int m = 0;
        if (p==0) m = louvainMove(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L, fa, fp);
        else      m = louvainMove(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L);
        l += m; ++p;
        if (m<=1 || p>=P) { louvainLookupCommunities(a, vcom); break; }
        // K N0 = y.order();
        y = louvainAggregate(vcs, vcout, y, vcom);
        // K N1 = y.order();
        // if (N1==N0) break;
        louvainLookupCommunities(a, vcom);
        PRINTFD("louvainSeq(): p=%d, l=%d, m=%d, Q=%f\n", p, l, m, modularity(y, M, R));
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
  return LouvainResult<K>(a, l, p, t);
}

#ifdef OPENMP
template <class G, class K, class FA, class FP>
auto louvainOmp(const G& x, const vector<K>* q, const LouvainOptions& o, FA fa, FP fp) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  double R = o.resolution;
  double D = o.passTolerance;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
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
    double Q0 = modularityOmp(x, M, R);
    G y = duplicate(x);
    fillValueOmpU(vcom, K());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    mark([&]() {
      louvainVertexWeightsOmp(vtot, x);
      if (q) louvainInitializeFromOmp(vcom, ctot, x, vtot, *q);
      else   louvainInitializeOmp(vcom, ctot, x, vtot);
      copyValuesOmpW(a, vcom);
      for (l=0, p=0; M>0 && p<P;) {
        int m = 0;
        if (p==0) m = louvainMoveOmp(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L, fa, fp);
        else      m = louvainMoveOmp(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L);
        l += m; ++p;
        if (m<=1 || p>=P) { louvainLookupCommunitiesOmp(a, vcom); break; }
        // K N0 = y.order();
        y = louvainAggregateOmp(vcs, vcout, y, vcom);
        // K N1 = y.order();
        // if (N1==N0) break;
        louvainLookupCommunitiesOmp(a, vcom);
        LOGD("louvainOmp(): p=%d, l=%d, m=%d, Q=%f\n", p, l, m, modularityOmp(y, M, R));
        double Q = D? modularityOmp(y, M, R) : 0;
        if (D && Q-Q0<=D) break;
        fillValueOmpU(vcom, K());
        fillValueOmpU(vtot, W());
        fillValueOmpU(ctot, W());
        louvainVertexWeightsOmp(vtot, y);
        louvainInitializeOmp(vcom, ctot, y, vtot);
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
#endif


template <class G, class K, class FA>
inline auto louvainSeq(const G& x, const vector<K>* q, const LouvainOptions& o, FA fa) {
  auto fp = [](auto u) {};
  return louvainSeq(x, q, o, fa, fp);
}
template <class G, class K>
inline auto louvainSeq(const G& x, const vector<K>* q, const LouvainOptions& o) {
  auto fa = [](auto u) { return true; };
  return louvainSeq(x, q, o, fa);
}

#ifdef OPENMP
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
#endif




// LOUVAIN-STATIC
// --------------

template <class G, class K>
inline auto louvainStaticSeq(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  return louvainSeq(x, q, o);
}

#ifdef OPENMP
template <class G, class K>
inline auto louvainStaticOmp(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  return louvainOmp(x, q, o);
}
#endif




// LOUVAIN DYNAMIC DELTA-SCREENING
// -------------------------------

template <class G, class K, class V>
inline auto louvainDynamicDeltaScreeningSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  size_t S = x.span();
  double R = o.resolution;
  double M = edgeWeight(x)/2;
  const vector<K>& vcom = *q;
  vector<V> vtot(S), ctot(S);
  louvainVertexWeights(vtot, x);
  louvainCommunityWeights(ctot, x, vcom, vtot);
  auto vaff = louvainAffectedVerticesDeltaScreening(x, deletions, insertions, vcom, vtot, ctot, M, R);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return louvainSeq(x, q, o, fa);
}

#ifdef OPENMP
template <class G, class K, class V>
inline auto louvainDynamicDeltaScreeningOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  size_t S = x.span();
  double R = o.resolution;
  double M = edgeWeightOmp(x)/2;
  const vector<K>& vcom = *q;
  vector<W> vtot(S), ctot(S);
  louvainVertexWeightsOmp(vtot, x);
  louvainCommunityWeightsOmp(ctot, x, vcom, vtot);
  auto vaff = louvainAffectedVerticesDeltaScreening(x, deletions, insertions, vcom, vtot, ctot, M, R);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return louvainOmp(x, q, o, fa);
}
#endif




// LOUVAIN DYNAMIC FRONTIER
// ------------------------

template <class G, class K, class V>
inline auto louvainDynamicFrontierSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  size_t S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = louvainAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return louvainSeq(x, q, o, fa, fp);
}

#ifdef OPENMP
template <class G, class K, class V>
inline auto louvainDynamicFrontierOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  const vector<K>& vcom = *q;
  auto vaff = louvainAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return louvainOmp(x, q, o, fa, fp);
}
#endif
