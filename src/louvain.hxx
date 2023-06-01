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
  double aggregationTolerance;
  double toleranceDecline;
  int    maxIterations;
  int    maxPasses;

  LouvainOptions(int repeat=1, double resolution=1, double tolerance=1e-2, double aggregationTolerance=0.8, double toleranceDecline=100, int maxIterations=20, int maxPasses=10) :
  repeat(repeat), resolution(resolution), tolerance(tolerance), aggregationTolerance(aggregationTolerance), toleranceDecline(toleranceDecline), maxIterations(maxIterations), maxPasses(maxPasses) {}
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




// LOUVAIN HASHTABLES
// ------------------

/**
 * Allocate a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 * @param S size of each hashtable
 */
template <class K, class W>
inline void louvainAllocateHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, size_t S) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    vcs[i]   = new vector<K>();
    vcout[i] = new vector<W>(S);
  }
}


/**
 * Free a number of hashtables.
 * @param vcs communities vertex u is linked to (updated)
 * @param vcout total edge weight from vertex u to community C (updated)
 */
template <class K, class W>
inline void louvainFreeHashtablesW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout) {
  size_t N = vcs.size();
  for (size_t i=0; i<N; ++i) {
    delete vcs[i];
    delete vcout[i];
  }
}




// LOUVAIN INITIALIZE
// ------------------

/**
 * Find the total edge weight of each vertex.
 * @param vtot total edge weight of each vertex (updated, should be initialized to 0)
 * @param x original graph
 */
template <class G, class W>
inline void louvainVertexWeightsW(vector<W>& vtot, const G& x) {
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      vtot[u] += w;
    });
  });
}

#ifdef OPENMP
template <class G, class W>
inline void louvainVertexWeightsOmpW(vector<W>& vtot, const G& x) {
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
inline void louvainCommunityWeightsW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    ctot[c] += vtot[u];
  });
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainCommunityWeightsOmpW(vector<W>& ctot, const G& x, const vector<K>& vcom, const vector<W>& vtot) {
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
inline void louvainInitializeW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
  x.forEachVertexKey([&](auto u) {
    vcom[u] = u;
    ctot[u] = vtot[u];
  });
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainInitializeOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot) {
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
inline void louvainInitializeFromW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  copyValuesW(vcom.data(), q.data(), min(q.size(), vcom.size()));
  louvainCommunityWeightsW(ctot, x, vcom, vtot);
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainInitializeFromOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, const vector<W>& vtot, const vector<K>& q) {
  copyValuesOmpW(vcom.data(), q.data(), min(q.size(), vcom.size()));
  louvainCommunityWeightsOmpW(ctot, x, vcom, vtot);
}
#endif




// LOUVAIN COMMUNITY VERTICES
// --------------------------

/**
 * Find the number of vertices in each community.
 * @param a number of vertices belonging to each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @returns number of communities
 */
template <class G, class K>
inline size_t louvainCountCommunityVerticesW(K *a, const G& x, const K *vcom) {
  size_t S = x.span();
  size_t n = 0;
  fillValueU(a, S, K());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    if (a[c]==0) ++n;
    ++a[c];
  });
  return n;
}
template <class G, class K>
inline size_t louvainCountCommunityVerticesW(vector<K>& a, const G& x, const vector<K>& vcom) {
  return louvainCountCommunityVerticesW(a.data(), x, vcom.data());
}


#ifdef OPENMP
template <class G, class K>
inline size_t louvainCountCommunityVerticesOmpW(K *a, const G& x, const K *vcom) {
  size_t S = x.span();
  size_t n = 0;
  fillValueOmpU(a, S, K());
  #pragma omp parallel for schedule(auto) reduction(+:n)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u], m = 0;
    #pragma omp atomic capture
    { m = a[c]; ++a[c]; }
    if (m==0) ++n;
  }
  return n;
}
template <class G, class K>
inline size_t louvainCountCommunityVerticesOmpW(vector<K>& a, const G& x, const vector<K>& vcom) {
  return louvainCountCommunityVerticesOmpW(a.data(), x, vcom.data());
}
#endif




/**
 * Find the vertices in each community.
 * @param co csr offsets for vertices belonging to each community (updated)
 * @param ce csr data vertices belonging to each community (updated)
 * @param cn number of vertices in each community (updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 */
template <class G, class K>
inline void louvainCommunityVerticesW(vector<K>& co, vector<K>& ce, vector<K>& cn, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  co[S] = exclusiveScanW(co, cn);
  fillValueU(cn, K());
  x.forEachVertexKey([&](auto u) {
    K c = vcom[u];
    K i = cn[c]++;
    ce[co[c] + i] = u;
  });
}


#ifdef OPENMP
template <class G, class K>
inline void louvainCommunityVerticesOmpW(vector<K>& co, vector<K>& ce, vector<K>& cn, vector<K>& bufk, const G& x, const vector<K>& vcom) {
  size_t S = x.span();
  co[S] = exclusiveScanOmpW(co, bufk, cn);
  fillValueOmpU(cn, K());
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    K c = vcom[u], i = 0;
    #pragma omp atomic capture
    { i = cn[c]; ++cn[c]; }
    ce[co[c] + i] = u;
  }
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
inline void louvainLookupCommunitiesU(vector<K>& a, const vector<K>& vcom) {
  for (auto& v : a)
    v = vcom[v];
}

#ifdef OPENMP
template <class K>
inline void louvainLookupCommunitiesOmpU(vector<K>& a, const vector<K>& vcom) {
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
inline void louvainScanCommunityW(vector<K>& vcs, vector<W>& vcout, K u, K v, V w, const vector<K>& vcom) {
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
inline void louvainScanCommunitiesW(vector<K>& vcs, vector<W>& vcout, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { louvainScanCommunityW<SELF>(vcs, vcout, u, v, w, vcom); });
}


/**
 * Clear communities scan data.
 * @param vcs total edge weight from vertex u to community C (updated)
 * @param vcout communities vertex u is linked to (updated)
 */
template <class K, class W>
inline void louvainClearScanW(vector<K>& vcs, vector<W>& vcout) {
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
inline void louvainChangeCommunityW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
  K d = vcom[u];
  ctot[d] -= vtot[u];
  ctot[c] += vtot[u];
  vcom[u] = c;
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainChangeCommunityOmpW(vector<K>& vcom, vector<W>& ctot, const G& x, K u, K c, const vector<W>& vtot) {
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
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @param fa is a vertex affected?
 * @param fp process vertices whose communities have changed
 * @returns iterations performed
 */
template <class G, class K, class W, class FC, class FA, class FP>
inline int louvainMoveW(vector<K>& vcom, vector<W>& ctot, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc, FA fa, FP fp) {
  int l = 0;
  for (; l<L;) {
    W el = W();
    x.forEachVertexKey([&](auto u) {
      if (!fa(u)) return;
      louvainClearScanW(vcs, vcout);
      louvainScanCommunitiesW(vcs, vcout, x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
      if (c)      { louvainChangeCommunityW(vcom, ctot, x, u, c, vtot); fp(u); }
      el += e;  // l1-norm
    });
    if (fc(el, l++)) break;
  }
  return l;
}

#ifdef OPENMP
template <class G, class K, class W, class FC, class FA, class FP>
inline int louvainMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc, FA fa, FP fp) {
  size_t S = x.span();
  int l = 0;
  for (; l<L;) {
    W el = W();
    #pragma omp parallel for schedule(auto) reduction(+:el)
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!fa(u)) continue;
      louvainClearScanW(*vcs[t], *vcout[t]);
      louvainScanCommunitiesW(*vcs[t], *vcout[t], x, u, vcom);
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (c)      { louvainChangeCommunityOmpW(vcom, ctot, x, u, c, vtot); fp(u); }
      el += e;  // l1-norm
    }
    if (fc(el, l++)) break;
  }
  return l;
}
#endif


template <class G, class K, class W, class FC>
inline int louvainMoveW(vector<K>& vcom, vector<W>& ctot, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc) {
  auto fa = [](auto u) { return true; };
  auto fp = [](auto u) {};
  return louvainMoveW(vcom, ctot, vcs, vcout, x, vtot, M, R, L, fc, fa, fp);
}

#ifdef OPENMP
template <class G, class K, class W, class FC>
inline int louvainMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc) {
  auto fa = [](auto u) { return true; };
  auto fp = [](auto u) {};
  return louvainMoveOmpW(vcom, ctot, vcs, vcout, x, vtot, M, R, L, fc, fa, fp);
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
 * @param co csr offsets for vertices belonging to each community
 * @param ce csr data vertices belonging to each community
 */
template <class G, class K, class W>
inline void louvainAggregateW(G& a, vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  size_t S = x.span();
  a.respan(S);
  for (K c=0; c<S; ++c) {
    K oc = co[c];
    K nc = co[c+1] - co[c];
    if (nc==0) continue;
    louvainClearScanW(vcs, vcout);
    for (K i=0; i<nc; ++i) {
      K u = ce[oc+i];
      louvainScanCommunitiesW<true>(vcs, vcout, x, u, vcom);
    }
    a.addVertex(c);
    for (auto d : vcs)
      a.addEdge(c, d, vcout[d]);
  }
  a.update();
}

#ifdef OPENMP
template <class G, class K, class W>
inline void louvainAggregateOmpW(G& a, vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  size_t S = x.span();
  a.respan(S);
  #pragma omp parallel for schedule(auto)
  for (K c=0; c<S; ++c) {
    int t = omp_get_thread_num();
    K oc = co[c];
    K nc = co[c+1] - co[c];
    if (nc==0) continue;
    louvainClearScanW(*vcs[t], *vcout[t]);
    for (K i=0; i<nc; ++i) {
      K u = ce[oc+i];
      louvainScanCommunitiesW<true>(*vcs[t], *vcout[t], x, u, vcom);
    }
    for (auto d : *vcs[t])
      a.addEdge(c, d, (*vcout[t])[d]);
  }
  updateOmpU(a);
}
#endif


template <class G, class K, class W>
inline auto louvainAggregate(vector<K>& vcs, vector<W>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  G a; louvainAggregateW(a, vcs, vcout, x, vcom, co, ce);
  return a;
}

#ifdef OPENMP
template <class G, class K, class W>
inline auto louvainAggregateOmp(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, const G& x, const vector<K>& vcom, const vector<K>& co, const vector<K>& ce) {
  G a; louvainAggregateOmpW(a, vcs, vcout, x, vcom, co, ce);
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
template <class B, class G, class K, class V, class W>
inline auto louvainAffectedVerticesDeltaScreeningW(vector<K>& vcs, vector<W>& vcout, vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R=1) {
  fillValueU(vertices,    B());
  fillValueU(neighbors,   B());
  fillValueU(communities, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  for (size_t i=0; i<insertions.size();) {
    K u = get<0>(insertions[i]);
    louvainClearScanW(vcs, vcout);
    for (; i<insertions.size() && get<0>(insertions[i])==u; ++i) {
      K v = get<1>(insertions[i]);
      V w = get<2>(insertions[i]);
      if (vcom[u] == vcom[v]) continue;
      louvainScanCommunityW(vcs, vcout, u, v, w, vcom);
    }
    auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, vcs, vcout, M, R);
    if (e<=0) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[c] = 1;
  }
  x.forEachVertexKey([&](auto u) {
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  });
}


#ifdef OPENMP
template <class B, class G, class K, class V, class W>
inline auto louvainAffectedVerticesDeltaScreeningOmpW(vector<vector<K>*>& vcs, vector<vector<W>*>& vcout, vector<B>& vertices, vector<B>& neighbors, vector<B>& communities, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R=1) {
  size_t S = x.span();
  size_t D = deletions.size();
  size_t I = insertions.size();
  fillValueOmpU(vertices,    B());
  fillValueOmpU(neighbors,   B());
  fillValueOmpU(communities, B());
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
    neighbors[u] = 1;
    communities[vcom[v]] = 1;
  }
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int t = omp_get_thread_num();
    K  u0 = I>0? get<0>(insertions[0]) : 0;
    for (size_t i=0, n=0; i<I;) {
      K u = get<0>(insertions[i]);
      if (u!=u0) { ++n; u0 = u; }
      if (n % T != t) { ++i; continue; }
      louvainClearScanW(*vcs[t], *vcout[t]);
      for (; i<I && get<0>(insertions[i])==u; ++i) {
        K v = get<1>(insertions[i]);
        V w = get<2>(insertions[i]);
        if (vcom[u] == vcom[v]) continue;
        louvainScanCommunityW(*vcs[t], *vcout[t], u, v, w, vcom);
      }
      auto [c, e] = louvainChooseCommunity(x, u, vcom, vtot, ctot, *vcs[t], *vcout[t], M, R);
      if (e<=0) continue;
      vertices[u]  = 1;
      neighbors[u] = 1;
      communities[c] = 1;
    }
  }
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    if (neighbors[u]) x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; });
    if (communities[vcom[u]]) vertices[u] = 1;
  }
}
#endif




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
template <class B, class G, class K, class V>
inline void louvainAffectedVerticesFrontierW(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueU(vertices, B());
  for (const auto& [u, v] : deletions) {
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  for (const auto& [u, v, w] : insertions) {
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}


#ifdef OPENMP
template <class B, class G, class K, class V>
inline void louvainAffectedVerticesFrontierOmpW(vector<B>& vertices, const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& vcom) {
  fillValueOmpU(vertices, B());
  size_t D = deletions.size();
  size_t I = insertions.size();
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<D; ++i) {
    K u = get<0>(deletions[i]);
    K v = get<1>(deletions[i]);
    if (vcom[u] != vcom[v]) continue;
    vertices[u]  = 1;
  }
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<I; ++i) {
    K u = get<0>(insertions[i]);
    K v = get<1>(insertions[i]);
    if (vcom[u] == vcom[v]) continue;
    vertices[u]  = 1;
  }
}
#endif




// LOUVAIN
// -------

template <class G, class K, class FM, class FA, class FP>
auto louvainSeq(const G& x, const vector<K>* q, const LouvainOptions& o, FM fm, FA fa, FP fp) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t S = x.span();
  double M = edgeWeight(x)/2;
  vector<K> vcom(S), vcs, a(S);
  vector<W> vtot(S), ctot(S), vcout(S);
  vector<K> co(S+1), ce(S), cn(S);
  float tm = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el < E; };
    G y; y.respan(S);
    fillValueU(vcom, K());
    fillValueU(vtot, W());
    fillValueU(ctot, W());
    mark([&]() {
      tm = measureDuration(fm);
      louvainVertexWeightsW(vtot, x);
      if (q) louvainInitializeFromW(vcom, ctot, x, vtot, *q);
      else   louvainInitializeW(vcom, ctot, x, vtot);
      for (l=0, p=0; M>0 && p<P;) {
        int m = 0;
        if (p==0) m = louvainMoveW(vcom, ctot, vcs, vcout, x, vtot, M, R, L, fc, fa, fp);
        else      m = louvainMoveW(vcom, ctot, vcs, vcout, y, vtot, M, R, L, fc);
        if (p==0) copyValuesW(a, vcom);
        else      louvainLookupCommunitiesU(a, vcom);
        l += m; ++p;
        if (m<=1 || p>=P) break;
        const G& g = p<=1? x : y;
        size_t gn = g.order();
        size_t yn = louvainCountCommunityVerticesW(cn, g, vcom);
        if (double(yn)/gn >= o.aggregationTolerance) break;
        louvainCommunityVerticesW(co, ce, cn, g, vcom);
        y = louvainAggregate(vcs, vcout, g, vcom, co, ce);
        fillValueU(vcom, K());
        fillValueU(vtot, W());
        fillValueU(ctot, W());
        louvainVertexWeightsW(vtot, y);
        louvainInitializeW(vcom, ctot, y, vtot);
        E /= o.toleranceDecline;
      }
    });
  }, o.repeat);
  return LouvainResult<K>(a, l, p, t, tm);
}

#ifdef OPENMP
template <class G, class K, class FM, class FA, class FP>
auto louvainOmp(const G& x, const vector<K>* q, const LouvainOptions& o, FM fm, FA fa, FP fp) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  int    T = omp_get_max_threads();
  vector<K> vcom(S), a(S);
  vector<W> vtot(S), ctot(S);
  vector<K> co(S+1), ce(S), cn(S);
  vector<K> bufk(T);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  louvainAllocateHashtablesW(vcs, vcout, S);
  float tm = 0;
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el < E; };
    G y; y.respan(S);
    fillValueOmpU(vcom, K());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    mark([&]() {
      tm = measureDuration(fm);
      louvainVertexWeightsOmpW(vtot, x);
      if (q) louvainInitializeFromOmpW(vcom, ctot, x, vtot, *q);
      else   louvainInitializeOmpW(vcom, ctot, x, vtot);
      for (l=0, p=0; M>0 && p<P;) {
        int m = 0;
        if (p==0) m = louvainMoveOmpW(vcom, ctot, vcs, vcout, x, vtot, M, R, L, fc, fa, fp);
        else      m = louvainMoveOmpW(vcom, ctot, vcs, vcout, y, vtot, M, R, L, fc);
        if (p==0) copyValuesW(a, vcom);
        else      louvainLookupCommunitiesOmpU(a, vcom);
        l += m; ++p;
        if (m<=1 || p>=P) break;
        const G& g = p<=1? x : y;
        size_t gn = g.order();
        size_t yn = louvainCountCommunityVerticesOmpW(cn, g, vcom);
        if (double(yn)/gn >= o.aggregationTolerance) break;
        louvainCommunityVerticesOmpW(co, ce, cn, bufk, g, vcom);
        y = louvainAggregateOmp(vcs, vcout, g, vcom, co, ce);
        fillValueOmpU(vcom, K());
        fillValueOmpU(vtot, W());
        fillValueOmpU(ctot, W());
        louvainVertexWeightsOmpW(vtot, y);
        louvainInitializeOmpW(vcom, ctot, y, vtot);
        E /= o.toleranceDecline;
      }
    });
  }, o.repeat);
  louvainFreeHashtablesW(vcs, vcout);
  return LouvainResult<K>(a, l, p, t, tm);
}
#endif




// LOUVAIN-STATIC
// --------------

template <class G, class K>
inline auto louvainStaticSeq(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  auto fm = []() {};
  auto fa = [](auto u) { return true; };
  auto fp = [](auto u) {};
  return louvainSeq(x, q, o, fm, fa, fp);
}

#ifdef OPENMP
template <class G, class K>
inline auto louvainStaticOmp(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  auto fm = []() {};
  auto fa = [](auto u) { return true; };
  auto fp = [](auto u) {};
  return louvainOmp(x, q, o, fm, fa, fp);
}
#endif




// LOUVAIN DYNAMIC DELTA-SCREENING
// -------------------------------

template <class FLAG=char, class G, class K, class V>
inline auto louvainDynamicDeltaScreeningSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = x.span();
  double R = o.resolution;
  double M = edgeWeight(x)/2;
  const vector<K>& vcom = *q;
  vector<V> vtot(S), ctot(S);
  vector<K> vcs; vector<W> vcout(S);
  vector<B> vertices(S), neighbors(S), communities(S);
  louvainVertexWeightsW(vtot, x);
  louvainCommunityWeightsW(ctot, x, vcom, vtot);
  auto fm = [&]() { louvainAffectedVerticesDeltaScreeningW(vcs, vcout, vertices, neighbors, communities, x, deletions, insertions, vcom, vtot, ctot, M, R); vcs.clear(); vcout.clear(); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [](auto u) {};
  return louvainSeq(x, q, o, fm, fa, fp);
}


#ifdef OPENMP
template <class FLAG=char, class G, class K, class V>
inline auto louvainDynamicDeltaScreeningOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = FLAG;
  size_t S = x.span();
  double R = o.resolution;
  double M = edgeWeightOmp(x)/2;
  int    T = omp_get_max_threads();
  const vector<K>& vcom = *q;
  vector<W> vtot(S), ctot(S);
  vector<B> vertices(S), neighbors(S), communities(S);
  vector<vector<K>*> vcs(T);
  vector<vector<W>*> vcout(T);
  louvainAllocateHashtablesW(vcs, vcout, S);
  louvainVertexWeightsOmpW(vtot, x);
  louvainCommunityWeightsOmpW(ctot, x, vcom, vtot);
  auto fm = [&]() { louvainAffectedVerticesDeltaScreeningOmpW(vcs, vcout, vertices, neighbors, communities, x, deletions, insertions, vcom, vtot, ctot, M, R); louvainFreeHashtablesW(vcs, vcout); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [](auto u) {};
  return louvainOmp(x, q, o, fm, fa, fp);
}
#endif




// LOUVAIN DYNAMIC FRONTIER
// ------------------------

template <class FLAG=char, class G, class K, class V>
inline auto louvainDynamicFrontierSeq(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  B = FLAG;
  size_t S = x.span();
  const vector<K>& vcom = *q;
  vector<B> vertices(S);
  auto fm = [&]() { louvainAffectedVerticesFrontierW(vertices, x, deletions, insertions, vcom); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; }); };
  return louvainSeq(x, q, o, fm, fa, fp);
}


#ifdef OPENMP
template <class FLAG=char, class G, class K, class V>
inline auto louvainDynamicFrontierOmp(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  using  B = FLAG;
  size_t S = x.span();
  const vector<K>& vcom = *q;
  vector<B> vertices(S);
  auto fm = [&]() { louvainAffectedVerticesFrontierOmpW(vertices, x, deletions, insertions, vcom); };
  auto fa = [&](auto u) { return vertices[u]==B(1); };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vertices[v] = 1; }); };
  return louvainOmp(x, q, o, fm, fa, fp);
}
#endif
