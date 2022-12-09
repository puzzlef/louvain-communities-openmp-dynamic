#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "properties.hxx"
#include "duplicate.hxx"
#include "modularity.hxx"
#include "louvain.hxx"

using std::tuple;
using std::vector;
using std::min;




// LOUVAIN-SEQ
// -----------

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




// LOUVAIN-SEQ-STATIC
// ------------------

template <class G, class K>
inline auto louvainSeqStatic(const G& x, const vector<K>* q=nullptr, const LouvainOptions& o={}) {
  return louvainSeq(x, q, o);
}




// LOUVAIN-SEQ-DYNAMIC-DELTA-SCREENING
// -----------------------------------

template <class G, class K, class V>
inline auto louvainSeqDynamicDeltaScreening(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
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




// LOUVAIN-SEQ-DYNAMIC-FRONTIER
// ----------------------------

template <class G, class K, class V>
inline auto louvainSeqDynamicFrontier(const G& x, const vector<tuple<K, K>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>* q, const LouvainOptions& o={}) {
  size_t S = x.span();
  const vector<K>& vcom = *q;
  auto vaff = louvainAffectedVerticesFrontier(x, deletions, insertions, vcom);
  auto fa = [&](auto u) { return vaff[u]==true; };
  auto fp = [&](auto u) { x.forEachEdgeKey(u, [&](auto v) { vaff[v] = true; }); };
  return louvainSeq(x, q, o, fa, fp);
}
