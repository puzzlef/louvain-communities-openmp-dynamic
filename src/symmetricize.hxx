#pragma once




// SYMMETRICIZE
// ------------

template <bool CHK=false, class G>
void symmetricizeU(G& a) {
  a.forEachVertexKey([&](auto u) {
    if (!CHK) a.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
    else a.forEachEdge(u, [&](auto v, auto w) { if (!a.hasEdge(v, u)) a.addEdge(v, u, w); });
  });
  a.correct(CHK);
}

template <class H, class G>
void symmetricizeW(H& a, const G& x) {
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      a.addEdge(u, v, w);
      a.addEdge(v, u, w);
    });
  });
  a.correct();
}

template <class G>
auto symmetricize(const G& x) {
  G a; symmetricizeW(a, x);
  return a;
}
