#pragma once
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"

#ifdef OPENMP
#include <omp.h>
#endif

using std::tuple;
using std::string;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::move;
using std::max;
using std::getline;




// READ MTX
// --------

inline size_t readMtxHeader(istream& s, bool& symmetric, size_t& rows, size_t& cols, size_t& size) {
  string line, h0, h1, h2, h3, h4;
  // Skip past the comments and read the graph type.
  while (true) {
    getline(s, line);
    if (line.find('%')!=0) break;
    if (line.find("%%")!=0) continue;
    istringstream sline(line);
    sline >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  if (h1!="matrix" || h2!="coordinate") return 0;
  symmetric = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  istringstream sline(line);
  sline >> rows >> cols >> size;
  return max(rows, cols);
}


template <class G>
inline void readMtxW(G& a, istream& s, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  bool symmetric; size_t rows, cols, size;
  size_t n = readMtxHeader(s, symmetric, rows, cols, size);
  if (n==0) return;
  PERFORMI( auto t0 = timeNow() );
  // Add all vertices first.
  a.reserve(n+1);
  for (size_t u=1; u<=n; ++u)
    a.addVertex(K(u));
  PERFORMI( auto t1 = timeNow() );
  // Then we add the edges.
  string line;
  while (getline(s, line)) {
    size_t u, v; double w = 1;
    istringstream sline(line);
    if (!(sline >> u >> v)) break;
    if (weighted) sline >> w;
    a.addEdge(K(u), K(v), E(w));
    if (symmetric) a.addEdge(K(v), K(u), E(w));
  }
  PERFORMI( auto t2 = timeNow() );
  // Apply graph update.
  a.update();
  PERFORMI( auto t3 = timeNow() );
  PERFORMI( float dvertices = duration(t0, t1) );
  PERFORMI( float dedges    = duration(t1, t2) );
  PERFORMI( float dupdate   = duration(t2, t3) );
  LOGI("readMtxW(): vertices=%.1fms, edges=%.1fms, update=%.1fms\n", dvertices, dedges, dupdate);
}
template <class G>
inline void readMtxW(G& a, const char *pth, bool weighted=false) {
  ifstream f(pth);
  readMtxW(a, f, weighted);
}


#ifdef OPENMP
template <class G>
inline void readMtxOmpW(G& a, istream& s, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  bool symmetric; size_t rows, cols, size;
  size_t n = readMtxHeader(s, symmetric, rows, cols, size);
  if (n==0) return;
  PERFORMI( auto t0 = timeNow() );
  // Add all vertices first.
  a.reserve(n+1);
  for (size_t u=1; u<=n; ++u)
    a.addVertex(u);
  PERFORMI( auto t1 = timeNow() );
  // Then we add the edges.
  const int THREADS = omp_get_max_threads();
  const int LINES   = 131072;
  const size_t CHUNK_SIZE = 1024;
  vector<string> lines(LINES);
  vector<tuple<size_t, size_t, double>> edges(LINES);
  PERFORMI( float dread  = 0 );
  PERFORMI( float dparse = 0 );
  PERFORMI( float dedges = 0 );
  while (true) {
    PERFORMI( auto t2 = timeNow() );
    // Read several lines from the stream.
    int READ = 0;
    for (int i=0; i<LINES; ++i, ++READ)
      if (!getline(s, lines[i])) break;
    if (READ==0) break;
    PERFORMI( auto t3 = timeNow() );
    // Parse lines using multiple threads.
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i=0; i<READ; ++i) {
      char *line = (char*) lines[i].c_str();
      size_t u = strtoull(line, &line, 10);
      size_t v = strtoull(line, &line, 10);
      double w = weighted? strtod(line, &line) : 0;
      edges[i] = {u, v, w? w : 1};
    }
    PERFORMI( auto t4 = timeNow() );
    // Add edges to the graph.
    #pragma omp parallel
    {
      for (int i=0; i<READ; ++i) {
        const auto& [u, v, w] = edges[i];
        addEdgeOmpU(a, K(u), K(v), E(w));
        if (symmetric) addEdgeOmpU(a, K(v), K(u), E(w));
      }
    }
    PERFORMI( auto t5 = timeNow() );
    PERFORMI( dread  += duration(t2, t3) );
    PERFORMI( dparse += duration(t3, t4) );
    PERFORMI( dedges += duration(t4, t5) );
  }
  PERFORMI( auto t6 = timeNow() );
  updateOmpU(a);
  PERFORMI( auto t7 = timeNow() );
  PERFORMI( float dvertices = duration(t0, t1) );
  PERFORMI( float dupdate   = duration(t6, t7) );
  LOGI("readMtxOmpW(): vertices=%.1fms, read=%.1fms, parse=%.1fms, edges=%.1fms, update=%.1fms\n", dvertices, dread, dparse, dedges, dupdate);
}
template <class G>
inline void readMtxOmpW(G& a, const char *pth, bool weighted=false) {
  ifstream f(pth);
  readMtxOmpW(a, f, weighted);
}
#endif



// WRITE MTX
// ---------

template <class G>
inline void writeMtx(ostream& a, const G& x) {
  a << "%%MatrixMarket matrix coordinate real asymmetric\n";
  a << x.order() << " " << x.order() << " " << x.size() << "\n";
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      a << u << " " << v << " " << w << "\n";
    });
  });
}
template <class G>
inline void writeMtx(string pth, const G& x) {
  ofstream f(pth);
  writeMtx(f, x);
  f.close();
}
