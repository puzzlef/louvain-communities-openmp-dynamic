#pragma once
#include <algorithm>
#include <vector>
#include "_debug.hxx"
#include "_algorithm.hxx"

using std::vector;
using std::fill;




// 2D
// --

template <class T>
using vector2d = vector<vector<T>>;




// SIZE
// ----

template <class T>
inline size_t size(const vector<T>& x) {
  return x.size();
}




// GATHER VALUES
// -------------

template <class T, class J, class TA, class FM>
inline void gatherValues(const T *x, const J& is, TA *a, FM fm) {
  ASSERT(x && a);
  size_t j = 0;
  for (auto i : is)
    a[j++] = TA(fm(x[i]));
}
template <class T, class J, class TA>
inline void gatherValues(const T *x, const J& is, TA *a) {
  ASSERT(x && a);
  auto fm = [](auto v) { return v; };
  gatherValues(x, is, a, fm);
}
template <class T, class J, class TA, class FM>
inline void gatherValues(const vector<T>& x, const J& is, vector<TA>& a, FM fm) {
  gatherValues(x.data(), is, a.data(), fm);
}
template <class T, class J, class TA>
inline void gatherValues(const vector<T>& x, const J& is, vector<TA>& a) {
  gatherValues(x.data(), is, a.data());
}

template <class T, class J, class TA, class FM>
inline void gatherValuesW(TA *a, const T *x, const J& is, FM fm) {
  ASSERT(a && x);
  gatherValues(x, is, a, fm);
}
template <class T, class J, class TA>
inline void gatherValuesW(TA *a, const T *x, const J& is) {
  ASSERT(a && x);
  gatherValues(x, is, a);
}
template <class T, class J, class TA, class FM>
inline void gatherValuesW(vector<TA>& a, const vector<T>& x, const J& is, FM fm) {
  gatherValues(x, is, a, fm);
}
template <class T, class J, class TA>
inline void gatherValuesW(vector<TA>& a, const vector<T>& x, const J& is) {
  gatherValues(x, is, a);
}




// SCATTER VALUES
// --------------

template <class T, class J, class TA, class FM>
inline void scatterValues(const T *x, const J& is, TA *a, FM fm) {
  ASSERT(x && a);
  size_t j = 0;
  for (auto i : is)
    a[i] = fm(x[j++]);
}
template <class T, class J, class TA>
inline void scatterValues(const T *x, const J& is, TA *a) {
  ASSERT(x && a);
  auto fm = [](auto v) { return v; };
  scatterValues(x, is, a, fm);
}
template <class T, class J, class TA, class FM>
inline void scatterValues(const vector<T>& x, const J& is, vector<TA>& a, FM fm) {
  scatterValues(x.data(), is, a.data(), fm);
}
template <class T, class J, class TA>
inline void scatterValues(const vector<T>& x, const J& is, vector<TA>& a) {
  scatterValues(x.data(), is, a.data());
}

template <class T, class J, class TA, class FM>
inline void scatterValuesW(TA *a, const T *x, const J& is, FM fm) {
  ASSERT(a && x);
  scatterValues(x, is, a, fm);
}
template <class T, class J, class TA>
inline void scatterValuesW(TA *a, const T *x, const J& is) {
  ASSERT(a && x);
  scatterValues(x, is, a);
}
template <class T, class J, class TA, class FM>
inline void scatterValuesW(vector<TA>& a, const vector<T>& x, const J& is, FM fm) {
  scatterValues(x, is, a, fm);
}
template <class T, class J, class TA>
inline void scatterValuesW(vector<TA>& a, const vector<T>& x, const J& is) {
  scatterValues(x, is, a);
}




// COPY VALUES
// -----------

template <class T, class TA>
inline size_t copyValues(const T *x, TA *a, size_t N) {
  ASSERT(x && a);
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
  return N;
}
template <class T, class TA>
inline size_t copyValues(const vector<T>& x, vector<TA>& a) {
  return copyValues(x.data(), a.data(), x.size());
}
template <class T, class TA>
inline size_t copyValues(const vector<T>& x, vector<TA>& a, size_t i, size_t N) {
  return copyValues(x.data()+i, a.data()+i, N);
}




// FILL VALUE
// ----------

template <class T, class V>
inline void fillValueU(T *a, size_t N, const V& v) {
  ASSERT(a);
  fill(a, a+N, v);
}
template <class T, class V>
inline void fillValueU(vector<T>& a, const V& v) {
  fill(a.begin(), a.end(), v);
}
template <class T, class V>
inline void fillValueU(vector<T>& a, size_t i, size_t N, const V& v) {
  fill(a.begin()+i, a.begin()+i+N, v);
}
