#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include "_debug.hxx"

using std::vector;
using std::abs;
using std::max;
using std::fill;




// VECTOR 2D
// ---------

template <class T>
using vector2d = vector<vector<T>>;




// GATHER VALUES
// -------------

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[j++] = TA(fm(x[i]));
}
template <class TA, class TX, class IS>
inline void gatherValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesW(a.data(), x.data(), is);
}

template <class IS>
inline void gatherValuesW(vector<bool>& a, const vector<bool>& x, const IS& is) {
  size_t j = 0;
  for (auto i : is)
    a[j++] = x[i];
}


#ifdef OPENMP
template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[j] = TA(fm(x[is[j]]));
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesOmpW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesOmpW(a.data(), x.data(), is);
}
#endif




// SCATTER VALUES
// --------------

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t j = 0;
  for (auto i : is)
    a[i] = TA(fm(x[j++]));
}
template <class TA, class TX, class IS>
inline void scatterValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesW(a.data(), x.data(), is);
}


#ifdef OPENMP
template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] = TA(fm(x[j]));
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesOmpW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesOmpW(a.data(), x.data(), is);
}
#endif




// FILL VALUE
// ----------

template <class T>
inline void fillValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  fill(a, a+N, v);
}
template <class T>
inline void fillValueU(vector<T>& a, const T& v) {
  fill(a.begin(), a.end(), v);
}


#ifdef OPENMP
template <class T>
inline void fillValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = v;
}
template <class T>
inline void fillValueOmpU(vector<T>& a, const T& v) {
  fillValueOmpU(a.data(), a.size(), v);
}
#endif




// COPY VALUES
// -----------

template <class TA, class TX>
inline void copyValuesW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

template <class TA, class TX>
inline void copyValuesW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesW(a.data(), x.data(), x.size());
}


#ifdef OPENMP
template <class TA, class TX>
inline void copyValuesOmpW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

template <class TA, class TX>
inline void copyValuesOmpW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesOmpW(a.data(), x.data(), x.size());
}
#endif




// MULTIPLY VALUES
// ---------------

template <class TA, class TX, class TY>
inline void multiplyValuesW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

template <class TA, class TX, class TY>
inline void multiplyValuesW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesW(a.data(), x.data(), y.data(), x.size());
}
template <class TA, class TX, class TY>
inline void multiplyValuesW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y, size_t i, size_t N) {
  multiplyValuesW(a.data()+i, x.data()+i, y.data()+i, N);
}


#ifdef OPENMP
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesOmpW(a.data(), x.data(), y.data(), x.size());
}
template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y, size_t i, size_t N) {
  multiplyValuesOmpW(a.data()+i, x.data()+i, y.data()+i, N);
}
#endif




// L1-NORM
// -------

template <class TX, class TA=TX>
inline TA l1Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

template <class TX, class TA=TX>
inline TA l1Norm(const vector<TX>& x, TA a=TA()) {
  return l1Norm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l1Norm(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l1Norm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1Norm(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA l1Norm(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return l1Norm(x.data()+i, y.data()+i, N, a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA l1NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

template <class TX, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, TA a=TA()) {
  return l1NormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l1NormOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormOmp(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return l1NormOmp(x.data()+i, y.data()+i, N, a);
}
#endif




// L2-NORM
// -------

template <class TX, class TA=TX>
inline TA l2Norm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

template <class TX, class TA=TX>
inline TA l2Norm(const vector<TX>& x, TA a=TA()) {
  return l2Norm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l2Norm(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l2Norm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2Norm(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA l2Norm(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return l2Norm(x.data()+i, y.data()+i, N, a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA l2NormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

template <class TX, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, TA a=TA()) {
  return l2NormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l2NormOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(+:a)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormOmp(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return l2NormOmp(x.data()+i, y.data()+i, N, a);
}
#endif




// LI-NORM
// -------

template <class TX, class TA=TX>
inline TA liNorm(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

template <class TX, class TA=TX>
inline TA liNorm(const vector<TX>& x, TA a=TA()) {
  return liNorm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA liNorm(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA liNorm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNorm(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA liNorm(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return liNorm(x.data()+i, y.data()+i, N, a);
}


#ifdef OPENMP
template <class TX, class TA=TX>
inline TA liNormOmp(const TX *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

template <class TX, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, TA a=TA()) {
  return liNormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA liNormOmp(const TX *x, const TY *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto) reduction(max:a)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormOmp(x.data(), y.data(), x.size(), a);
}
template <class TX, class TY, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, const vector<TY>& y, size_t i, size_t N, TA a=TA()) {
  return liNormOmp(x.data()+i, y.data()+i, N, a);
}
#endif
