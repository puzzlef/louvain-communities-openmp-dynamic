#pragma once
#include <cstddef>
#include <omp.h>
#include "_debug.hxx"




// BELONGS
// -------
// Check if work belongs to current thread.

template <class K>
inline bool belongsOmp(K key, int thread, int THREADS) {
  const K CHUNK_SIZE = 1024;
  K chunk = key / CHUNK_SIZE;
  return chunk % THREADS == thread;
}
template <class K>
inline bool belongsOmp(K key) {
  int thread  = omp_get_thread_num();
  int THREADS = omp_get_num_threads();
  return belongsOmp(key, thread, THREADS);
}




// FILL VALUE
// ----------

template <class T, class V>
void fillValueOmpU(T *a, size_t N, const V& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = v;
}
template <class T, class V>
inline void fillValueOmpU(vector<T>& a, const V& v) {
  fillValueOmpU(a.data(), a.size(), v);
}
template <class T, class V>
inline void fillValueOmpU(vector<T>& a, size_t i, size_t N, const V& v) {
  fillValueOmpU(a.data()+i, N, v);
}




// COPY VALUES
// -----------

template <class T, class TA>
inline size_t copyValuesOmp(const T *x, TA *a, size_t N) {
  ASSERT(x && a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
  return N;
}
template <class T, class TA>
inline size_t copyValuesOmp(const vector<T>& x, vector<TA>& a) {
  return copyValuesOmp(x.data(), a.data(), x.size());
}
template <class T, class TA>
inline size_t copyValuesOmp(const vector<T>& x, vector<TA>& a, size_t i, size_t N) {
  return copyValuesOmp(x.data()+i, a.data()+i, N);
}


template <class T, class TA>
inline size_t copyValuesOmpW(TA *a, const T *x, size_t N) {
  ASSERT(a && x);
  return copyValuesOmp(x, a, N);
}
template <class T, class TA>
inline size_t copyValuesOmpW(vector<TA>& a, const vector<T>& x) {
  return copyValuesOmp(x, a);
}
template <class T, class TA>
inline size_t copyValuesOmpW(vector<TA>& a, const vector<T>& x, size_t i, size_t N) {
  return copyValuesOmp(x, a, i, N);
}
