#pragma once
#include <vector>
#include <unordered_map>

using std::vector;
using std::unordered_map;




// VALUES OF
// ---------

/**
 * Get values of an unordered map as a vector.
 * @param x unordered map
 */
template <class K, class V, class FM>
inline auto valuesOf(const unordered_map<K, V>& x, FM fm) {
  using W = decltype(fm(K(), V()));
  vector<W> a;
  for (const auto& [k, v] : x)
    a.push_back(fm(k, v));
  return a;
}
