//------------------------------------------------------------------------------
//
//  Common includes and defines for permutations
//
//------------------------------------------------------------------------------
//
// Permutation support extensively uses standard library. Most common uses
// collected here
//
//------------------------------------------------------------------------------

#ifndef PERMCOMMON_GUARD_
#define PERMCOMMON_GUARD_

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using std::array;
using std::cerr;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;
using std::cout;
using std::endl;
using std::exception;
using std::find;
using std::forward;
using std::initializer_list;
using std::iota;
using std::logic_error;
using std::make_pair;
using std::make_reverse_iterator;
using std::make_tuple;
using std::map;
using std::min_element;
using std::move;
using std::next;
using std::ostream;
using std::ostream_iterator;
using std::pair;
using std::prev;
using std::queue;
using std::rbegin;
using std::rend;
using std::reverse;
using std::rotate;
using std::runtime_error;
using std::set;
using std::sort;
using std::string;
using std::stringstream;
using std::swap;
using std::transform;
using std::vector;

template <typename T> struct reversion_wrapper { T &iterable; };

template <typename T> auto begin(reversion_wrapper<T> w) {
  return rbegin(w.iterable);
}

template <typename T> auto end(reversion_wrapper<T> w) {
  return rend(w.iterable);
}

template <typename T> reversion_wrapper<T> reverse(T &&iterable) {
  return {iterable};
}

template <typename F, typename... Args>
auto duration(F &&func, Args &&... args) {
  auto start = steady_clock::now();
  forward<decltype(func)>(func)(forward<Args>(args)...);
  return duration_cast<milliseconds>(steady_clock::now() - start);
}

#endif
