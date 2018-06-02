//------------------------------------------------------------------------------
//
//  Group generators for finitely generated ones
//
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

#ifndef GROUPGENS_GUARD_
#define GROUPGENS_GUARD_

#include "perms.hpp"

using permutations::PermLoop;
using permutations::Permutation;

namespace groupgens {

// generators for permutation group
template <typename T> using gens_t = vector<Permutation<T>>;

// cyclic group is like {{(1 2 3 4 5)}}
template <typename Domain> auto cyclic_gens() {
  vector<Domain> init(Domain::fin - Domain::start + 1);
  iota(init.begin(), init.end(), Domain::start);
  PermLoop<Domain> cyclic(init.begin(), init.end());
  gens_t<Domain> retval{{cyclic}};
  return retval;
}

// symmetric group is like {{(1 2)}, {(1 3)}, {(1 4)}, {(1, 5)}}
template <typename Domain> auto symmetric_gens() {
  assert(Domain::fin - Domain::start > 1);
  gens_t<Domain> retval;
  for (typename Domain::type i = 1; i < Domain::fin - Domain::start + 1; ++i)
    retval.push_back({{Domain::start, Domain::start + i}});
  return retval;
}

// symmetric group with minimal set is like {{(1 2)}, {(1 2 3 4 5)}}
template <typename Domain> auto min_symmetric_gens() {
  assert(Domain::fin - Domain::start > 1);
  gens_t<Domain> retval = cyclic_gens<Domain>();
  PermLoop<Domain> pair = {Domain::start, Domain::start + 1};
  retval.push_back({pair});
  return retval;
}

// alternating group is like {{(1 2 3)}, {(1 2 4)}, {(1 2 5)}}
template <typename Domain> auto alternating_gens() {
  assert(Domain::fin - Domain::start > 2);
  gens_t<Domain> retval;
  for (typename Domain::type i = 2; i < Domain::fin - Domain::start + 1; ++i)
    retval.push_back({{Domain::start, Domain::start + 1, Domain::start + i}});
  return retval;
}

} // namespace groupgens

#endif
