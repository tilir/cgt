//------------------------------------------------------------------------------
//
// HCGT, page 78. Simplest and far from ideal orbit calculation
//
//------------------------------------------------------------------------------

#include <cassert>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include "perms.hpp"

using std::cout;
using std::copy;
using std::endl;
using std::map;
using std::ostream_iterator;
using std::set;

template<typename T>
void orbit(T num, T DomainStart, T DomainEnd, vector<Permutation<T>> gens) {
  assert (!gens.empty());
  set<T> orbit;
  set<T> next = {num};
  map<T, Permutation<T>> reps;
  Permutation<T> id (DomainStart, DomainEnd);
  reps[num] = id;

  while (!next.empty()) {
    set<T> tmp {};
    orbit.insert(next.begin(), next.end());
    for (const auto &elem : next)
      for (auto gen : gens) {
        auto newelem = gen.apply(elem); // elem ^ gen
        if (orbit.count(newelem) == 0) {
          tmp.insert(newelem);
          reps[newelem] = product(reps[elem], gen); // T[delta] * gen
        }
      }
    next.swap(tmp);
  }

  cout << "Orbit for " << num << " is: [ ";
  copy(orbit.begin(), orbit.end(), ostream_iterator<T>(cout, " "));
  cout << "]" << endl;
  cout << "Coset representatives:" << endl;
  for (auto &g : reps)
    cout << g.first << ": " << g.second << endl;
}

int
main () {
  unsigned DomainStart = 1;
  unsigned DomainEnd = 6;

#ifdef ALT
  Permutation<unsigned> g1 (DomainStart, DomainEnd, {{1, 2},{3, 4}});
  Permutation<unsigned> g2 (DomainStart, DomainEnd, {{2, 3},{4, 5}});

  for (unsigned i = DomainStart; i != DomainEnd; ++i)
    orbit(i, DomainStart, DomainEnd, {g1, g2});
#endif

#ifdef SYM
  Permutation<unsigned> g3 (DomainStart, DomainEnd, {{1, 2}});
  Permutation<unsigned> g4 (DomainStart, DomainEnd, {{1, 2, 3, 4, 5}});

  for (unsigned i = DomainStart; i != DomainEnd; ++i)
    orbit(i, DomainStart, DomainEnd, {g3, g4});
#endif

  unsigned Domain7Start = 1;
  unsigned Domain7End = 8;

  Permutation<unsigned> g5 (Domain7Start, Domain7End, {{1, 3, 7},{2, 5}});
  Permutation<unsigned> g6 (Domain7Start, Domain7End, {{3, 4, 6, 7}});
  orbit(1u, Domain7Start, Domain7End, {g5, g6});
}
