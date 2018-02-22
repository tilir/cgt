//------------------------------------------------------------------------------
//
// 2. Better orbit calculation with stabilizer added
//
//------------------------------------------------------------------------------

#include <iostream>
#include <iterator>
#include <map>
#include <vector>
#include "perms.hpp"

using std::cout;
using std::copy;
using std::endl;
using std::make_pair;
using std::map;
using std::pair;
using std::ostream_iterator;
using std::vector;

template<typename T>
using orbit_t = map<T, Permutation<T>>;

template<typename T>
using gens_t = vector<Permutation<T>>;

template<typename T> auto
orbit_stabilizer(T num, T DomainStart, T DomainEnd, gens_t<T> gens) {
  orbit_t<T> Delta;
  orbit_t<T> DeltaNext;
  Permutation<T> id (DomainStart, DomainEnd);
  DeltaNext[num] = id;
  gens_t<T> Stabilizer;

  while (!DeltaNext.empty()) {
    orbit_t<T> tmp {};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto [elem, curgen] : DeltaNext)
      for (auto gen : gens) {
        auto newelem = gen.apply(elem); // elem ^ gen
        auto newgen = product(curgen, gen);
        if (Delta.find(newelem) == Delta.end()) {
          tmp[newelem] = newgen;
        }
        else {
          auto stabgen = product(newgen, invert(Delta[newelem]));
          Stabilizer.push_back(stabgen);
        }
      }
    DeltaNext.swap(tmp);
  }

  return make_pair(Delta, Stabilizer);
}

int
main () {
  unsigned DomainStart = 1;
  unsigned DomainEnd = 8;

  unsigned interesting = 1;

  Permutation<unsigned> g1 (DomainStart, DomainEnd, {{1, 3, 7},{2, 5}});
  Permutation<unsigned> g2 (DomainStart, DomainEnd, {{3, 4, 6, 7}});
  auto [orb, stab] = orbit_stabilizer(interesting, DomainStart, DomainEnd, {g1, g2});

  cout << "Orbit for " << interesting << " is: " << endl;
  for (auto [elem, gen] : orb)
    cout << elem << " " << gen << endl;

  cout << "Stabilizer generators are: " << endl;
  for (auto gen : stab)
    cout << gen << endl;
}
