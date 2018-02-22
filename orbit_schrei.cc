//------------------------------------------------------------------------------
//
// 3. Better orbit calculation with schreier vector
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
using std::pair;
using std::ostream_iterator;
using std::vector;

template<typename T>
using orbit_t = set<T>;

template<typename T>
using schreier_t = vector<int>;

template<typename T>
using gens_t = vector<Permutation<T>>;

template<typename T> auto
orbit_shreier(T num, T DomainStart, T DomainEnd, gens_t<T> gens) {
  orbit_t<T> Delta, DeltaNext;
  schreier_t<T> Schreier(DomainEnd - DomainStart);
  Schreier[num - DomainStart] = -1;
  DeltaNext.insert(num);
  while (!DeltaNext.empty()) {
    orbit_t<T> tmp {};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto elem : DeltaNext) {      
      for (size_t genidx = 0; genidx != gens.size(); ++genidx) {
        auto gen = gens[genidx];
        auto newelem = gen.apply(elem);
        if (Delta.count(newelem) == 0) {
          tmp.insert(newelem);
          Schreier[newelem - DomainStart] = genidx + 1;
        }
      }
    }
    DeltaNext.swap(tmp);
  }

  return make_pair(Delta, Schreier);
}

template<typename T> auto
ubeta (T num, T orbelem, T DomainStart, T DomainEnd, 
       schreier_t<T> schrei, gens_t<T> gens) {
  Permutation<T> res (DomainStart, DomainEnd);
  auto k = schrei[orbelem - DomainStart];
  if (k == 0)
    return make_pair(res, false);
  
  while (k != -1) {
    assert (k != 0); // we can not normally go away from orbit unwinding
    res.lmul(gens[k-1]);
    auto invgen = invert(gens[k-1]);
    orbelem = invgen.apply(orbelem);
    k = schrei[orbelem - DomainStart];
  }
  return make_pair(res, true);
}

int
main () {
  unsigned DomainStart = 1;
  unsigned DomainEnd = 8;

  unsigned interesting = 1;

  Permutation<unsigned> g1 (DomainStart, DomainEnd, {{1, 3, 7},{2, 5}});
  Permutation<unsigned> g2 (DomainStart, DomainEnd, {{3, 4, 6, 7}});
  auto [orb, schrei] = orbit_shreier(interesting, DomainStart, DomainEnd, {g1, g2});

  cout << "Orbit for " << interesting << " is: " << endl;
  for (auto elem : orb)
    cout << elem << " " 
         << ubeta(interesting, elem, DomainStart, DomainEnd, schrei, {g1, g2}).first 
         << endl;

  cout << "Schreier vector is: ";
  for (auto v : schrei)
    cout << v << " ";
  cout << endl;
}
