//------------------------------------------------------------------------------
//
//  Group theoretical algorithms
//
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

#include "perms.hpp"

#define SHOWCOSETREPS

template<typename T>
using gens_t = vector<Permutation<T>>;

template<typename T>
using compact_orbit_t = set<T>;

template<typename T>
using orbit_t = map<T, Permutation<T>>;

// shreier vector definition:
// v[num] = -1
// v[elt not in num orbit] = 0
// v[newelem] = i, where i is (#newgen + 1)
using shreier_t = vector<int>;

template<typename GenIt>
void print_simple(GenIt gbeg, GenIt gend) {
  for (auto git = gbeg; git != gend; ++git)
    cout << *git << endl;
}

template<typename OrbIt>
void print_orb(OrbIt obeg, OrbIt oend) {
  for (auto oit = obeg; oit != oend; ++oit)
    cout << oit->first << ": " << oit->second << endl;
  cout << endl;
}

// ref: HCGT, page 78
// calculates orbit for 'a' in compact form (i.e. as a set)
template<typename T, typename RandIt>
auto orbit(T num, RandIt gensbeg, RandIt gensend) {
  compact_orbit_t<T> Delta; 
  vector<T> DeltaNext{num};

  while (!DeltaNext.empty()) {
    vector<T> tmp {};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto elem : DeltaNext)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto newelem = igen->apply(elem);
        if (Delta.count(newelem) == 0)
          tmp.push_back(newelem);
      }
    DeltaNext.swap(tmp);
  }

  return Delta;
}

// ref: HCGT, page 79
// calculates orbit for 'a' in form { elt => u(elt) }, so, that a^u(elt) = elt
// also return generating set Y for stabilizer
template<typename T, typename RandIt>
auto orbit_stab(T num, RandIt gensbeg, RandIt gensend) {
  orbit_t<T> Delta, DeltaNext {{ num, {} }};
  gens_t<T> Y;

  while (!DeltaNext.empty()) {
    orbit_t<T> tmp {};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto&& [elem, curgen] : DeltaNext)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto newelem = igen->apply(elem); // beta ^ x
        auto newgen = product(curgen, *igen); // u_beta * x
        auto deltait = Delta.find(newelem); // u_(beta ^ x)
        if (deltait == Delta.end())
          tmp.insert({newelem, newgen});
        else
          Y.push_back(product(newgen, invert(deltait->second)));
      }
    DeltaNext.swap(tmp);
  }

  return make_pair(Delta, Y);
}

// ref: HCGT, page 80
// calculates orbit for 'a' in compact form
// returns pair of orbit and shreier vector
template<typename T, typename RandIt> 
auto orbit_shreier(T num, RandIt gensbeg, RandIt gensend) {
  compact_orbit_t<T> Delta; 
  vector<T> DeltaNext {num};
  shreier_t v(T::fin - T::start + 1);
  v[num - T::start] = -1;

  while (!DeltaNext.empty()) {
    vector<T> tmp {};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto elem : DeltaNext)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto genidx = igen - gensbeg;
        auto newelem = igen->apply(elem);
        if (Delta.count(newelem) == 0) {
          tmp.push_back(newelem);
          v[newelem - T::start] = genidx + 1;
        }
      }    
    DeltaNext.swap(tmp);
  }

  return make_pair(Delta, v);
}
