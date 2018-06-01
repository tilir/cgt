//------------------------------------------------------------------------------
//
// Performance tests for groups, etc
//
//------------------------------------------------------------------------------
// g++ --std=c++17 -Wfatal-errors perftests.cc -O3 -DNDEBUG -S -fno-exceptions
// -fno-rtti

#include "groups.hpp"
#include "idomain.hpp"

using orbits::DirectOrbit;
using orbits::ShreierOrbit;
using permutations::PermLoop;
using namespace groups;

template <template <class> class Orb, typename T> bool perftest_orbit(T elt) {
  auto result = true;
  auto gens = symmetric_gens<T>();
  Orb<T> orbit(elt, gens.begin(), gens.end());
  for (auto &&beta : orbit) {
    auto u_beta = orbit.ubeta(beta);
    result = result && (u_beta.apply(elt) == beta);
  }
  return result;
}

#ifndef DORBC
#define DORBC 1
#endif

#ifndef SORBC
#define SORBC 1
#endif

// given DORBC = 1, DORBS grows as: 1000 => 0.74, 1500 => 1.7, 2000 => 5.2
#ifndef DORBS
#define DORBS 2000
#endif

// given SORBC = 1, SORBS grows as: 200 => 0.63, 300 => 2.4, 400 => 8.8
// reason is we have long and increasing unwinding in ubeta function
// say with gens: (1, 2) and (1, 2, 3, 4, 5, 6, 7, 8, 9)
// going from 1 to 9 involves 9 applications and sv looks like:
// [-1, 1, 2, 2, 2, 2, 2, 2, 2]
#ifndef SORBS
#define SORBS 400
#endif

int main() {
  // some cache warmup
  UnsignedDomain<1, 1000> elt = 1;
  perftest_orbit<DirectOrbit>(elt);
  UnsignedDomain<1, 200> selt = 1;
  perftest_orbit<ShreierOrbit>(selt);

  bool res = true;
  cout << "direct orbit: ";
  auto tdorb = duration([&] {
    for (int x = 1; x <= DORBC; ++x) {
      UnsignedDomain<1, DORBS> elt = x;
      res = res && perftest_orbit<DirectOrbit>(elt);
    }
  });
  cout << tdorb.count() << ", " << res << endl;

  res = true;
  cout << "shreier orbit: ";
  auto tsorb = duration([&] {
    for (int x = 1; x <= SORBC; ++x) {
      UnsignedDomain<1, SORBS> elt = x;
      res = res && perftest_orbit<ShreierOrbit>(elt);
    }
  });
  cout << tsorb.count() << ", " << res << endl;
}
