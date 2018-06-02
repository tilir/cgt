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
using namespace groupgens;

//------------------------------------------------------------------------------
//
// 01: minimal generating set for symmetric group
//
//------------------------------------------------------------------------------

#ifndef DORBC_01
#define DORBC_01 1
#endif

#ifndef SORBC_01
#define SORBC_01 1
#endif

// given DORBC = 1, DORBS grows as: 1000 => 0.74, 1500 => 1.7, 2000 => 5.2
#ifndef DORBS_01
#define DORBS_01 2000
#endif

// given SORBC = 1, SORBS grows as: 200 => 0.63, 300 => 2.4, 400 => 8.8
// reason is we have long and increasing unwinding in ubeta function
// say with gens: (1, 2) and (1, 2, 3, 4, 5, 6, 7, 8, 9)
// going from 1 to 9 involves 9 applications and sv looks like:
// [-1, 1, 2, 2, 2, 2, 2, 2, 2]
#ifndef SORBS_01
#define SORBS_01 400
#endif

template <template <class> class Orb, typename T>
bool perftest_orbit_01(T elt) {
  auto result = true;
  auto gens = min_symmetric_gens<T>();
  Orb<T> orbit(elt, gens.begin(), gens.end());
  for (auto &&beta : orbit) {
    auto u_beta = orbit.ubeta(beta);
    result = result && (u_beta.apply(elt) == beta);
  }
  return result;
}

//------------------------------------------------------------------------------
//
// 02: shallow generating set for symmetric group
//
//------------------------------------------------------------------------------

#ifndef DORBC_02
#define DORBC_02 1
#endif

#ifndef SORBC_02
#define SORBC_02 1
#endif

// given DORBC = 1, DORBS grows as: 400 => 0.4, 800 => 2.5, 1000 => 4.8
// direct orbit significantly slows down compared to test01
#ifndef DORBS_02
#define DORBS_02 1000
#endif

// given SORBC = 1, SORBS grows as: 400 => 0.6, 800 => 4.6, 1000 => 8.7
// shreier orbit significantly speeds up compared to test01
// in shallow case ubeta chains are smaller
#ifndef SORBS_02
#define SORBS_02 800
#endif

template <template <class> class Orb, typename T>
bool perftest_orbit_02(T elt) {
  auto result = true;
  auto gens = symmetric_gens<T>();
  Orb<T> orbit(elt, gens.begin(), gens.end());
  for (auto &&beta : orbit) {
    auto u_beta = orbit.ubeta(beta);
    result = result && (u_beta.apply(elt) == beta);
  }
  return result;
}

int main() {
  // some cache warmup
  UnsignedDomain<1, 1000> elt = 1;
  perftest_orbit_01<DirectOrbit>(elt);
  UnsignedDomain<1, 200> selt = 1;
  perftest_orbit_01<ShreierOrbit>(selt);

  bool res = true;

// test 01: minimal generating set for symmetric group
#ifndef NOTEST_01
  cout << "direct orbit: ";
  auto tdorb_01 = duration([&] {
    for (int x = 1; x <= DORBC_01; ++x) {
      UnsignedDomain<1, DORBS_01> elt = x;
      res = res && perftest_orbit_01<DirectOrbit>(elt);
    }
  });
  cout << tdorb_01.count() << ", " << res << endl;

  res = true;
  cout << "shreier orbit: ";
  auto tsorb_01 = duration([&] {
    for (int x = 1; x <= SORBC_01; ++x) {
      UnsignedDomain<1, SORBS_01> elt = x;
      res = res && perftest_orbit_01<ShreierOrbit>(elt);
    }
  });
  cout << tsorb_01.count() << ", " << res << endl;
#endif

// test 02: shallow generating set for symmetric group
#ifndef NOTEST_02
  cout << "direct orbit: ";
  auto tdorb_02 = duration([&] {
    for (int x = 1; x <= DORBC_02; ++x) {
      UnsignedDomain<1, DORBS_02> elt = x;
      res = res && perftest_orbit_02<DirectOrbit>(elt);
    }
  });
  cout << tdorb_02.count() << ", " << res << endl;

  res = true;
  cout << "shreier orbit: ";
  auto tsorb_02 = duration([&] {
    for (int x = 1; x <= SORBC_02; ++x) {
      UnsignedDomain<1, SORBS_02> elt = x;
      res = res && perftest_orbit_02<ShreierOrbit>(elt);
    }
  });
  cout << tsorb_02.count() << ", " << res << endl;
#endif
}
