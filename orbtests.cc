//------------------------------------------------------------------------------
//
// Tests for orbits
//
//------------------------------------------------------------------------------

#include "idomain.hpp"
#include "orbits.hpp"

using namespace orbits;

#define orbit_check(cond, orb)                                                 \
  do {                                                                         \
    if (!(cond)) {                                                             \
      cerr << orb << endl;                                                     \
      throw std::logic_error(#cond);                                           \
    }                                                                          \
  } while (0)

template <typename T> using gens_t = vector<Permutation<T>>;

template <template <class, class> class Orb, typename T, typename GenIT,
          typename OrbIt>
void do_test_simple_orbit(T elt, GenIT gbeg, GenIT gend, OrbIt refbeg,
                          OrbIt refend) {
  // getting orbit and checking elements
  Orb<T, GenIT> orbit(elt, gbeg, gend);
  for (auto rit = refbeg; rit != refend; ++rit)
    orbit_check(orbit.contains(*rit), orbit);

    // TODO: test stabilizer
#if 0
  for (auto s : stab)
    orbit_check(s.apply(elt) == elt, orbit);
#endif

  for (auto &&beta : orbit) {
    auto u_beta = orbit.ubeta(beta);
    orbit_check(u_beta.apply(elt) == beta, orbit);
  }
}

template <template <class, class> class Orb> int test_simple_orbit() {
  cout << "Simple orbit tests" << endl;
  using UD5 = UnsignedDomain<1, 5>;

  // most of tests below have full orbit as reference result
  set<UD5> ref{1, 2, 3, 4, 5};

  // cyclic group
  UD5 celt = 1;
  gens_t<UD5> cgens{{{1, 5, 4, 3, 2}}};
  do_test_simple_orbit<Orb>(celt, cgens.begin(), cgens.end(), ref.begin(),
                            ref.end());

  // alternating group
  UD5 aelt = 2;
  gens_t<UD5> agens{{{1, 2, 3}}, {{1, 2, 3, 4, 5}}};
  do_test_simple_orbit<Orb>(aelt, agens.begin(), agens.end(), ref.begin(),
                            ref.end());

  // symmetric group
  UD5 selt = 3;
  vector<Permutation<UD5>> sgens{{{1, 2, 3, 4, 5}}, {{1, 2}}};
  do_test_simple_orbit<Orb>(selt, sgens.begin(), sgens.end(), ref.begin(),
                            ref.end());

  // slightly harder case
  UD5 delt = 4;
  set<UD5> ref2{3, 4, 5};
  vector<Permutation<UD5>> dgens{{{1, 2}}, {{3, 4, 5}}};
  do_test_simple_orbit<Orb>(delt, dgens.begin(), dgens.end(), ref2.begin(),
                            ref2.end());

  // isolated case
  UD5 ielt = 5;
  set<UD5> ref3{5};
  vector<Permutation<UD5>> igens{{{1, 2, 3, 4}}, {{1, 2}}};
  do_test_simple_orbit<Orb>(ielt, igens.begin(), igens.end(), ref3.begin(),
                            ref3.end());

  return 0;
}

int main() {
  try {
    test_simple_orbit<DirectOrbit>();
    test_simple_orbit<ShreierOrbit>();
  } catch (exception &e) {
    cerr << "Failed: " << e.what() << endl;
    exit(-1);
  }
  cout << "Passed" << endl;
}
