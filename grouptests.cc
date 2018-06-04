//------------------------------------------------------------------------------
//
// Tests for groups
//
//------------------------------------------------------------------------------

#include "groups.hpp"
#include "idomain.hpp"

using orbits::DirectOrbit;
using orbits::ShreierOrbit;
using namespace groupgens;
using namespace groups;

#define simple_check(cond)                                                     \
  do {                                                                         \
    if (!(cond))                                                               \
      throw std::logic_error(#cond);                                           \
  } while (0)

int test_primitive_blocks() {
  cout << "Primitive block tests" << endl;
  using UD6 = UnsignedDomain<1, 6>;
  vector<Permutation<UD6>> gens{{{1, 2, 3, 4, 5, 6}}, {{2, 6}, {3, 5}}};
  auto bs1 = primitive_blocks(UD6{1}, UD6{3}, gens.begin(), gens.end());

  vector<vector<UD6>> ref1 = {{1, 3, 5}, {2, 4, 6}};
  simple_check(bs1 == ref1);

  auto bs2 = primitive_blocks(UD6{1}, UD6{4}, gens.begin(), gens.end());

  vector<vector<UD6>> ref2 = {{1, 4}, {2, 5}, {3, 6}};
  simple_check(bs2 == ref2);

  return 0;
}

template <template <class> class OrbT> int test_strip() {
  cout << "Strip tests" << endl;
  using UD5 = UnsignedDomain<1, 5>;
  using RandIt = typename gens_t<UD5>::iterator;

  Permutation<UD5> e{};
  Permutation<UD5> a{{1, 2, 4, 3}};
  Permutation<UD5> b{{1, 2, 5, 4}};
  Permutation<UD5> a2 = product(a, a);
  Permutation<UD5> ab = product(a, b);
  Permutation<UD5> a3 = product(a2, a);

  gens_t<UD5> S1 = {a, b};
  map<UD5, Permutation<UD5>> Delta1 = {
      {1, e}, {2, a}, {4, a2}, {5, ab}, {3, a3}};

  Permutation<UD5> c = product(b, invert(a));
  Permutation<UD5> d = product(a2, b);
  Permutation<UD5> cd = product(c, d);

  vector<Permutation<UD5>> S2 = {c, d};
  map<UD5, Permutation<UD5>> Delta2 = {{2, e}, {5, c}, {3, d}, {4, cd}};

  set<Permutation<UD5>> all;
  all_elements(S1.begin(), S1.end(), std::inserter(all, all.end()));

  vector<UD5> Bref = {1, 2};
  gensets_t<UD5> Sref = {S1, S2};
  vector<map<UD5, Permutation<UD5>>> DeltaRef = {Delta1, Delta2};

  // got the same from shreier-sims
  auto[B, S, Delta] = shreier_sims<OrbT>(S1.begin(), S1.end());

  simple_check(B == Bref);
  simple_check(Delta.size() == DeltaRef.size());

  for (size_t i = 0; i < Delta.size(); ++i) {
    simple_check(Delta[i].size() == DeltaRef[i].size());
    for (auto beta : Delta[i]) {
      auto u_beta = Delta[i].ubeta(beta);

      // we do not need DeltaRef[i][beta] to be exactly u_beta
      // really we need only B[i] -> beta to be correct
      simple_check(DeltaRef[i][beta].apply(B[i]) == u_beta.apply(B[i]));
    }
  }

// TODO: a lot of excessive gens in S. Cleanup by HCGT routine?
#if 0
  cout << "S" << endl;
  for (auto &&gens : S) {
    cout << " --- " << endl;
    for (auto &&g : gens)
      cout << g << endl;
  }

  cout << "Sref" << endl;
  for (auto &&gens : Sref) {
    cout << " --- " << endl;
    for (auto &&g : gens)
      cout << g << endl;
  }
#endif

  // test that all element in group stripped to end()
  for (auto &&x : all) {
    auto res = strip(x, B.begin(), B.end(), Delta.begin());
    simple_check(res.first == x.id() && res.second == B.end());
  }

  // test that all element NOT in group stripped to something other than end()
  gens_t<UD5> sgens = symmetric_gens<UD5>();
  set<Permutation<UD5>> allsym;
  all_elements(sgens.begin(), sgens.end(), std::inserter(allsym, allsym.end()));

  for (auto &&x : allsym) {
    if (all.count(x) != 0)
      continue;
    auto res = strip(x, B.begin(), B.end(), Delta.begin());
    simple_check(res.first != x.id() || res.second != B.end());
  }

  return 0;
}

template <template <class> class OrbT> int test_shreier_sims() {
  cout << "Schreier-Sims tests" << endl;

  using UD5 = UnsignedDomain<1, 5>;
  using RandIt = typename gens_t<UD5>::iterator;

  // worst-case test: symmetric group
  auto sgens = symmetric_gens<UD5>();
  auto[B, S, Delta] = shreier_sims<OrbT>(sgens.begin(), sgens.end());
  simple_check(B.size() == 4);

  // now we can get order of group as product of Deltas
  size_t gorder = 1;
  for (auto &&d : Delta)
    gorder *= d.size();

  // Sym(5) size is 120
  simple_check(gorder == 120);

  // cyclic group has size 5
  gens_t<UD5> cgens = cyclic_gens<UD5>();
  auto[B2, S2, Delta2] = shreier_sims<OrbT>(cgens.begin(), cgens.end());
  simple_check(B2.size() == 1);
  simple_check(Delta2.size() == 1);
  simple_check(Delta2[0].size() == 5);

  // alternating group
  gens_t<UD5> agens = alternating_gens<UD5>();
  auto[BA, SA, DeltaA] = shreier_sims<OrbT>(agens.begin(), agens.end());
  simple_check(BA.size() == 3);

  gorder = 1;
  for (auto &&d : DeltaA)
    gorder *= d.size();

  // Alt(5) size is 60
  simple_check(gorder == 60);

  set<Permutation<UD5>> allalt;
  all_elements(agens.begin(), agens.end(), std::inserter(allalt, allalt.end()));

  for (auto &&x : allalt) {
    auto res = strip(x, BA.begin(), BA.end(), DeltaA.begin());
    simple_check(res.first == x.id() && res.second == BA.end());
  }

  // Some special group
  Permutation<UD5> a{{1, 2, 4, 3}};
  Permutation<UD5> b{{1, 2, 5, 4}};
  gens_t<UD5> xgens = {a, b};
  auto[BX, SX, DeltaX] = shreier_sims<OrbT>(xgens.begin(), xgens.end());

  gorder = 1;
  for (auto &&d : DeltaX)
    gorder *= d.size();

  simple_check(gorder == 20);

  auto xrand = random_init(xgens.begin(), xgens.end());

  for (int x = 0; x < 10; ++x) {
    auto elt = xrand();
    auto res = strip(elt, BX.begin(), BX.end(), DeltaX.begin());
    simple_check(res.first == elt.id() && res.second == BX.end());
  }

  return 0;
}

int main() {
  try {
    test_primitive_blocks();
    test_strip<DirectOrbit>();
    test_shreier_sims<DirectOrbit>();

    test_strip<ShreierOrbit>();
    test_shreier_sims<ShreierOrbit>();
  } catch (exception &e) {
    cout << "Failed: " << e.what() << endl;
    exit(-1);
  }
  cout << "Passed" << endl;
}
