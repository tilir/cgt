//------------------------------------------------------------------------------
//
// Tests for permutation loops
//
//------------------------------------------------------------------------------

#include <iterator>

#include "idomain.hpp"
#include "permloops.hpp"
#include "perms.hpp"
#include "groups.hpp"

using std::back_inserter;

int test_loops() {
  cout << "Loop tests" << endl;
  using CD = CharDomain<'a', 'd'>;
  PermLoop<CD> p1{'a', 'c', 'd'};
  vector<char> vp2{'d', 'a', 'c'};
  PermLoop<CD> p2(vp2.begin(), vp2.end());
  assert(p1 == p2);
  p2.inverse();
  assert(p1 != p2);
  assert(p1.smallest() == 'a');
  assert(p1.contains('d'));
  assert(p1.apply('a') == 'c');

  vector<CD> initial{'a', 'b', 'c', 'd'};
  vector<CD> permuted{'c', 'b', 'd', 'a'};
  vector<CD> v = initial;
  p1.apply(v.begin(), v.end());
  assert(v == permuted);
  p2.apply(v.begin(), v.end());
  assert(v == initial);

  return 0;
}

int test_create_loops() {
  cout << "Create loop tests" << endl;
  using CD = CharDomain<'a', 'd'>;
  vector<PermLoop<CD>> loops1;
  PermLoop<CD> p1{'a', 'c', 'd'};
  vector<CD> permuted{'c', 'b', 'd', 'a'};
  create_loops(permuted.begin(), permuted.end(), back_inserter(loops1));
  assert(loops1.size() == 2);
  assert(loops1[0].size() == 3);
  assert(loops1[1].size() == 1);
  assert(loops1[0] == p1);

  vector<PermLoop<CD>> loops2;
  vector<CD> unpermuted{'a', 'b', 'c', 'd'};
  create_loops(unpermuted.begin(), unpermuted.end(), back_inserter(loops2));
  assert(loops2.size() == 4);
  assert(loops2[0].size() == 1);
  assert(loops2[1].size() == 1);
  assert(loops2[2].size() == 1);
  assert(loops2[3].size() == 1);

  vector<PermLoop<CD>> loops3;
  vector<CD> shifted{'d', 'a', 'b', 'c'};
  create_loops(shifted.begin(), shifted.end(), back_inserter(loops3));
  assert(loops3.size() == 1);
  assert(loops3[0].size() == 4);

  using UD = UnsignedDomain<1, 9>;
  vector<PermLoop<UD>> loops4;
  vector<UD> numbers{9, 2, 3, 1, 7, 6, 8, 5, 4};
  create_loops(numbers.begin(), numbers.end(), back_inserter(loops4));
  assert(loops4.size() == 5);
  assert(loops4[0].size() == 3);
  assert(loops4[1].size() == 1);
  assert(loops4[2].size() == 1);
  assert(loops4[3].size() == 3);
  assert(loops4[4].size() == 1);

  return 0;
}

int test_simplify_loops() {
  cout << "Simplify loops tests" << endl;
  using UD = UnsignedDomain<1, 3>;
  vector<PermLoop<UD>> in_loops0{{1, 2}, {2, 3}};
  vector<PermLoop<UD>> out_loops0;
  simplify_loops(in_loops0.begin(), in_loops0.end(), back_inserter(out_loops0));
  assert(out_loops0.size() == 1);
  PermLoop<UD> ref0{1, 3, 2};
  assert(out_loops0[0] == ref0);

  using UD4 = UnsignedDomain<1, 4>;
  vector<PermLoop<UD4>> in_loops1{{1, 3, 2}, {1, 2, 4}, {1, 4, 3, 2}};
  vector<PermLoop<UD4>> out_loops1;
  simplify_loops(in_loops1.begin(), in_loops1.end(), back_inserter(out_loops1));
  assert(out_loops1.size() == 3);

  using CG = CharDomain<'a', 'g'>;

  vector<PermLoop<CG>> in_loops{{'a', 'c', 'f', 'g'},
                                {'b', 'c', 'd'},
                                {'a', 'e', 'd'},
                                {'f', 'a', 'd', 'e'},
                                {'b', 'g', 'f', 'a', 'e'}};

  vector<PermLoop<CG>> out_loops;
  simplify_loops(in_loops.begin(), in_loops.end(), back_inserter(out_loops));

  PermLoop<CG> s1 = {'a', 'd', 'g'};
  PermLoop<CG> s2 = {'b', 'c', 'e'};
  PermLoop<CG> s3 = {'f'};

  assert(out_loops[0] == s1);
  assert(out_loops[1] == s2);
  assert(out_loops[2] == s3);

  return 0;
}

int test_simple_perms() {
  cout << "Simple perm tests" << endl;
  using UD3 = UnsignedDomain<1, 3>;

  Permutation<UD3> e;
  Permutation<UD3> g1{{1, 2}};
  Permutation<UD3> g2{{2, 3}};

  assert(e == g1.id());

  auto ne = product(e, e);
  assert(ne == e);

  auto n1 = product(e, g1);
  assert(n1 == g1);

  auto b1 = product(n1, g1);
  assert(b1 == e);

  auto n2 = product(e, g2);
  assert(n2 == g2);

  auto b2 = product(n2, g2);
  assert(b2 == e);

  auto n3 = product(n1, g2);
  Permutation<UD3> refn3{{1, 3, 2}};
  assert(n3 == refn3);

  auto b3 = product(n3, g2);
  assert(b3 == g1);

  auto n4 = product(n2, g1);
  Permutation<UD3> refn4{{1, 2, 3}};
  assert(n3 != n4);
  assert(n4 == refn4);

  auto b4 = product(n4, g1);
  assert(b4 == g2);

  auto n5 = product(n4, g2);
  Permutation<UD3> refn5{{1, 3}};
  assert(n5 == refn5);

  auto n6 = product(n3, g1);
  assert(n5 == n6);

  using UD5 = UnsignedDomain<1, 5>;

  Permutation<UD5> e6;
  Permutation<UD5> g3{{1, 2}, {3, 4, 5}};
  auto g3orig = g3;
  vector<UD5> initial{1, 2, 3, 4, 5};
  vector<UD5> permuted{2, 1, 4, 5, 3};
  auto v = initial;
  g3.apply(v.begin(), v.end());
  assert(v == permuted);
  g3.inverse();
  g3.apply(v.begin(), v.end());
  assert(v == initial);
  auto p1 = product(g3, g3orig);
  assert(p1 == e6);
  assert(product(g3orig, g3) == p1);
  g3 = g3orig;
  Permutation<UD5> g4({{1, 2, 3}, {4, 5}});
  Permutation<UD5> g5({{1, 3}, {2, 4, 5}});
  auto p34 = product(g3, g4);
  auto p45 = product(g4, g5);
  Permutation<UD5> refnp34{{1, 3, 5}};
  Permutation<UD5> refnp45{{1, 4, 2}};
  assert(p34 == refnp34);
  assert(p45 == refnp45);
  auto px = product(g3, p45);
  auto py = product(p34, g5);
  Permutation<UD5> refnpx{{2, 4, 5, 3}};
  assert(px == refnpx);
  assert(px == py);

  Permutation<UD5> u1 = {};
  Permutation<UD5> u2 = {{1, 2}};
  Permutation<UD5> u3 = {{1, 3}};
  Permutation<UD5> x1 = {{1, 2}};
  Permutation<UD5> x2 = {{2, 3}};
  assert(product(product(u3, x1), u3) == x2);
  assert(product(product(u2, x2), u3) == x2);
  assert(product(product(u3, x2), u2) == x2);

  return 0;
}

template <typename T, typename GenIT, typename OrbIt>
void do_test_simple_orbit(T elt, GenIT gbeg, GenIT gend, OrbIt refbeg,
                          OrbIt refend) {
  // 1. getting orbit
  auto simple_orbit = orbit(elt, gbeg, gend);
  for (auto rit = refbeg; rit != refend; ++rit)
    if (simple_orbit.count(*rit) == 0) {
      cout << "Orbit not equal to reference" << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Orbit is: " << endl;
      print_simple(simple_orbit.begin(), simple_orbit.end());
      cout << "Ref is: " << endl;
      print_simple(refbeg, refend);
    }

  // 2. getting orbit with ubeta and stabilizer generators
  auto[ orb, stab ] = orbit_stab(elt, gbeg, gend);
  assert(orb.size() == simple_orbit.size());
  for (auto so : simple_orbit)
    if (orb.find(so) == orb.end()) {
      cout << "Orbit-stabilizer not equal to reference" << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Orbit is: " << endl;
      print_orb(orb.begin(), orb.end());
      cout << "Ref is: " << endl;
      print_simple(refbeg, refend);
    }
  for (auto s : stab)
    if (s.apply(elt) != elt) {
      cout << "incorrect stabilizer generator " << s << " not stabilizes "
           << elt << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Stabilizer generators are: " << endl;
      print_simple(stab.begin(), stab.end());
    }
  for (auto && [ beta, u_beta ] : orb)
    if (u_beta.apply(elt) != beta) {
      cout << " incorrect ubeta " << u_beta << " for elt " << elt << " at beta "
           << beta << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Orbit is: " << endl;
      print_orb(orb.begin(), orb.end());
    }

  // 3. getting orbit with shreier and stabilizer generators
  auto[ orb_s, v ] = orbit_shreier(elt, gbeg, gend);
  assert(orb_s.size() == simple_orbit.size());
  for (auto so : simple_orbit) {
    if (orb_s.find(so) == orb_s.end()) {
      cout << "Orbit-stabilizer not equal to reference" << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Orbit is: " << endl;
      print_orb(orb.begin(), orb.end());
      cout << "Ref is: " << endl;
      print_simple(refbeg, refend);
    }
    auto ueltp = ubeta(so, gbeg, v.begin());
    assert(ueltp.second);
    auto &uelt = ueltp.first;
    if (uelt.apply(elt) != so) {
      cout << " incorrect uelt " << uelt << " for elt " << elt << " at so "
           << so << endl;
      cout << "Group generators are: " << endl;
      print_simple(gbeg, gend);
      cout << "Orbit is: " << endl;
      print_orb(orb.begin(), orb.end());
    }
  }
}

int test_simple_orbit() {
  cout << "Simple orbit tests" << endl;
  using UD5 = UnsignedDomain<1, 5>;

  // most of tests below have full orbit as reference result
  set<UD5> ref{1, 2, 3, 4, 5};

  // cyclic group
  UD5 celt = 1;
  gens_t<UD5> cgens{{{1, 5, 4, 3, 2}}};
  do_test_simple_orbit(celt, cgens.begin(), cgens.end(), ref.begin(),
                       ref.end());

  // alternating group
  UD5 aelt = 2;
  gens_t<UD5> agens{{{1, 2, 3}}, {{1, 2, 3, 4, 5}}};
  do_test_simple_orbit(aelt, agens.begin(), agens.end(), ref.begin(),
                       ref.end());

  // symmetric group
  UD5 selt = 3;
  vector<Permutation<UD5>> sgens{{{1, 2, 3, 4, 5}}, {{1, 2}}};
  do_test_simple_orbit(selt, sgens.begin(), sgens.end(), ref.begin(),
                       ref.end());

  // slightly harder case
  UD5 delt = 4;
  set<UD5> ref2{3, 4, 5};
  vector<Permutation<UD5>> dgens{{{1, 2}}, {{3, 4, 5}}};
  do_test_simple_orbit(delt, dgens.begin(), dgens.end(), ref2.begin(),
                       ref2.end());

  // isolated case
  UD5 ielt = 5;
  set<UD5> ref3{5};
  vector<Permutation<UD5>> igens{{{1, 2, 3, 4}}, {{1, 2}}};
  do_test_simple_orbit(ielt, igens.begin(), igens.end(), ref3.begin(),
                       ref3.end());

  return 0;
}

int test_primitive_blocks() {
  cout << "Primitive block tests" << endl;
  using UD6 = UnsignedDomain<1, 6>;
  vector<Permutation<UD6>> gens{{{1, 2, 3, 4, 5, 6}}, {{2, 6}, {3, 5}}};
  auto bs1 = primitive_blocks(UD6{1}, UD6{3}, gens.begin(), gens.end());

  vector<vector<UD6>> ref1 = {{1, 3, 5}, {2, 4, 6}};
  assert(bs1 == ref1);

  auto bs2 = primitive_blocks(UD6{1}, UD6{4}, gens.begin(), gens.end());

  vector<vector<UD6>> ref2 = {{1, 4}, {2, 5}, {3, 6}};
  assert(bs2 == ref2);

  return 0;
}

int test_strip() {
  cout << "Strip tests" << endl;
  using UD5 = UnsignedDomain<1, 5>;

  Permutation<UD5> e{};
  Permutation<UD5> a{{1, 2, 4, 3}};
  Permutation<UD5> b{{1, 2, 5, 4}};
  Permutation<UD5> a2 = product(a, a);
  Permutation<UD5> ab = product(a, b);
  Permutation<UD5> a3 = product(a2, a);

  gens_t<UD5> S1 = {a, b};
  orbit_t<UD5> Delta1 = {{1, e}, {2, a}, {4, a2}, {5, ab}, {3, a3}};

  Permutation<UD5> c = product(b, invert(a));
  Permutation<UD5> d = product(a2, b);
  Permutation<UD5> cd = product(c, d);

  vector<Permutation<UD5>> S2 = {c, d};
  orbit_t<UD5> Delta2 = {{2, e}, {5, c}, {3, d}, {4, cd}};

  set<Permutation<UD5>> all;
  all_elements(S1.begin(), S1.end(), std::inserter(all, all.end()));

  vector<UD5> Bref = {1, 2};
  genset_t<UD5> Sref = {S1, S2};
  delta_t<UD5> DeltaRef = {Delta1, Delta2};

  // got the same from shreier-sims
  auto[ B, S, Delta ] = shreier_sims(S1.begin(), S1.end());

  assert(B == Bref);
  assert(Delta == DeltaRef);

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
    assert(res.first == x.id() && res.second == B.end());
  }

  // test that all element NOT in group stripped to something other than end()
  gens_t<UD5> GSYM = {{{1, 2}}, {{1, 2, 3, 4, 5}}};
  set<Permutation<UD5>> allsym;
  all_elements(GSYM.begin(), GSYM.end(), std::inserter(allsym, allsym.end()));

  for (auto &&x : allsym) {
    if (all.count(x) != 0)
      continue;
    auto res = strip(x, B.begin(), B.end(), Delta.begin());
    assert(res.first != x.id() || res.second != B.end());
  }

  return 0;
}

int test_shreier_sims() {
  cout << "Schreier-Sims tests" << endl;

  using UD5 = UnsignedDomain<1, 5>;

  // worst-case test: symmetric group
  Permutation<UD5> a{{1, 2}};
  Permutation<UD5> b{{1, 2, 3, 4, 5}};
  gens_t<UD5> sgens = {a, b};
  auto[ B, S, Delta ] = shreier_sims(sgens.begin(), sgens.end());
  auto[ Bv, Sv, DeltaV, V] = shreier_sims_sv(sgens.begin(), sgens.end());

  assert(B.size() == 4);
  assert(Bv.size() == 4);

  // now we can get order of group as product of Deltas
  size_t gorder = 1;
  for (auto &&d : Delta)
    gorder *= d.size();

  // Sym(5) size is 120
  assert(gorder == 120);

  // cyclic group has size 5
  Permutation<UD5> c{{1, 2, 3, 4, 5}};
  gens_t<UD5> cgens = {c};
  auto[ B2, S2, Delta2 ] = shreier_sims(cgens.begin(), cgens.end());
  assert(B2.size() == 1);
  assert(Delta2.size() == 1);
  assert(Delta2[0].size() == 5);

  // alternating group
  gens_t<UD5> agens{{{1, 2, 3}}, {{1, 2, 3, 4, 5}}};
  auto[ BA, SA, DeltaA ] = shreier_sims(agens.begin(), agens.end());
  assert(BA.size() == 3);

  gorder = 1;
  for (auto &&d : DeltaA)
    gorder *= d.size();

  // Alt(5) size is 60
  assert(gorder == 60);

  set<Permutation<UD5>> allsym;
  all_elements(agens.begin(), agens.end(), std::inserter(allsym, allsym.end()));

  for (auto &&x : allsym) {
    auto res = strip(x, BA.begin(), BA.end(), DeltaA.begin());
    assert(res.first == x.id() && res.second == BA.end());
  }

  return 0;
}

int main() {
  try {
    test_loops();
    test_create_loops();
    test_simplify_loops();
    test_simple_perms();
    test_simple_orbit();
    test_primitive_blocks();
    test_strip();
    test_shreier_sims();
  } catch (exception &e) {
    cout << "Failed: " << e.what() << endl;
    exit(-1);
  }
  cout << "Passed" << endl;
}
