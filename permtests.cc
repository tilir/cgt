//------------------------------------------------------------------------------
//
// Tests for permutation loops
//
//------------------------------------------------------------------------------

#include <iterator>

#include "idomain.hpp"
#include "permloops.hpp"
#include "perms.hpp"

using namespace permutations;

using std::back_inserter;

#define simple_check(cond)                                                     \
  do {                                                                         \
    if (!(cond))                                                               \
      throw std::logic_error(#cond);                                           \
  } while (0)

int test_loops() {
  cout << "Loop tests" << endl;
  using CD = CharDomain<'a', 'd'>;
  PermLoop<CD> p1{'a', 'c', 'd'};
  vector<char> vp2{'d', 'a', 'c'};
  PermLoop<CD> p2(vp2.begin(), vp2.end());
  simple_check(p1 == p2);
  p2.inverse();
  simple_check(p1 != p2);
  simple_check(p1.smallest() == 'a');
  simple_check(p1.contains('d'));
  simple_check(p1.apply('a') == 'c');

  vector<CD> initial{'a', 'b', 'c', 'd'};
  vector<CD> permuted{'c', 'b', 'd', 'a'};
  vector<CD> v = initial;
  p1.apply(v.begin(), v.end());
  simple_check(v == permuted);
  p2.apply(v.begin(), v.end());
  simple_check(v == initial);

  return 0;
}

int test_create_loops() {
  cout << "Create loop tests" << endl;
  using CD = CharDomain<'a', 'd'>;
  vector<PermLoop<CD>> loops1;
  PermLoop<CD> p1{'a', 'c', 'd'};
  vector<CD> permuted{'c', 'b', 'd', 'a'};
  create_loops(permuted.begin(), permuted.end(), back_inserter(loops1));
  simple_check(loops1.size() == 2);
  simple_check(loops1[0].size() == 3);
  simple_check(loops1[1].size() == 1);
  simple_check(loops1[0] == p1);

  vector<PermLoop<CD>> loops2;
  vector<CD> unpermuted{'a', 'b', 'c', 'd'};
  create_loops(unpermuted.begin(), unpermuted.end(), back_inserter(loops2));
  simple_check(loops2.size() == 4);
  simple_check(loops2[0].size() == 1);
  simple_check(loops2[1].size() == 1);
  simple_check(loops2[2].size() == 1);
  simple_check(loops2[3].size() == 1);

  vector<PermLoop<CD>> loops3;
  vector<CD> shifted{'d', 'a', 'b', 'c'};
  create_loops(shifted.begin(), shifted.end(), back_inserter(loops3));
  simple_check(loops3.size() == 1);
  simple_check(loops3[0].size() == 4);

  using UD = UnsignedDomain<1, 9>;
  vector<PermLoop<UD>> loops4;
  vector<UD> numbers{9, 2, 3, 1, 7, 6, 8, 5, 4};
  create_loops(numbers.begin(), numbers.end(), back_inserter(loops4));
  simple_check(loops4.size() == 5);
  simple_check(loops4[0].size() == 3);
  simple_check(loops4[1].size() == 1);
  simple_check(loops4[2].size() == 1);
  simple_check(loops4[3].size() == 3);
  simple_check(loops4[4].size() == 1);

  return 0;
}

int test_simplify_loops() {
  cout << "Simplify loops tests" << endl;
  using UD = UnsignedDomain<1, 3>;
  vector<PermLoop<UD>> in_loops0{{1, 2}, {2, 3}};
  vector<PermLoop<UD>> out_loops0;
  simplify_loops(in_loops0.begin(), in_loops0.end(), back_inserter(out_loops0));
  simple_check(out_loops0.size() == 1);
  PermLoop<UD> ref0{1, 3, 2};
  simple_check(out_loops0[0] == ref0);

  using UD4 = UnsignedDomain<1, 4>;
  vector<PermLoop<UD4>> in_loops1{{1, 3, 2}, {1, 2, 4}, {1, 4, 3, 2}};
  vector<PermLoop<UD4>> out_loops1;
  simplify_loops(in_loops1.begin(), in_loops1.end(), back_inserter(out_loops1));
  simple_check(out_loops1.size() == 3);

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

  simple_check(out_loops[0] == s1);
  simple_check(out_loops[1] == s2);
  simple_check(out_loops[2] == s3);

  return 0;
}

int test_simple_perms() {
  cout << "Simple perm tests" << endl;
  using UD3 = UnsignedDomain<1, 3>;

  Permutation<UD3> e;
  Permutation<UD3> g1{{1, 2}};
  Permutation<UD3> g2{{2, 3}};

  simple_check(e == g1.id());

  auto ne = product(e, e);
  simple_check(ne == e);

  auto n1 = product(e, g1);
  simple_check(n1 == g1);

  auto b1 = product(n1, g1);
  simple_check(b1 == e);

  auto n2 = product(e, g2);
  simple_check(n2 == g2);

  auto b2 = product(n2, g2);
  simple_check(b2 == e);

  auto n3 = product(n1, g2);
  Permutation<UD3> refn3{{1, 3, 2}};
  simple_check(n3 == refn3);

  auto b3 = product(n3, g2);
  simple_check(b3 == g1);

  auto n4 = product(n2, g1);
  Permutation<UD3> refn4{{1, 2, 3}};
  simple_check(n3 != n4);
  simple_check(n4 == refn4);

  auto b4 = product(n4, g1);
  simple_check(b4 == g2);

  auto n5 = product(n4, g2);
  Permutation<UD3> refn5{{1, 3}};
  simple_check(n5 == refn5);

  auto n6 = product(n3, g1);
  simple_check(n5 == n6);

  using UD5 = UnsignedDomain<1, 5>;

  Permutation<UD5> e6;
  Permutation<UD5> g3{{1, 2}, {3, 4, 5}};
  auto g3orig = g3;
  vector<UD5> initial{1, 2, 3, 4, 5};
  vector<UD5> permuted{2, 1, 4, 5, 3};
  auto v = initial;
  g3.apply(v.begin(), v.end());
  simple_check(v == permuted);
  g3.inverse();
  g3.apply(v.begin(), v.end());
  simple_check(v == initial);
  auto p1 = product(g3, g3orig);
  simple_check(p1 == e6);
  simple_check(product(g3orig, g3) == p1);
  g3 = g3orig;
  Permutation<UD5> g4({{1, 2, 3}, {4, 5}});
  Permutation<UD5> g5({{1, 3}, {2, 4, 5}});
  auto p34 = product(g3, g4);
  auto p45 = product(g4, g5);
  Permutation<UD5> refnp34{{1, 3, 5}};
  Permutation<UD5> refnp45{{1, 4, 2}};
  simple_check(p34 == refnp34);
  simple_check(p45 == refnp45);
  auto px = product(g3, p45);
  auto py = product(p34, g5);
  Permutation<UD5> refnpx{{2, 4, 5, 3}};
  simple_check(px == refnpx);
  simple_check(px == py);

  Permutation<UD5> u1 = {};
  Permutation<UD5> u2 = {{1, 2}};
  Permutation<UD5> u3 = {{1, 3}};
  Permutation<UD5> x1 = {{1, 2}};
  Permutation<UD5> x2 = {{2, 3}};
  simple_check(product(product(u3, x1), u3) == x2);
  simple_check(product(product(u2, x2), u3) == x2);
  simple_check(product(product(u3, x2), u2) == x2);

  return 0;
}

int test_powers() {
  using UD5 = UnsignedDomain<1, 5>;
  Permutation<UD5> e;
  Permutation<UD5> g1{{1, 2, 3, 4, 5}};
  auto g2 = product(g1, g1);
  auto g3 = product(g2, g1);
  auto g4 = product(g3, g1);
  auto g5 = product(g4, g1);
  simple_check(g5 == e);
  simple_check(e == perm_pow(g1, 0));
  simple_check(g1 == perm_pow(g1, 1));
  simple_check(g2 == perm_pow(g1, 2));
  simple_check(g3 == perm_pow(g1, 3));
  simple_check(g4 == perm_pow(g1, 4));
  simple_check(e == perm_pow(g1, 5));
  simple_check(g4 == perm_pow(g1, -1));
  simple_check(g3 == perm_pow(g1, -2));
  simple_check(g2 == perm_pow(g1, -3));
  simple_check(g1 == perm_pow(g1, -4));

  return 0;
}

int main() {
  try {
    test_loops();
    test_create_loops();
    test_simplify_loops();
    test_simple_perms();
    test_powers();
  } catch (exception &e) {
    cout << "Failed: " << e.what() << endl;
    exit(-1);
  }
  cout << "Passed" << endl;
}
