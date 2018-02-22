//------------------------------------------------------------------------------
//
// Tests for permutations and permutation loops
//
//------------------------------------------------------------------------------

#include "permcommon.hpp"
#include "permloops.hpp"
#include "perms.hpp"

int
test_loops() {
  cout << "Loop tests" << endl;
  PermLoop<char> p1 {'a', 'c', 'd'};
  vector<char> vp2 {'d', 'a', 'c'};
  PermLoop<char> p2(vp2.begin(), vp2.end());
  cout << p1 << " == " << p2 << endl;
  assert (p1 == p2);
  p2.inverse();
  cout << p1 << " != " << p2 << endl;
  assert (p1 != p2);
  assert(p1.smallest() == 'a');
  assert(p1.contains('d'));
  assert(p1.apply('a') == 'c');
  
  vector<char> initial {'a', 'b', 'c', 'd'};
  vector<char> permuted {'c', 'b', 'd', 'a'};
  vector<char> v = initial;
  p1.apply('a', 'e', v);
  assert (v == permuted);
  p2.apply('a', 'e', v);
  assert (v == initial);

  auto loops1 = create_loops('a', 'e', permuted);
  cout << "Loops from [c b d a]: ";
  for (const auto& l : loops1)
    cout << l;
  cout << endl;
  assert (loops1.size() == 2);
  assert (loops1[0].size() == 3);
  assert (loops1[1].size() == 1);
  assert (loops1[0] == p1);

  vector<PermLoop<char>> loops2 {
    {'a', 'c', 'f', 'g'}, 
    {'b', 'c', 'd'},
    {'a', 'e', 'd'},
    {'f', 'a', 'd', 'e'},
    {'b', 'g', 'f', 'a', 'e'}
  };

  cout << "Loops to simplify: ";
  for (const auto& l : loops2)
    cout << l;
  cout << endl;

  simplify_loops('a', 'h', loops2);

  cout << "Loops were simplified to: ";
  for (const auto& l : loops2)
    cout << l;
  cout << endl;
}

int
test_simple_perms() {
  cout << "Simple perm tests" << endl;

  Permutation<int> e(1, 4);
  Permutation<int> g1(1, 4, {{1, 2}});
  Permutation<int> g2(1, 4, {{2, 3}});

  cout << "e = " << e << endl;
  cout << "g1 = " << g1 << endl;
  cout << "g2 = " << g2 << endl;

  auto ne = product(e, e);
  assert (ne == e);

  auto n1 = product(e, g1);
  assert (n1 == g1);

  auto b1 = product(n1, g1);
  assert (b1 == e);

  auto n2 = product(e, g2);
  assert (n2 == g2);

  auto b2 = product(n2, g2);
  assert (b2 == e);

  auto n3 = product(n1, g2);
  cout << n1 << " * " << g2 << " = " << n3 << endl;

  auto b3 = product(n3, g2);
  cout << n3 << " * " << g2 << " = " << b3 << endl;
  assert (b3 == g1);

  auto n4 = product(n2, g1);
  cout << n2 << " * " << g1 << " = " << n4 << endl;
  assert (n3 != n4);

  auto b4 = product(n4, g1);
  cout << n4 << " * " << g1 << " = " << b4 << endl;
  assert (b4 == g2);

  auto n5 = product(n4, g2);
  cout << "e * g2 * g1 * g2 = " << n5 << endl;

  auto n6 = product(n3, g1);
  cout << "e * g1 * g2 * g1 = " << n6 << endl;

  assert (n5 == n6);

  Permutation<int> e6(1, 6);
  Permutation<int> g3(1, 6, {{1, 2},{3, 4, 5}});
  auto g3orig = g3;
  cout << "Permutation: " << g3 << endl;
  vector<int> initial {1, 2, 3, 4, 5};
  vector<int> permuted {2, 1, 4, 5, 3};
  auto v = initial;
  g3.apply(v);
  assert (v == permuted);
  g3.inverse();
  cout << "Inversed permutation: " << g3 << endl;
  g3.apply(v);
  assert (v == initial);
  auto p1 = product(g3, g3orig);
  cout << "Product of perm and inverse " << p1 << endl;
  assert (p1 == e6);
  assert (product(g3orig, g3) == p1);
  g3 = g3orig;
  Permutation<int> g4(1, 6, {{1, 2, 3},{4, 5}});
  Permutation<int> g5(1, 6, {{1, 3},{2, 4, 5}});
  auto p34 = product(g3, g4);
  auto p45 = product(g4, g5);
  cout << g3 << " * " << g4 << " = " << p34 << endl;  
  cout << g4 << " * " << g5 << " = " << p45 << endl;  
  auto px = product(g3, p45);
  auto py = product(p34, g5);
  cout << "g3 * (g4 * g5) = " << px << endl;
  cout << "(g3 * g4) * g5) = " << py << endl;
  assert (px == py);
}

int
main()
{
  try {
    test_loops();
    test_simple_perms();
  }
  catch(runtime_error &e) {
    cout << "Failed: " << e.what() << endl;
    exit(-1);
  }
  cout << "Passed" << endl;
}
