//------------------------------------------------------------------------------
//
// HCGT, page 84. Primitive blocks for group
//
//------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <vector>
#include <queue>
#include "perms.hpp"

using std::cout;
using std::endl;
using std::make_pair;
using std::map;
using std::min;
using std::multimap;
using std::pair;
using std::swap;
using std::queue;
using std::vector;

template<typename T>
using classes_t = map<T, size_t>;

template<typename T>
using representatives_t = map<size_t, T>;

template<typename T>
using gens_t = vector<Permutation<T>>;

// good idea: consider classes over T to control domain
// then it will be: C num1, C num2, gens_t<C> gens
// with C::begin() and C::end()

template<typename T> auto
primitive_blocks(T num1, T num2, T DomainStart, T DomainEnd, gens_t<T> gens) {
  classes_t<T> classes;
  representatives_t<T> reps;
  queue<T> q;
  assert(num1 != num2);
  classes[num1] = 0;
  classes[num2] = 0;
  reps[0] = num1;
  q.push(num2);
  size_t classnum = 1; 
  
  for (T elem = DomainStart; elem != DomainEnd; ++elem) {
    if ((elem == num1) || (elem == num2)) continue;
    classes[elem] = classnum;
    reps[classnum] = elem;
    classnum += 1;
  }

  while (!q.empty()) {
    auto gamma = q.front();
    q.pop();
    for (auto gen: gens) {
      auto delta = reps[classes[gamma]];
      auto c1 = classes[gen.apply(gamma)];
      auto c2 = classes[gen.apply(delta)];
      auto kappa = reps[c1];
      auto lambda = reps[c2];
      if (kappa != lambda) {
        if (c1 > c2) {
          swap(c1, c2); // now c1 < c2
          swap(kappa, lambda);
        }
        for (auto& c: classes)
          if (c.second == c2)
            c.second = c1; // TODO: algorithm for this loop?
        reps[c1] = kappa; // because we might swap
        q.push(lambda);
      }
    }
  }
  return classes;
}

int
main()
{
  unsigned DomainStart = 1;
  unsigned DomainEnd = 7;
  unsigned num1 = 1;
  unsigned num2 = 3; // put 4 for block system #2

  Permutation<unsigned> g1 (DomainStart, DomainEnd, {{1, 2, 3, 4, 5, 6}});
  Permutation<unsigned> g2 (DomainStart, DomainEnd, {{2, 6}, {3, 5}});

  auto classes = primitive_blocks(num1, num2, DomainStart, DomainEnd, {g1, g2});

  // for (auto c: classes)
  //  cout << c.first << " " << c.second << endl;

  vector<pair<unsigned, size_t>> classelems (classes.begin(), classes.end());
  sort(classelems.begin(), classelems.end(), [](auto a, auto b) { return a.second < b.second; });
  auto current_class = classelems.front().second;
  cout << "Class: " << current_class << endl;
  for (auto elt: classelems) {
    if (elt.second != current_class) {
      cout << endl;
      current_class = elt.second;
      cout << "Class: " << current_class << endl;
    }
    cout << elt.first << " ";
  }
  cout << endl;
}
