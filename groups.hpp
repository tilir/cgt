//------------------------------------------------------------------------------
//
//  Group theoretical algorithms
//
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

#include "perms.hpp"

#define SHOWCOSETREPS

// generators for permutation group
template<typename T>
using gens_t = vector<Permutation<T>>;

// orbit as simple set of orbit elements
template<typename T>
using compact_orbit_t = set<T>;

// orbit as map { elt => u(elt) }, so, that a^u(elt) = elt
template<typename T>
using orbit_t = map<T, Permutation<T>>;

// primitivity blocks. Index in vector is #of class, vector of elements is class
template<typename T>
using classes_t = vector<vector<T>>;

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

template<typename BlockIt>
void print_block(BlockIt gbeg, BlockIt gend) {
  cout << "[";
  for (auto git = gbeg; git != gend; ++git)
    cout << *git << ((git != prev(gend)) ? " " : "");
  cout << "]";
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
  set<Permutation<T>> Y;

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
          Y.insert(product(newgen, invert(deltait->second)));
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

// ref: HCGT, page 80
// restores U-beta from orbit element and shreier vector
template<typename T, typename RandIt, typename ShreiIt> 
auto ubeta (T num, T orbelem, RandIt gensbeg, RandIt gensend,
            ShreiIt shrbeg, ShreiIt shrend) {
  Permutation<T> res {};
  auto k = shrbeg[orbelem - T::start];
  if (k == 0)
    return make_pair(res, false);
  
  while (k != -1) {
    assert (k != 0); // we can not normally go away from orbit unwinding
    res.lmul(gensbeg[k-1]);
    auto invgen = invert(gensbeg[k-1]);
    orbelem = invgen.apply(orbelem);
    k = shrbeg[orbelem - T::start];
  }
  return make_pair(res, true);
}

// ref: HCGT, page 84
// primitive block system for given transitive group action
template<typename T, typename RandIt> 
auto primitive_blocks(T num1, T num2, RandIt gensbeg, RandIt gensend) {  
  map<T, size_t> classes;
  map<size_t, T> reps;
  queue<T> q;
  assert(num1 != num2);
  classes[num1] = 0;
  classes[num2] = 0;
  reps[0] = num1;
  q.push(num2);
  size_t classnum = 1; 

  for (auto elem = T::start; elem <= T::fin; ++elem) {
    if ((elem == num1) || (elem == num2)) continue;
    classes[elem] = classnum;
    reps[classnum] = elem;
    classnum += 1;
  }

  while (!q.empty()) {
    auto gamma = q.front();
    q.pop();
    for (auto igen = gensbeg; igen != gensend; ++igen) {
      auto delta = reps[classes[gamma]];
      auto c1 = classes[igen->apply(gamma)];
      auto c2 = classes[igen->apply(delta)];
      auto kappa = reps[c1];
      auto lambda = reps[c2];
      if (kappa != lambda) {
        if (c1 > c2) {
          swap(c1, c2); // now c1 < c2
          swap(kappa, lambda);
        }
        for (auto& c: classes)
          if (c.second == c2)
            c.second = c1;
        reps[c1] = kappa; // because we might swap
        q.push(lambda);
      }
    }
  }

  // TODO: this lengthy final part only transforms map { a => b } to 
  //       vector of vector of a for different b-s.
  //       I believe, I can do it better.

  classes_t<T> outcv;

  vector<pair<T, size_t>> classelems (classes.begin(), classes.end());
  sort(classelems.begin(), classelems.end(), 
       [](auto a, auto b) { return a.second < b.second; });
  auto current_class = classelems.front().second;
  outcv.push_back({});
  for (auto elt: classelems) {
    if (elt.second != current_class) {
      current_class = elt.second;
      outcv.push_back({});
    }
    outcv[outcv.size() - 1].push_back(elt.first);
  }

  return outcv;
}
