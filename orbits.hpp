//------------------------------------------------------------------------------
//
//  Orbits
//
//------------------------------------------------------------------------------
//
// Orbit of element a over group G is all (distinct) a^g
//
// orbit concept requires:
// 1. begin and end for elements
// 2. ubeta: exact group element for b from orbit, such as a^ubeta = b
// 3. contains: if orbit contains element
// 4. size: size of orbit
// 5. dump: pretty-print orbit
// 6. extend_orbit: extends orbit (takes extended generators set)
//
//------------------------------------------------------------------------------

#ifndef ORBITS_GUARD__
#define ORBITS_GUARD__

#include "groupgens.hpp"
#include "perms.hpp"

using permutations::Permutation;

namespace orbits {

// orbit, storing all generators along with elements
template <typename T> class DirectOrbit {
  using iter_t = typename map<T, Permutation<T>>::iterator;

  T elt_;
  map<T, Permutation<T>> orb_;
  set<Permutation<T>> gens_;

  struct PartialIt {
    iter_t cur_;
    auto operator*() { return cur_->first; }
    auto operator-> () { return &cur_->first; }
    auto operator++() {
      ++cur_;
      return *this;
    }
    auto operator++(int) {
      PartialIt tmp(*this);
      ++cur_;
      return tmp;
    }
    bool operator==(const PartialIt &rhs) { return cur_ == rhs.cur_; }
    bool operator!=(const PartialIt &rhs) { return cur_ != rhs.cur_; }
  };

public:
  using type = T;

  template <typename GenIter>
  DirectOrbit(T num, GenIter gensbeg, GenIter gensend) : elt_(num) {
    orb_.insert({num, {}});
    gens_.insert(gensbeg, gensend);
    extend_orbit();
  }

  void extend_orbit();
  void extend_orbit(const Permutation<T> &newgen) {
    auto oldsize = gens_.size();
    gens_.insert(newgen);
    if (oldsize != gens_.size())
      extend_orbit();
  }
  auto begin() { return PartialIt{orb_.begin()}; }
  auto end() { return PartialIt{orb_.end()}; }
  auto size() const { return orb_.size(); }
  bool contains(const T &x) const { return orb_.find(x) != orb_.end(); }
  auto ubeta(T x) { return orb_[x]; }
  ostream &dump(ostream &os);
};

template <typename T> void DirectOrbit<T>::extend_orbit() {
  map<T, Permutation<T>> next = orb_;
  while (!next.empty()) {
    map<T, Permutation<T>> tmp{};
    for (auto &&gen : gens_)
      for (auto && [ elem, curgen ] : next)
        if (auto newelem = gen.apply(elem); orb_.find(newelem) == orb_.end())
          tmp.emplace(newelem, product(curgen, gen));
    next.swap(tmp);
    orb_.insert(next.begin(), next.end());
  }
}

template <typename T> ostream &DirectOrbit<T>::dump(ostream &os) {
  os << "[ ";
  for (auto &&oit : orb_)
    os << oit.first << ": " << oit.second << " ";
  os << "]";
  return os;
}

template <typename T> ostream &operator<<(ostream &os, DirectOrbit<T> d) {
  return d.dump(os);
}

// shreier vector definition:
// v[num] = -1
// v[elt not in num orbit] = 0
// v[newelem] = i, where i is (#newgen + 1)

// orbit, internally storing shreier vectors
template <typename T> class ShreierOrbit {
  T elt_;
  set<T> orb_;
  vector<int> v_;
  vector<Permutation<T>> gens_;

  // also storing inverse generators to speed up critical ubeta section
  vector<Permutation<T>> invgens_;

public:
  using type = T;

  template <typename GenIter>
  ShreierOrbit(T num, GenIter gensbeg, GenIter gensend) : elt_(num) {
    gens_.assign(gensbeg, gensend);
    transform(gens_.begin(), gens_.end(), back_inserter(invgens_),
              [](auto x) { return x.inverse(); });
    orb_.insert(num);
    v_.resize(T::fin - T::start + 1);
    v_[num - T::start] = -1;
    extend_orbit();
  }

  void extend_orbit();
  void extend_orbit(const Permutation<T> &newgen) {
    auto oldsize = gens_.size();
    if (find(gens_.begin(), gens_.end(), newgen) == gens_.end()) {
      gens_.push_back(newgen);
      invgens_.push_back(invert(newgen));
      extend_orbit();
    }
  }
  auto begin() { return orb_.begin(); }
  auto end() { return orb_.end(); }
  bool contains(const T &x) { return orb_.find(x) != orb_.end(); }
  auto size() const { return orb_.size(); }
  auto ubeta(T orbelem);
  ostream &dump(ostream &os);
};

template <typename T> void ShreierOrbit<T>::extend_orbit() {
  set<T> next = orb_;
  while (!next.empty()) {
    set<T> tmp;
    for (auto elem : next)
      for (auto igen = gens_.begin(); igen != gens_.end(); ++igen)
        if (auto newelem = igen->apply(elem); orb_.count(newelem) == 0) {
          auto genidx = igen - gens_.begin();
          tmp.insert(newelem);
          v_[newelem - T::start] = genidx + 1;
        }
    next.swap(tmp);
    orb_.insert(next.begin(), next.end());
  }
}

template <typename T> auto ShreierOrbit<T>::ubeta(T orbelem) {
  assert(orbelem >= T::start);
  assert(orbelem <= T::fin);
  Permutation<T> res{};
  auto k = v_[orbelem - T::start];
  if (k == 0)
    return res;

  while (k != -1) {
    assert(k != 0); // we can not normally go away from orbit unwinding
    const auto &gen = gens_[k - 1];
    res.lmul(gen);
    const auto &invgen = invgens_[k - 1];
    auto neworbelem = invgen.apply(orbelem);
    assert(orbelem != neworbelem); // we can not loop forever
    orbelem = neworbelem;
    k = v_[orbelem - T::start];
  }
  return res;
}

template <typename T> ostream &ShreierOrbit<T>::dump(ostream &os) {
  os << "[ ";
  for (auto &&oit : orb_)
    os << oit << ": " << ubeta(oit) << " ";
  os << "]";
  return os;
}

template <typename T> ostream &operator<<(ostream &os, ShreierOrbit<T> d) {
  return d.dump(os);
}

} // namespace orbits

#endif
