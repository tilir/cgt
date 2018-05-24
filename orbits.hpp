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
//
//------------------------------------------------------------------------------

#ifndef ORBITS_GUARD__
#define ORBITS_GUARD__

#include "perms.hpp"

using permutations::Permutation;

namespace orbits {

// orbit, storing all generators along with elements
template <typename T,
          typename GenIter = typename vector<Permutation<T>>::iterator>
class DirectOrbit {
  using iter_t = typename map<T, Permutation<T>>::iterator;

  T elt_;
  GenIter gbeg_, gend_;
  map<T, Permutation<T>> orb_;
  set<Permutation<T>> stab_;

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
  DirectOrbit() = default;
  DirectOrbit(T num, GenIter gensbeg, GenIter gensend);
  auto begin() { return PartialIt{orb_.begin()}; }
  auto end() { return PartialIt{orb_.end()}; }
  auto size() const { return orb_.size(); }
  bool contains(const T &x) { return orb_.find(x) != orb_.end(); }
  auto ubeta(T x) { return orb_[x]; }
  ostream &dump(ostream &os);
};

template <typename T, typename GenIter>
DirectOrbit<T, GenIter>::DirectOrbit(T num, GenIter gensbeg, GenIter gensend)
    : elt_(num), gbeg_(gensbeg), gend_(gensend) {
  map<T, Permutation<T>> next{{num, {}}};
  while (!next.empty()) {
    map<T, Permutation<T>> tmp{};
    orb_.insert(next.begin(), next.end());
    for (auto && [ elem, curgen ] : next)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto newelem = igen->apply(elem);     // beta ^ x
        auto newgen = product(curgen, *igen); // u_beta * x
        auto deltait = orb_.find(newelem);    // u_(beta ^ x)
        if (deltait == orb_.end())
          tmp.insert({newelem, newgen});
        else
          stab_.insert(product(newgen, invert(deltait->second)));
      }
    next.swap(tmp);
  }
}

template <typename T, typename GenIter>
ostream &DirectOrbit<T, GenIter>::dump(ostream &os) {
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
//
// std::array from T::start to T::fin is also possible,
// but I am worrying about stack overflows
using shreier_t = vector<int>;

// orbit, internally storing shreier vectors
template <typename T,
          typename GenIter = typename vector<Permutation<T>>::iterator>
class ShreierOrbit {
  T elt_;
  GenIter gbeg_, gend_;
  set<T> orb_;
  shreier_t v_;

public:
  ShreierOrbit() = default;
  ShreierOrbit(T num, GenIter gensbeg, GenIter gensend);
  auto begin() { return orb_.begin(); }
  auto end() { return orb_.end(); }
  bool contains(const T &x) { return orb_.find(x) != orb_.end(); }
  auto size() const { return orb_.size(); }
  auto ubeta(T orbelem);
  ostream &dump(ostream &os);
};

template <typename T, typename GenIter>
ShreierOrbit<T, GenIter>::ShreierOrbit(T num, GenIter gensbeg, GenIter gensend)
    : elt_(num), gbeg_(gensbeg), gend_(gensend) {
  vector<T> next{num};
  v_.resize(T::fin - T::start + 1);
  v_[num - T::start] = -1;

  while (!next.empty()) {
    vector<T> tmp{};
    orb_.insert(next.begin(), next.end());
    for (auto elem : next)
      for (auto igen = gensbeg; igen != gensend; ++igen)
        if (auto newelem = igen->apply(elem); orb_.count(newelem) == 0) {
          auto genidx = igen - gensbeg;
          tmp.push_back(newelem);
          v_[newelem - T::start] = genidx + 1;
        }
    next.swap(tmp);
  }
}

template <typename T, typename GenIter>
auto ShreierOrbit<T, GenIter>::ubeta(T orbelem) {
  assert(orbelem >= T::start);
  assert(orbelem <= T::fin);
  Permutation<T> res{};
  auto k = v_[orbelem - T::start];
  if (k == 0)
    return res;

  while (k != -1) {
    assert(k != 0); // we can not normally go away from orbit unwinding
    res.lmul(gbeg_[k - 1]);
    auto invgen = invert(gbeg_[k - 1]);
    orbelem = invgen.apply(orbelem);
    k = v_[orbelem - T::start];
  }
  return res;
}

template <typename T, typename GenIter>
ostream &ShreierOrbit<T, GenIter>::dump(ostream &os) {
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