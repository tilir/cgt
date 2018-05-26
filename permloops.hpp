//------------------------------------------------------------------------------
//
//  Permutation loops
//
//------------------------------------------------------------------------------
//
// Permutation loop like (a c d) encodes permutation
//
// a b c d e
// c b d a e
//
// i.e. fixes umentioned points and permutes a --> c --> d --> a
//
// Loops (a c d), (d a c) and (c d a) are equivalent. First is normal, i.e.
// it has smallest element first and serves as class representative
//
// Main assumptions are:
// * loop is non-empty
// * all elements in loop are unique
// * first element is smallest one
//
// use -DCHECKS to build with expensive runtime checks of this assumptions on
// creation and on any modification
//
//------------------------------------------------------------------------------

#ifndef PERMLOOPS_GUARD_
#define PERMLOOPS_GUARD_

#include "permcommon.hpp"

namespace permutations {

//------------------------------------------------------------------------------
//
// Permloop template
//
//------------------------------------------------------------------------------

template <typename T> class PermLoop {
  vector<T> loop_;

  // traits
public:
  using value_type = T;

  // ctors/dtors
public:
  // from initializer list
  PermLoop(initializer_list<T> ls);

  // from begin-end range (vector, list, etc)
  template <typename FwdIter> PermLoop(FwdIter b, FwdIter e);

  // Modifiers
public:
  // add element to permutation loop
  void add(T x);

  // inverse permutation loop:
  // (a b c) to (a c b), etc
  void inverse();

  // Getters
public:
  // get smallest element in loop
  T smallest() const { return loop_.front(); }

  // if loop is primitive e.g. like (5)
  bool is_primitive() const { return loop_.size() < 2; }

  // return true if loop contains given element
  bool contains(T x) const {
    return (find(loop_.begin(), loop_.end(), x) != loop_.end());
  }

  // apply loop to given element
  T apply(T x) const;

  // permute on table with loop
  // this is more effective then element-wise application
  template <typename RandIt> void apply(RandIt tbeg, RandIt tend) const;

  // true if loops are equal
  bool equals(const PermLoop &rhs) const { return loop_ == rhs.loop_; }

  // lexicographical less-than
  bool less(const PermLoop &rhs) const {
    size_t sz = rhs.loop_.size();
    if (loop_.size() != sz)
      return loop_.size() < sz;
    for (size_t i = 0; i != sz; ++i)
      if (loop_[i] != rhs.loop_[i])
        return loop_[i] < rhs.loop_[i];
    return false;
  }

  // size of loop
  size_t size() const { return loop_.size(); }

  // Serialization and dumps
public:
  // dump to given stream
  void dump(ostream &buffer) const;

  // return as string
  string to_string() const {
    stringstream buffer;
    dump(buffer);
    return buffer.str();
  }

  // Service functions
private:
  // CHECK postcondition consistency after modification
  void check();

  // roll to canonical: smallest element first
  void reroll() {
    auto smallest = min_element(loop_.begin(), loop_.end());
    rotate(loop_.begin(), smallest, loop_.end());
  }
};

//------------------------------------------------------------------------------
//
// Standalone operations
//
//------------------------------------------------------------------------------

template <typename T>
bool operator==(const PermLoop<T> &lhs, const PermLoop<T> &rhs) {
  return lhs.equals(rhs);
}

template <typename T>
bool operator<(const PermLoop<T> &lhs, const PermLoop<T> &rhs) {
  return lhs.less(rhs);
}

template <typename T>
bool operator!=(const PermLoop<T> &lhs, const PermLoop<T> &rhs) {
  return !operator==(lhs, rhs);
}

template <typename T> ostream &operator<<(ostream &os, const PermLoop<T> &rhs) {
  rhs.dump(os);
  return os;
}

// creates array of loops from permutation given by table
// say: [d, c, e, g, b, f, a] with type CharDomain<a, g>
// gives: [(a, d, g), (b, c, e), (f)]
template <typename T, typename RandIt = typename vector<T>::iterator,
          typename OutIt = typename vector<PermLoop<T>>::iterator>
void create_loops(RandIt tbeg, RandIt tend, OutIt lbeg);

// products all input loops over [start, fin) to minimize them
// for example:
// (a, c, f, g) (b, c, d) (a, e, d) (f, a, d, e) (b, g, f, a, e)
// simplifies to:
// (a, d, g) (b, c, e) (f)
// see TAOCP, Alg. 1.3.3B
template <typename RandIt, typename OutIt>
void simplify_loops(RandIt tbeg, RandIt tend, OutIt lbeg);

//------------------------------------------------------------------------------
//
// Ctors/dtors
//
//------------------------------------------------------------------------------

template <typename T>
PermLoop<T>::PermLoop(initializer_list<T> ls) : loop_(ls) {
  reroll();
#ifdef CHECKS
  check();
#endif
}

// from begin-end range (vector, list, etc)
template <typename T>
template <typename FwdIter>
PermLoop<T>::PermLoop(FwdIter b, FwdIter e) : loop_(b, e) {
  reroll();
#ifdef CHECKS
  check();
#endif
}

//------------------------------------------------------------------------------
//
// Modifiers
//
//------------------------------------------------------------------------------

template <typename T> void PermLoop<T>::add(T x) {
  loop_.push_back(x);
  reroll();
#ifdef CHECKS
  check();
#endif
}

template <typename T> void PermLoop<T>::inverse() {
  if (loop_.size() < 3)
    return;
  reverse(next(loop_.begin()), loop_.end());
#ifdef CHECKS
  check();
#endif
}

//------------------------------------------------------------------------------
//
// Selectors
//
//------------------------------------------------------------------------------

template <typename T> T PermLoop<T>::apply(T x) const {
  auto it = find(loop_.begin(), loop_.end(), x);
  if (it == loop_.end())
    return x;
  auto nxt = next(it);
  if (nxt == loop_.end())
    return *loop_.begin();
  return *nxt;
}

template <typename T>
template <typename RandIt>
void PermLoop<T>::apply(RandIt tbeg, RandIt tend) const {
  assert(T::start < T::fin);
  assert(loop_.front() >= T::start);
  assert(loop_.back() <= T::fin);
  assert(tend - tbeg == T::fin - T::start + 1);
  size_t nxt = loop_.front() - T::start;
  T tmp = tbeg[nxt];
  for (auto l : loop_) {
    size_t prev = nxt;
    nxt = l - T::start;
    if (l == loop_.front())
      continue;
    tbeg[prev] = tbeg[nxt];
  }
  tbeg[nxt] = tmp;
}

//------------------------------------------------------------------------------
//
// Dump and printing
//
//------------------------------------------------------------------------------

template <typename T> void PermLoop<T>::dump(ostream &os) const {
  for (const auto &t : loop_) {
    if (t == loop_.front())
      os << "(";
    os << t;
    if (t == loop_.back())
      os << ")";
    else
      os << " ";
  }
}

//------------------------------------------------------------------------------
//
// Service functions
//
//------------------------------------------------------------------------------

template <typename T> void PermLoop<T>::check() {
  // check any elements
  if (loop_.empty())
    throw logic_error("PermLoop shall be non-empty");

  // check unique elements
  set<T> uniqs(loop_.begin(), loop_.end());
  if (uniqs.size() != loop_.size())
    throw logic_error("PermLoop elements shall be unique");

  // check minimal element is first
  auto smallest = min_element(loop_.begin(), loop_.end());
  if (*smallest != loop_.front())
    throw("Unnormalized state detected");
}

//------------------------------------------------------------------------------
//
// Standalone operations
//
//------------------------------------------------------------------------------

template <typename RandIt, typename OutIt>
void create_loops(RandIt tbeg, RandIt tend, OutIt lbeg) {
  using T = typename std::decay<decltype(*tbeg)>::type;
  using OutT = typename OutIt::container_type::value_type;
  vector<bool> marked(tend - tbeg, false);

  for (auto mit = marked.begin(); mit != marked.end();
       mit = find(mit + 1, marked.end(), false)) {
    auto relt = mit - marked.begin();
    auto t = static_cast<T>(T::start + relt);
    OutT perm = {t};
    marked[relt] = true;
    for (T nxt = tbeg[relt]; nxt != t; nxt = tbeg[nxt - T::start]) {
      perm.add(nxt);
      marked[nxt - T::start] = true;
    }
    *lbeg = perm;
    lbeg++;
  }
}

template <typename RandIt, typename OutIt>
void simplify_loops(RandIt tbeg, RandIt tend, OutIt lbeg) {
  using T = typename std::decay<decltype(*tbeg)>::type::value_type;
  vector<T> table(T::fin - T::start + 1, T::start);
  iota(table.begin(), table.end(), T::start);
  for (auto loopit = make_reverse_iterator(tend);
       loopit != make_reverse_iterator(tbeg); ++loopit)
    loopit->apply(table.begin(), table.end());
  create_loops(table.begin(), table.end(), lbeg);
}

} // namespace permutations

#endif
