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

//------------------------------------------------------------------------------
//
// Permloop template
//
//------------------------------------------------------------------------------

template <typename T>
class PermLoop {
  vector<T> loop_;

// ctors/dtors
public:
  // from initializer list
  PermLoop (initializer_list <T> ls);

  // from begin-end range (vector, list, etc)
  template<typename FwdIter>
  PermLoop (FwdIter b, FwdIter e);

// Modifiers
public:
  // add element to permutation loop
  void add (T x);

  // inverse permutation loop:
  // (a b c) to (a c b), etc
  void inverse();

// Getters
public:
  // get smallest element in loop
  T smallest() const { return loop_.front(); }

  // return true if loop contains given element
  bool contains (T x) const {
    return (find(loop_.begin(), loop_.end(), x) != loop_.end());
  }

  // apply loop to given element
  T apply (T x) const;

  // permute [start..fin) table with loop
  // this is more effective then element-wise application
  void apply(T start, T fin, vector<T>& table) const;

  // true if loops are equal
  bool equals (const PermLoop& rhs) const {
    return loop_ == rhs.loop_;
  }

  // size of loop
  size_t size() const { return loop_.size(); }

// Serialization and dumps
public:
  // dump to given stream
  void dump(ostream& buffer) const;

  // return as string
  string to_string() const {
    stringstream buffer;
    dump(buffer);
    return buffer.str();
  }
 
  // return as vector
  vector<T> to_vector() const { return loop_; }

// Service functions
private:
  // CHECK postcondition consistency after modification
  void check();

  // roll to canonical: smallest element first
  void reroll() {
    auto smallest = min_element (loop_.begin(), loop_.end());
    rotate(loop_.begin(), smallest, loop_.end());
  }
};

//------------------------------------------------------------------------------
//
// Standalone operations
//
//------------------------------------------------------------------------------

template <typename T>
bool operator == (const PermLoop<T>& lhs, const PermLoop<T>& rhs) {
  return lhs.equals(rhs);
}

template <typename T>
bool operator != (const PermLoop<T>& lhs, const PermLoop<T>& rhs) {
  return !operator==(lhs, rhs);
}

template <typename T>
ostream& operator<<(ostream& os, const PermLoop<T>& rhs) {
  rhs.dump(os);
  return os;
}

// creates array of loops from permutation given by table
// say: a, g, [d, c, e, g, b, f, a] 
// gives: [(a, d, g), (b, c, e), (f)]
template <typename T>
vector<PermLoop<T>> create_loops(T start, T fin, const vector<T>& table);

// products all input loops over [start, fin) to minimize them
// for example: 
// (a, c, f, g) (b, c, d) (a, e, d) (f, a, d, e) (b, g, f, a, e)
// simplifies to:
// (a, d, g) (b, c, e) (f)
// see TAOCP, Alg. 1.3.3B
template <typename T>
void simplify_loops (T start, T fin, vector<PermLoop<T>> &input);

//------------------------------------------------------------------------------
//
// Ctors/dtors
//
//------------------------------------------------------------------------------

template <typename T>
PermLoop<T>::PermLoop(initializer_list <T> ls) : loop_(ls) {
  reroll();
#ifdef CHECKS
  check();
#endif
}

// from begin-end range (vector, list, etc)
template <typename T>
template<typename FwdIter>
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

template <typename T>
void PermLoop<T>::add (T x) {
  loop_.push_back(x);
  reroll();
#ifdef CHECKS
  check();
#endif
}

template <typename T>
void PermLoop<T>::inverse() {
  if (loop_.size() < 3) return;
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

template <typename T>
T PermLoop<T>::apply (T x) const {
  auto it = find(loop_.begin(), loop_.end(), x);
  if (it == loop_.end()) return x;
  auto nxt = next(it);
  if (nxt == loop_.end()) return *loop_.begin();
  return *nxt;
}

template <typename T>
void PermLoop<T>::apply(T start, T fin, vector<T>& table) const {
  assert(start < fin);
  assert(loop_.front() >= start);
  assert(loop_.back() < fin);
  assert(table.size() == fin - start);
  size_t nxt = loop_.front() - start;
  T tmp = table[nxt];
  for (auto l : loop_) {
    size_t prev = nxt;
    nxt = l - start;
    if (l == loop_.front()) continue;
    table[prev] = table[nxt];
  }
  table[nxt] = tmp;
}

//------------------------------------------------------------------------------
//
// Dump and printing
//
//------------------------------------------------------------------------------

template <typename T>
void PermLoop<T>::dump(ostream& os) const {
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

template <typename T>
void PermLoop<T>::check() {
  // check any elements
  if (loop_.empty())
    throw logic_error("PermLoop shall be non-empty");
  
  // check unique elements
  set<T> uniqs(loop_.begin(), loop_.end());
  if (uniqs.size() != loop_.size())
    throw logic_error("PermLoop elements shall be unique");

  // check minimal element is first
  auto smallest = min_element (loop_.begin(), loop_.end());
  if (*smallest != loop_.front())
    throw ("Unnormalized state detected");
}

//------------------------------------------------------------------------------
//
// Standalone operations
//
//------------------------------------------------------------------------------

template <typename T>
vector<PermLoop<T>> create_loops(T start, T fin, const vector<T>& table) {
  assert(start < fin);
  assert(table.size() == fin - start);
  vector<PermLoop<T>> retval;
  vector<bool> marked(fin - start, false);
  vector<T> idtable(fin - start);
  iota(idtable.begin(), idtable.end(), start);

  for (auto t: idtable) {
    if (marked[t - start]) continue;
    vector<T> targ = {t};
    marked[t - start] = true;
    T nxt = table[t - start];
    while (nxt != t) {
      targ.push_back(nxt);
      marked[nxt - start] = true;
      nxt = table[nxt - start];
    }
    PermLoop<T> perm(targ.begin(), targ.end());
    retval.push_back(perm);
  }

  return retval;
}

template <typename T>
void simplify_loops (T start, T fin, vector<PermLoop<T>> &input) {
  vector<T> table(fin - start);
  iota(table.begin(), table.end(), start);
  for (auto loop : reverse(input))
    loop.apply(start, fin, table);
  input = create_loops(start, fin, table);
}

#endif