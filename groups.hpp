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
template <typename T> using gens_t = vector<Permutation<T>>;

// orbit as simple set of orbit elements
template <typename T> using compact_orbit_t = set<T>;

// orbit as map { elt => u(elt) }, so, that a^u(elt) = elt
template <typename T> using orbit_t = map<T, Permutation<T>>;

// primitivity blocks. Index in vector is #of class, vector of elements is class
template <typename T> using classes_t = vector<vector<T>>;

// generating set
template <typename T> using genset_t = vector<gens_t<T>>;

// orbits of BSGS
template <typename T> using delta_t = vector<orbit_t<T>>;

// shreier vector definition:
// v[num] = -1
// v[elt not in num orbit] = 0
// v[newelem] = i, where i is (#newgen + 1)
using shreier_t = vector<int>;

template <typename GenIt> void print_simple(GenIt gbeg, GenIt gend) {
  for (auto git = gbeg; git != gend; ++git)
    cout << *git << endl;
}

template <typename OrbIt> void print_orb(OrbIt obeg, OrbIt oend) {
  for (auto oit = obeg; oit != oend; ++oit)
    cout << oit->first << ": " << oit->second << endl;
  cout << endl;
}

template <typename BlockIt> void print_block(BlockIt gbeg, BlockIt gend) {
  cout << "[";
  for (auto git = gbeg; git != gend; ++git)
    cout << *git << ((git != prev(gend)) ? " " : "");
  cout << "]";
}

// all elements from group
// will likely explode in egneral case, but useful for small tests
template <typename RandIt, typename OutIt>
size_t all_elements(RandIt gensbeg, RandIt gensend, OutIt results) {
  auto id = gensbeg->id();
  set<decltype(id)> next{id};
  set<decltype(id)> total;
  while (!next.empty()) {
    set<decltype(id)> tmp;
    total.insert(next.begin(), next.end());
    for (auto &&elem : next)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto newelem = product(elem, *igen);
        if (total.count(newelem) == 0)
          tmp.insert(newelem);
      }
    next.swap(tmp);
  }
  for (auto &&elt : total)
    *results++ = elt;
  return total.size();
}

// ref: HCGT, page 78
// calculates orbit for 'a' in compact form (i.e. as a set)
template <typename T, typename RandIt>
auto orbit(T num, RandIt gensbeg, RandIt gensend) {
  compact_orbit_t<T> Delta;
  vector<T> DeltaNext{num};

  while (!DeltaNext.empty()) {
    vector<T> tmp{};
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
template <typename T, typename RandIt>
auto orbit_stab(T num, RandIt gensbeg, RandIt gensend) {
  orbit_t<T> Delta, DeltaNext{{num, {}}};
  set<Permutation<T>> Y;

  while (!DeltaNext.empty()) {
    orbit_t<T> tmp{};
    Delta.insert(DeltaNext.begin(), DeltaNext.end());
    for (auto && [ elem, curgen ] : DeltaNext)
      for (auto igen = gensbeg; igen != gensend; ++igen) {
        auto newelem = igen->apply(elem);     // beta ^ x
        auto newgen = product(curgen, *igen); // u_beta * x
        auto deltait = Delta.find(newelem);   // u_(beta ^ x)
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
template <typename T, typename RandIt>
auto orbit_shreier(T num, RandIt gensbeg, RandIt gensend) {
  compact_orbit_t<T> Delta;
  vector<T> DeltaNext{num};
  shreier_t v(T::fin - T::start + 1);
  v[num - T::start] = -1;

  while (!DeltaNext.empty()) {
    vector<T> tmp{};
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
template <typename T, typename RandIt, typename ShreiIt>
auto ubeta(T orbelem, RandIt gensbeg, RandIt gensend, ShreiIt shrbeg,
           ShreiIt shrend) {
  Permutation<T> res{};
  auto k = shrbeg[orbelem - T::start];
  if (k == 0)
    return make_pair(res, false);

  while (k != -1) {
    assert(k != 0); // we can not normally go away from orbit unwinding
    res.lmul(gensbeg[k - 1]);
    auto invgen = invert(gensbeg[k - 1]);
    orbelem = invgen.apply(orbelem);
    k = shrbeg[orbelem - T::start];
  }
  return make_pair(res, true);
}

// ref: HCGT, page 84
// primitive block system for given transitive group action
template <typename T, typename RandIt>
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
    if ((elem == num1) || (elem == num2))
      continue;
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
        for (auto &c : classes)
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

  vector<pair<T, size_t>> classelems(classes.begin(), classes.end());
  sort(classelems.begin(), classelems.end(),
       [](auto a, auto b) { return a.second < b.second; });
  auto current_class = classelems.front().second;
  outcv.push_back({});
  for (auto elt : classelems) {
    if (elt.second != current_class) {
      current_class = elt.second;
      outcv.push_back({});
    }
    outcv[outcv.size() - 1].push_back(elt.first);
  }

  return outcv;
}

// ref: HCGT, page 89
// takes generator g, sets Base and Delta* from shreier_sims below
// returns h (possibly id). If h is not id, then g not in G.
// If h is id, then g = uk * uk-1 * ... * u1
template <typename Perm, typename BaseIt, typename DeltaIt>
auto strip(Perm g, BaseIt bstart, BaseIt bfin, DeltaIt dstart) {
  auto h = g;
  auto dit = dstart;
  for (auto bit = bstart; bit != bfin; ++bit, ++dit) {
    auto beta = h.apply(*bit);
    auto di = dit->find(beta);
    if (di == dit->end())
      return make_pair(h, bit);
    assert(di->second.apply(*bit) == beta);
    h.rmul(invert(di->second)); // h = h * (ui)^(-1)
  }
  return make_pair(h, bfin);
}

// ref: HCGT, page 91
// takes generators g[0] .. g[s]
// returns (B, Delta*) where
// B = {b[1] .. b[k]} is base set (i.e. none of g[.] fixes all of B)
// Define: Gi = Stab(G, b[1], .. b[i-1]), G1 = G, G2 fixes b1, etc...
// S is strong generating set {S1 .. Sk} if Si == Gi
// Delta* is set of orbits Delta[i] as map { elt => gen }
//        each Delta[i] is orbit of b[i] in <S[i]>
template <typename RandIt> auto shreier_sims(RandIt gensbeg, RandIt gensend) {
  using T = typename RandIt::value_type::value_type;
  genset_t<T> S;
  vector<T> B;
  delta_t<T> DeltaStar;

  // in terms of book, S1 = S
  S.emplace_back(gensbeg, gensend);

  // looking for first base element. It shall not be fixed by all generators
  for (auto b = T::start; b <= T::fin; ++b)
    if (find_if(gensbeg, gensend, [b](const auto &elt) {
          return elt.apply(b) == b;
        }) == gensend) {
      B.push_back(b);
      break;
    }

  if (B.empty())
    throw logic_error("Domain for Schreier-Sims shall have at least one element"
                      " not fixed by all generators");

  size_t k = B.size();
  assert(k == S.size());

  // TODO: we may extend algorith to start from non-void B, but not now
  assert(k == 1);

  auto[ Delta1, G1 ] = orbit_stab(B[k - 1], S[k - 1].begin(), S[k - 1].end());
  DeltaStar.push_back(Delta1);

  size_t i = k;

  while (i >= 1) {
    bool globalcont = false;

    for (auto && [ beta, ubeta ] : DeltaStar[i - 1]) {
      assert(ubeta.apply(B[i - 1]) == beta);
      for (auto &&x : S[i - 1]) {
        auto ub_x = product(ubeta, x);
        auto u_bx = DeltaStar[i - 1][x.apply(beta)];
        if (ub_x != u_bx) {
          bool y = true;
          auto newgen = product(ub_x, invert(u_bx));
          auto[ h, itj ] = strip(newgen, B.begin(), B.end(), DeltaStar.begin());
          size_t j = itj - B.begin() + 1;
          if ((itj != B.end()) || (h != h.id()))
            y = false;
          if ((itj == B.end()) && (h != h.id())) {
            // looking for elt, moved by h
            auto itnonprim =
                find_if(h.rbegin(), h.rend(),
                        [](const auto &elt) { return !elt.is_primitive(); });
            assert(itnonprim != h.rend());
            auto gamma = itnonprim->smallest();

            if (find(B.begin(), B.end(), gamma) != B.end())
              throw logic_error("Can not add duplicating gamma");

            B.push_back(gamma);
            S.push_back({});
            DeltaStar.push_back({});
            k = k + 1;
          }

          // check do we need to update Delta
          if (y == false) {
            for (size_t l = i; l <= j; ++l) {
              assert(l <= S.size());
              S[l - 1].push_back(h);
              auto[ Deltal, Gl ] =
                  orbit_stab(B[l - 1], S[l - 1].begin(), S[l - 1].end());

              DeltaStar[l - 1] = Deltal;
            }
            i = j;
            globalcont = true;
            break;
          }
        }
      }
      if (globalcont)
        break;
    }

    if (globalcont)
      continue;

    i -= 1;
  }

  return make_pair(B, DeltaStar);
}

