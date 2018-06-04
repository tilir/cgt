//------------------------------------------------------------------------------
//
//  Group theoretical algorithms
//
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

#ifndef GROUPS_GUARD_
#define GROUPS_GUARD_

#include "orbits.hpp"

using groupgens::gens_t;

namespace groups {

// generating sets for chain of groups
template <typename T> using gensets_t = vector<gens_t<T>>;

// primitivity blocks. Index in vector is #of class, vector of elements is class
template <typename T> using classes_t = vector<vector<T>>;

// random element via product replacement algorithm
// ref: HCGT page 71
// ref: https://en.wikipedia.org/wiki/Nielsen_transformation
// r is size of generators vector
// n is amount of burn-in session
// passing r < 2 means "pick automaticaly"
// returns function, returning random group elements
// x = random_init(gbeg, gend, 10); auto t1 = x(), t2 = x(), t3 = x();
template <typename RandIt>
auto random_init(RandIt gensbeg, RandIt gensend, size_t r = 0, size_t n = 10);

// all elements from group
// will likely explode in general case, but useful for small tests
template <typename RandIt, typename OutIt>
size_t all_elements(RandIt gensbeg, RandIt gensend, OutIt results);

// ref: HCGT, page 84
// primitive block system for given transitive group action
// really returns classes_t<T>
template <typename T, typename RandIt>
auto primitive_blocks(T num1, T num2, RandIt gensbeg, RandIt gensend);

// ref: HCGT, page 89
// takes generator g, sets Base and Delta* from shreier_sims below
// returns h (possibly id). If h is not id, then g not in G.
// If h is id, then g = uk * uk-1 * ... * u1
template <typename Perm, typename BaseIt, typename DeltaIt>
auto strip(Perm g, BaseIt bstart, BaseIt bfin, DeltaIt dstart);

// ref: HCGT, page 91
// takes generators g[0] .. g[s]
// returns (B, S, Delta*) where
// B = {b[1] .. b[k]} is base set (i.e. none of g[.] fixes all of B)
// Define: Gi = Stab(G, b[1], .. b[i-1]), G1 = G, G2 fixes b1, etc...
// S is strong generating set {S1 .. Sk} if Gi == <S[i]>
// Delta* is set of orbits. Each Delta[i] is orbit of b[i] in G[i]
template <template <class> class OrbT, typename RandIt>
auto shreier_sims(RandIt gensbeg, RandIt gensend);

//------------------------------------------------------------------------------
//
// implementation
//
//------------------------------------------------------------------------------

template <typename RandIt>
auto random_init(RandIt gensbeg, RandIt gensend, size_t r, size_t n) {
  static random_device rd;
  static mt19937 g(rd());
  using PermT = typename RandIt::value_type;
  vector<PermT> x;
  PermT x0;

  if (r < 2)
    r = max<size_t>(10, distance(gensbeg, gensend));

  while (x.size() < r)
    for (auto git = gensbeg; git != gensend; ++git)
      if (x.size() < r)
        x.push_back(*git);

  auto randget = [r, x, x0]() mutable {
    size_t s = g() % r;
    size_t t = g() % (r - 1);
    if (t >= s)
      t = (t + 1) % r;
    size_t b = g() % 2;
    int e = (g() % 2) * 2 - 1;

    if (b == 0) {
      x[s].rmul(perm_pow(x[t], e));
      x0.rmul(x[s]);
    } else {
      x[s].lmul(perm_pow(x[t], e));
      x0.lmul(x[s]);
    }

    return x0;
  };

  while (n-- > 0)
    randget();

  return randget;
}

template <typename RandIt, typename OutIt>
size_t all_elements(RandIt gensbeg, RandIt gensend, OutIt results) {
  auto id = gensbeg->id();
  set<decltype(id)> next{id};
  set<decltype(id)> total;
  while (!next.empty()) {
    set<decltype(id)> tmp;
    total.insert(next.begin(), next.end());
    for (auto &&elem : next)
      for (auto igen = gensbeg; igen != gensend; ++igen)
        if (auto newelem = product(elem, *igen); total.count(newelem) == 0)
          tmp.insert(newelem);
    next.swap(tmp);
  }
  for (auto &&elt : total)
    *results++ = elt;
  return total.size();
}

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

template <typename Perm, typename BaseIt, typename DeltaIt>
auto strip(Perm g, BaseIt bstart, BaseIt bfin, DeltaIt dstart) {
  auto h = g;
  auto dit = dstart;
  for (auto bit = bstart; bit != bfin; ++bit, ++dit) {
    auto beta = h.apply(*bit);
    if (!dit->contains(beta))
      return make_pair(h, bit);
    auto u_beta = dit->ubeta(beta);
    h.rmul(invert(u_beta)); // h = h * (ui)^(-1)
  }
  return make_pair(h, bfin);
}

// shreier-sims subroutine. Try to strip generator,
// then decide do we need to extend base and/or recalculate BSGS
template <typename Gen, typename BaseIt, typename DeltaIt>
auto try_newgen(Gen newgen, BaseIt bstart, BaseIt bfin, DeltaIt dstart) {
  using T = typename BaseIt::value_type;
  bool need_recalc_orbit = false;
  bool need_extend_base = false;
  T gamma;

  auto[h, itj] = strip(newgen, bstart, bfin, dstart);
  size_t newidx = itj - bstart;
  if ((itj == bfin) && (h != h.id())) {
    need_extend_base = true;
    // looking for elt, moved by h
    auto itnonprim = find_if(h.rbegin(), h.rend(), [](const auto &elt) {
      return !elt.is_primitive();
    });
    assert(itnonprim != h.rend());
    gamma = itnonprim->smallest();

    if (find(bstart, bfin, gamma) != bfin)
      throw logic_error("Can not add duplicating gamma");
  }

  if ((itj != bfin) || need_extend_base)
    need_recalc_orbit = true;

  return make_tuple(need_recalc_orbit, need_extend_base, gamma, newidx, h);
}

// HEART OF SHREIER-SIMS
// finding way to extend/update B and S by stripping current BSGS
template <typename BaseIt, typename SetIt, typename DeltaIt>
auto extend_base(size_t curidx, BaseIt Base, BaseIt BaseEnd, SetIt Gens,
                 DeltaIt Delta) {
  for (auto &&beta : Delta[curidx]) {
    auto u_beta = Delta[curidx].ubeta(beta);
    for (auto &&x : Gens[curidx]) {
      auto ub_x = product(u_beta, x);
      auto u_bx = Delta[curidx].ubeta(x.apply(beta));
      if (ub_x != u_bx) {
        auto newgen = product(ub_x, invert(u_bx));
        auto[recalc, extend, gamma, newidx, h] =
            try_newgen(newgen, Base, BaseEnd, Delta);
        if (recalc)
          return make_tuple(true, extend, newidx, gamma, h);
      }
    }
  }

  return make_tuple(false, false, size_t(0), typename BaseIt::value_type{},
                    Gens[curidx][0].id());
}

template <template <class> class OrbT, typename RandIt>
auto shreier_sims(RandIt gensbeg, RandIt gensend) {
  using T = typename RandIt::value_type::value_type;
  gensets_t<T> S;
  vector<T> B;
  vector<OrbT<T>> DeltaStar;

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

  // TODO: we may extend algorithm to start from non-void B, but not now
  assert(B.size() == S.size());
  assert(B.size() == 1);

  DeltaStar.emplace_back(B[0], S[0].begin(), S[0].end());

  int curidx = static_cast<int>(B.size() - 1);

  while (curidx != -1) {
    auto[succ, extend, newidx, gamma, h] =
        extend_base(curidx, B.begin(), B.end(), S.begin(), DeltaStar.begin());

    if (!succ) {
      curidx -= 1;
      continue;
    }

#ifndef NDEBUG
    if ((!extend && (newidx == S.size())) || (newidx > S.size()))
      throw runtime_error("Orbit extended beyond possible");
#endif

    for (size_t l = curidx; l <= newidx; ++l) {
      if (extend && (l == newidx)) {
        assert(newidx == S.size());
        B.push_back(gamma);
        S.push_back({h});
        DeltaStar.emplace_back(B[l], S[l].begin(), S[l].end());
        continue;
      }

      // extend orbit if required
      S[l].push_back(h);
      DeltaStar[l].extend_orbit(h);
    }

    curidx = newidx;
  }

  return make_tuple(B, S, DeltaStar);
}

} // namespace groups

#endif
