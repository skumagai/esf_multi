// -*- mode: c++; coding: utf-8; -*-

// hit_prob.cc - Computing hitting probabilities.

// Copyright (C) 2013 Seiji Kumagai

// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice (including the next
// paragraph) shall be included in all copies or substantial portions of the
// Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.


#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseLU>

#include "hit_prob.hh"
#include "state_space.hh"
#include "util.hh"

namespace esf {


using Matrix = Eigen::SparseMatrix<Value>;
using Vector = Eigen::SparseVector<Value>;
using Solver = Eigen::SparseLU<Matrix>;
using Eigen::VectorXi;
using Eigen::VectorXd;

using StateSpace::index_to_state;
using StateSpace::init_to_state;
using StateSpace::neighbors;
using StateSpace::state_to_index;
using StateSpace::state_to_init;
using StateSpace::total_state;


Value Params::mig_rate(Index i, Index j) {

  return mig[i + j * pop.size()];

}


Value Params::mut_rate(Index i) {

  return mut[i];

}


Value Params::pop_size(Index i) {

  return pop[i];

}


HitProb::HitProb(Init i, Params p)
    : init(i), params(p) {

  auto n = total_state(init);

  prob.reserve(n);

  demeprob.reserve(n * init.size());

  compute();

}


Value HitProb::get(Index idx, Index deme) {

  if (prob.empty()) {

    compute();

  }

  return prob[idx] * demeprob[init.size() * idx + deme];

}


Value HitProb::get(State state, Index deme) {

  return get(state_to_index(init, state), deme);

}


void HitProb::update(Params p) {

  params = p;

  compute();

}


void HitProb::compute() {

  auto dim = total_state(init);

  Matrix u(dim, dim);

  auto ndeme = init.size();

  u.reserve(VectorXi::Constant(dim, (2 * ndeme * (ndeme - 1) + 1) * dim));

  Value u_val, total;

  auto itrbase = demeprob.begin();

  State s;
  Init ii;
  for (Index i = 0; i < dim; ++i) {

    s = index_to_state(init, i);

    set_deme_specific_rate(state_to_init(s), dim, ndeme);

    total = 0.0;

    for (auto adj: neighbors(init, i)) {

      u_val = compute_u(s, adj);
      u.insert(state_to_index(init, adj), i) = u_val;

      total += u_val;

    }

    u.insert(i, i) = -total;

  }

  u.makeCompressed();

  Vector a(dim);
  a.insert(state_to_index(init, init_to_state(init))) = 1.0;

  Solver solver(u);
  if (solver.info() != Eigen::Success) {
    return;
  }

  auto x = -solver.solve(a).col(0);
  if (solver.info() != Eigen::Success) {
    return;
  }

  for (Index i = 0; i < dim; ++i) {
    prob.push_back(x(i));
  }

}


Value HitProb::compute_u(State from, State to) {

  Index size = init.size();

  Index i = 0;
  while (from[i] != to[i] + 1) {
    ++i;
  }

  Index j = i / size * size;
  while (from[j] != to[j] - 1) {
    ++j;
  }

  return from[i] * params.mig_rate(i % size, j % size);

}


void HitProb::set_deme_specific_rate(Init genes, Index dim, Index ndeme) {

  for (Index i = 0; i < dim; ++i) {

    for (Index j = 0; j < ndeme; ++j) {

      demeprob[i * ndeme + j] = 2.0 * params.pop_size(j) * binomial(genes[i], 2);

    }

  }

}


};
