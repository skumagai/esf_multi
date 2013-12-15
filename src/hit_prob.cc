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
  prob.reserve(init.size() * init.size());
  compute();

}


Value HitProb::get(Index idx) {

  if (prob.empty()) {

    compute();

  }

  return prob[idx];

}


Value HitProb::get(IndexList idx) {

  return get(state_to_index(init, idx));

}


void HitProb::update(Params p) {

  params = p;

  compute();

}


void HitProb::compute() {

  auto dim = total_state(init);

  Matrix u(dim, dim), v(dim, dim);

  auto ndeme = init.size();

  u.reserve(VectorXi::Constant(dim, (2 * ndeme * (ndeme - 1) + 1) * dim));
  v.reserve(VectorXi::Constant(dim, 1));

  Value u_val, v_val, total;

  State s;
  for (Index i = 0; i < dim; ++i) {

    s = index_to_state(init, i);

    v_val = compute_v(s);
    v.insert(i, i) = v_val;

    total = v_val;

    for (auto adj: neighbors(init, i)) {

      u_val = compute_u(s, adj);
      u.insert(state_to_index(init, adj), i) = u_val;

      total += u_val;

    }

    u.insert(i, i) = -total;

  }

  u.makeCompressed();
  v.makeCompressed();

  Vector a(dim);
  a.insert(state_to_index(init, init_to_state(init))) = 1.0;

  Solver solver(u);
  if (solver.info() != Eigen::Success) {
    return;
  }

  auto x = (-v * solver.solve(a)).col(0);
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


Value HitProb::compute_v(State state) {

  auto genes = state_to_init(state);

  Value data = 0.0;
  for (Index i = 0; i < genes.size(); ++i) {

    data += 2.0 * params.pop_size(i) * binomial(genes[i], 2);
    data += genes[i] * params.mut_rate(i);

  }

  return data;

}


};
