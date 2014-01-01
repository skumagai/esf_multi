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

#include <algorithm>
#include <iterator>
#include <numeric>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseLU>

#include "hit_prob.hh"
#include "param.hh"
#include "state.hh"
#include "util.hh"

namespace esf {

using ::std::accumulate;

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::SparseVector<double>;
using Solver = Eigen::SparseLU<Matrix>;
using Eigen::VectorXi;
using Eigen::VectorXd;


HitProb::HitProb(Init i, Param p)
    : init(i), param(p) {

  prob.reserve(init.size() * init.deme());

  compute();

}


double HitProb::get(Index idx, Index deme) {

  if (prob.empty()) {

    compute();

  }

  return prob[init.deme() * idx + deme];

}


double HitProb::get(State state, Index deme) {

  return get(state.id(), deme);

}


void HitProb::update(Param p) {

  param = p;

  compute();

}


void HitProb::compute() {

  auto dim = init.size();

  Matrix u(dim, dim);

  auto ndeme = init.deme();

  u.reserve(VectorXi::Constant(dim, (2 * ndeme * (ndeme - 1) + 1) * dim));

  vector<Index> coals;
  coals.reserve(ndeme * dim);

  for (Index i = 0; i < dim; ++i) {

    State s(init, i);

    Init ii(s);

    double total = 0.0;

    for (auto adj: s.neighbors()) {

      double u_val = compute_u(s, adj);
      u.insert(adj.id(), i) = u_val;

      total += u_val;

    }

    for (auto deme = 0; deme < ndeme; ++deme) {

      auto coal = static_cast<double>(binomial(ii[deme], 2));

      total += 2.0 * coal * param.pop_size(deme) + ii[deme] * param.mut_rate(deme);

      coals.push_back(coal);

    }

    u.insert(i, i) = -total;

  }

  u.makeCompressed();

  Vector a(dim);
  a.insert(State(init).id()) = 1.0;

  Solver solver(u);
  if (solver.info() != Eigen::Success) {

    return;

  }

  VectorXd x = -solver.solve(a);
  if (solver.info() != Eigen::Success) {

    return;

  }

  for (Index i = 0; i < dim; ++i) {

    for (Index j = 0; j < ndeme; ++j) {

      prob.push_back(x(i) * 2.0 * param.pop_size(j) * coals[i * ndeme + j]);

    }

  }

}


double HitProb::compute_u(State from, State to) {

  Index size = init.deme();

  Index i = 0;
  while (from[i] != to[i] + 1) {
    ++i;
  }

  Index j = i / size * size;
  while (from[j] != to[j] - 1) {
    ++j;
  }

  return from[i] * param.mig_rate(i % size, j % size);

}


};
