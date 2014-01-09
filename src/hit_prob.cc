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

using Matrix = Eigen::SparseMatrix<double, Eigen::ColMajor, Index>;
using Vector = Eigen::SparseVector<double, Eigen::ColMajor, Index>;
using Solver = Eigen::SparseLU<Matrix>;
using VectorXi = Eigen::Matrix<Index, Eigen::Dynamic, 1>;
using VectorXd = Eigen::Matrix<double, Eigen::Dynamic, 1>;


HitProb::HitProb(Init i, Param p)
    : m_init(i), m_param(p) {

  m_prob.reserve(unsign(m_init.dim() * m_init.deme()));

  compute();

}


double HitProb::get(Index idx,
                    Index deme) const {

  return m_prob[unsign(m_init.deme() * idx + deme)];

}


double HitProb::get(State state, Index deme) const {

  return get(state.id(), deme);

}


HitProb HitProb::update(Param p) const {

  return HitProb(m_init, p);

}


// Because Eigen only handles
void HitProb::compute() {

  auto dim = m_init.dim();

  Matrix u(dim, dim);

  auto ndeme = m_init.deme();

  u.reserve(VectorXi::Constant(dim, (2 * ndeme * (ndeme - 1) + 1) * dim));

  vector<double> coals;
  coals.reserve(unsign(ndeme * dim));

  for (decltype(dim) i = 0; i < dim; ++i) {

    State s(m_init, i);

    Init ii(s);

    double total = 0.0;

    for (auto adj: s.neighbors()) {

      double u_val = compute_u(s, adj);
      u.insert(adj.id(), i) = u_val;

      total += u_val;

    }

    for (decltype(ndeme) deme = 0; deme < ndeme; ++deme) {

      Index choice = 2;
      auto coal = binomial(ii[deme], choice);
      total += 2.0 * coal * m_param.pop_size(deme) + ii[deme] * m_param.mut_rate(deme);

      coals.push_back(coal);

    }

    u.insert(i, i) = -total;

  }

  u.makeCompressed();

  Vector a(dim);
  a.insert(State(m_init).id()) = 1.0;

  Solver solver(u);
  if (solver.info() != Eigen::Success) {

    return;

  }

  VectorXd x = -solver.solve(a);
  if (solver.info() != Eigen::Success) {

    return;

  }

  for (decltype(dim) i = 0; i < dim; ++i) {

    for (decltype(ndeme) j = 0; j < ndeme; ++j) {

      m_prob.push_back(x(i) * 2.0 * m_param.pop_size(j) * coals[unsign(i * ndeme + j)]);

    }

  }

}


double HitProb::compute_u(State from, State to) {

  Index size = m_init.deme();

  Index i = 0;
  while (from[i] != to[i] + 1) {
    ++i;
  }

  Index j = i / size * size;
  while (from[j] != to[j] - 1) {
    ++j;
  }

  return from[i] * m_param.mig_rate(i % size, j % size);

}


};
