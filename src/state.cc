// -*- mode: c++; coding: utf-8; -*-

// state.cc - State in structured ESF

// Copyright (C) 2013 Seiji Kumagai

// Permissi2013 is heSeiji Kumagaiby granted, free of charge, to any person obtaining a
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
#include <iostream>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include "init.hh"
#include "state.hh"
#include "typedef.hh"
#include "util.hh"

namespace esf {

using ::std::accumulate;
using ::std::copy;
using ::std::ostream;
using ::std::ostream_iterator;
using ::std::ptrdiff_t;
using ::std::size_t;
using ::std::transform;
using ::std::vector;

namespace {


vector<Index> state_dim(Init&);

vector<Index> generate_state(Index, Index, Index);

Index generate_id(vector<Index>&, Index, Index);


}


State::State(Init const& init)
    : m_init(init), m_data(expand_init()) {}


State::State(Init const& init, Index id)
    : m_init(init) {
  m_data = compute_state(id);
}


State::State(Init const& init, vector<Index> const& data)
    : m_init(init), m_data(data)) {}


vector<State> State::neighbors() const {

  auto deme = m_init.deme();

  vector<State> neighbors;

  auto size = sign(m_data.size());

  for (Index src = 0; src < size; ++src) {

    if (m_data[unsign(src)] != 0) {

      auto d = src / deme;

      for (auto tar = deme * d; tar < deme * (d + 1); ++tar) {

        if (src != tar) {

          State new_state(*this);

          new_state[src] -= 1;
          new_state[tar] += 1;

          neighbors.push_back(new_state);

        }

      }

    }

  }

  return neighbors;

}


Index State::id() const {

  return compute_id();

}


Index State::deme() const {

  return m_init.deme();

}


Index State::operator[](Index i) const {

  return m_data[unsign(i)];

}


Index& State::operator[](Index i) {

  return m_data[unsign(i)];

}


State::value_type State::compute_state(Index idx) {

  auto deme = m_init.deme();

  value_type data;

  data.reserve(unsign(deme * deme));

  auto idx_list = index_1_to_n(state_dim(m_init), idx);

  for (decltype(deme) i = 0; i < deme; ++i) {

    auto substate = generate_state(idx_list[unsign(i)], deme, m_init[i]);

    data.insert(data.end(), substate.begin(), substate.end());

  }

  return data;

}


Index State::compute_id() {

  auto deme = m_init.deme();

  vector<Index> accum(unsign(deme));

  for (decltype(deme) i = 0; i < deme; ++i) {

    accum[unsign(i)] = generate_id(m_data, deme * i, deme * i + deme);

  }

  return index_n_to_1(state_dim(m_init), accum);

}


State::value_type State::expand_init() {

  auto deme = m_init.deme();

  value_type data;

  data.resize(unsign(deme * deme));

  for (decltype(deme) i = 0; i < deme; ++i) {

    data[unsign(i + i * deme)] = m_init[i];

  }

  return data;

}


State::iterator State::begin() {

  return m_data.begin();

}


State::const_iterator State::begin() const {

  return m_data.begin();

}


State::iterator State::end() {

  return m_data.end();

}


State::const_iterator State::end() const {

  return m_data.end();

}


bool operator==(State const& a, State const& b) {

  return a.m_init == b.m_init && a.m_data == b.m_data;

}


bool operator<(State const& a, State const& b) {

  return a.m_data < b.m_data;

}


State const operator+(State const& a, State const& b) {
  vector<Index> data{a.m_data};
  transform(data.begin(), data.end(), b.m_data.begin(), data.begin(), plus<Index>());
  return State{a.m_init + b.m_init, data};
}


namespace {


vector<Index> state_dim(Init& init) {

  auto deme = init.deme();

  vector<Index> dim(unsign(deme));

  transform(init.begin(), init.end(), dim.begin(),
            [deme](Index i)
            {
              return binomial(i + deme - 1, i);
            });

  return dim;

}

vector<Index> generate_state(Index idx, Index deme, Index gene) {

  Index offset = 0;

  vector<Index> state(unsign(deme));

  for (decltype(deme) i = 0; i < deme - 1; ++i) {

    Index j = 0;

    while (idx > 0 &&
           (offset = binomial(gene - j + deme - i - 2, gene - j)) <= idx) {

      idx -= offset;

      ++j;

    }

    gene -= j;

    state[unsign(i)] = j;

  }

  state[unsign(deme - 1)] = gene;

  return state;

}


Index generate_id(vector<Index>& state, Index b, Index e) {


  auto deme = e - b;

  ptrdiff_t begin = static_cast<ptrdiff_t>(b);
  ptrdiff_t end = static_cast<ptrdiff_t>(e);

  auto size = accumulate(state.begin() + begin, state.begin() + end, 0);

  Index idx = 0;

  for (decltype(deme) i = 0; i < deme - 1; ++i) {

    auto curr = state[unsign(b)];

    ++b;

    for (decltype(curr) j = 0; j < curr; ++j) {

      idx += binomial(size - j + deme - i - 2, size - j);

    }

    size -= curr;

  }

  return idx;

}


}


ostream& operator<<(ostream& str, State const& state) {

  str << "State(";
  copy(state.begin(), state.end(), ostream_iterator<Index>(str, ", "));
  str << "\b\b)";

  return str;

}


}
