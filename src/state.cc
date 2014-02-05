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
using ::std::plus;
using ::std::ptrdiff_t;
using ::std::size_t;
using ::std::transform;
using ::std::vector;

namespace {


State::value_type state_dim(Init const&);

State::value_type generate_state(esf_uint_t, esf_uint_t, esf_uint_t);

esf_uint_t generate_id(State::value_type const&, esf_uint_t, esf_uint_t);


}


State::State(Init const& init)
    : m_init(init), m_data(expand_init()) {}


State::State(Init const& init, esf_uint_t id)
    : m_init(init), m_data(init.deme() * init.deme()) {
  m_data = compute_state(id);
}


State::State(Init const& init, value_type const& data)
    : m_init(init), m_data(data) {}


State State::move(esf_uint_t src, esf_uint_t tar) const {
  value_type data{m_data};
  if (data[src] > 0) {
    --data[src];
  }
  ++data[tar];
  return State{m_init, data};
}


vector<State> State::neighbors() const {

  auto deme = m_init.deme();

  vector<State> neighbors;

  auto size = m_data.size();

  for (esf_uint_t src = 0; src < size; ++src) {
    if (m_data[src] != 0) {

      auto d = src / deme;

      for (auto tar = deme * d; tar < deme * (d + 1); ++tar) {
        if (src != tar) {
          neighbors.push_back(this->move(src, tar));
        }
      }
    }

  }

  return neighbors;

}


esf_uint_t State::id() const {
  return compute_id();
}


esf_uint_t State::deme() const {
  return m_init.deme();
}


esf_uint_t State::operator[](esf_uint_t i) const {
  return m_data[i];
}


State::value_type State::compute_state(esf_uint_t idx) {
  auto deme = m_init.deme();

  auto idx_list = index_1_to_n(state_dim(m_init), idx);

  value_type data;
  data.reserve(deme * deme);
  for (esf_uint_t i = 0; i < deme; ++i) {

    auto substate = generate_state(idx_list[i], deme, m_init[i]);
    data.insert(data.end(), substate.begin(), substate.end());

  }

  return data;
}


esf_uint_t State::compute_id() const {
  auto deme = m_init.deme();
  value_type accum(deme);
  for (esf_uint_t i = 0; i < deme; ++i) {
    accum[i] = generate_id(m_data, deme * i, deme * i + deme);
  }
  return index_n_to_1(state_dim(m_init), accum);
}


State::value_type State::expand_init() const {
  auto deme = m_init.deme();
  value_type data(deme * deme);
  for (decltype(deme) i = 0; i < deme; ++i) {
    data[i + i * deme] = m_init[i];
  }

  return data;
}


State::const_iterator State::begin() const {
  return m_data.begin();
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


State operator+(State const& a, State const& b) {
  State::value_type data{a.begin(), a.end()};
  transform(data.begin(), data.end(), b.begin(), data.begin(), plus<esf_uint_t>());
  return State{a.m_init + b.m_init, data};
}


namespace {


State::value_type state_dim(Init const& init) {
  esf_uint_t deme = init.deme();
  State::value_type dim(deme);
  transform(init.begin(), init.end(), dim.begin(),
            [deme](esf_uint_t i) {
              return binomial(i + deme - 1, i);
            });

  return dim;
}

State::value_type generate_state(esf_uint_t idx, esf_uint_t deme, esf_uint_t gene) {

  esf_uint_t offset = 0;

  State::value_type state(deme);

  for (esf_uint_t i = 0; i < deme - 1; ++i) {
    esf_uint_t j = 0, idx_old = idx;
    while (idx <= idx_old &&
           (offset = binomial(gene - j + deme - i - 2, gene - j)) <= idx) {
      idx_old = idx;
      idx -= offset;
      ++j;
    }
    gene -= j;
    state[i] = j;
  }

  state[deme - 1] = gene;
  return state;

}


esf_uint_t generate_id(State::value_type const& state, esf_uint_t b, esf_uint_t e) {
  auto deme = e - b;

  ptrdiff_t begin = static_cast<ptrdiff_t>(b);
  ptrdiff_t end = static_cast<ptrdiff_t>(e);

  auto size = accumulate(state.begin() + begin, state.begin() + end, 0U);

  esf_uint_t idx = 0;
  for (esf_uint_t i = 0; i < deme - 1; ++i) {
    auto curr = state[b];
    ++b;
    for (esf_uint_t j = 0; j < curr; ++j) {
      idx += binomial(size - j + deme - i - 2, size - j);
    }
    size -= curr;
  }

  return idx;
}


}


ostream& operator<<(ostream& str, State const& state) {
  str << "State(";
  copy(state.begin(), state.end(), ostream_iterator<esf_uint_t>(str, ", "));
  str << "\b\b)";

  return str;
}


}
