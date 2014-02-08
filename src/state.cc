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
using ::std::advance;
using ::std::copy;
using ::std::ostream;
using ::std::ostream_iterator;
using ::std::plus;
using ::std::ptrdiff_t;
using ::std::size_t;
using ::std::transform;
using ::std::vector;


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

  vector<esf_uint_t> tmp(deme);
  dimensions(m_init.begin(), m_init.end(), tmp.begin());
  vector<esf_uint_t> idx_list(deme);
  index_1_to_n(tmp.begin(), tmp.end(), idx_list.begin(), idx);


  value_type data;
  data.reserve(deme * deme);
  vector<esf_uint_t> substate(deme);
  for (esf_uint_t i = 0; i < deme; ++i) {
    distribute_n(substate.begin(), substate.end(), idx_list[i], m_init[i]);
    data.insert(data.end(), substate.begin(), substate.end());
  }

  return data;
}


esf_uint_t State::compute_id() const {
  auto deme = m_init.deme();
  auto step = sign(deme);
  value_type accum(deme);
  decltype(m_data)::const_iterator iterb = m_data.begin();
  decltype(m_data)::const_iterator itere = m_data.begin();
  advance(itere, step);
  for (esf_uint_t i = 0; i < deme; ++i, advance(iterb, step), advance(itere, step)) {
    accum[i] = convert_id<esf_uint_t>(iterb, itere);
  }
  vector<esf_uint_t> tmp(deme);
  dimensions(m_init.begin(), m_init.end(), tmp.begin());
  return index_n_to_1(tmp.begin(), tmp.end(), accum.begin());
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


ostream& operator<<(ostream& str, State const& state) {
  auto end = state.end();
  --end;

  str << "State(";
  copy(state.begin(), end, ostream_iterator<esf_uint_t>(str, ", "));
  str << *end << ')';

  return str;
}


}
