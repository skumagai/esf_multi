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


#include <utility>
#include <vector>

#include "init.hh"
#include "typedef.hh"
#include "util.hh"

namespace esf {


State::State()
    : m_id(0), m_deme(0), m_init(new Init) {}


State::State(State const& other)
    : m_id(other.m_id), m_deme(other.m_deme), m_init(new Init(other.m_init)) {}


State::State(State&& other)
    : m_id(other.m_id), m_deme(other.m_deme), m_init(nullptr) {

  m_init = other.m_init;

  other.m_init = nullptr;

}


State& operator=(State const& other) {

  if (this != &other) {

    delete m_init;

    m_id = other.m_id;

    m_deme = other.m_deme;

    m_init = new Init(other);

  }

  return *this;

}


State& operator=(State&& other) {

  delete m_init;

  m_id = other.m_id;

  m_deme = other.m_deme;

  m_init = other.m_init;

  other.m_init = nullptr;

  return *this;
}


State::State(::std::vector<Index> const& data)
    : m_data(data), m_id(compute_id()), m_deme(compute_deme(data)) {}


State::State(Init const& init)
    : m_data(expand_init(init)), m_id(compute_id(init)), m_deme(init.deme()) {}


State::State(Init const& init, Init id) {}


State::~State() {

  delete m_init;

}


Index State::id() const {

  return m_id;

}


void State::compute_state(Init const& init, Index id, Index deme) {

  auto idx_list = index_1_to_n(state_dim(init), id);

  m_data.reserve(deme * deme);

  for (Index i = 0; i < deme; ++i) {

    auto substate = generate_state(idx_list[i], deme, init[i]);

    m_data.insert(m_data.end(), substate.begin(), substate.end());

  }

}


void State::compute_id(Init const& init) {

  ;

}


void State::expand_init(Init const& init) {

  auto deme = init.deme();

  m_data.resize(deme * deme);

  for (Index i = 0; i < ndeme; ++i) {

    m_data[i + i * deme] = init[i];

  }

}




}
