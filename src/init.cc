// -*- mode: c++; coding: utf-8; -*-

// init.cc - Initial state

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
#include <cassert>
#include <functional>
#include <numeric>

#include "afs.hh"
#include "init.hh"
#include "state.hh"
#include "typedef.hh"
#include "util.hh"

namespace esf {

using ::std::accumulate;
using ::std::plus;
using ::std::transform;


Init::Init(value_type const& data)
    : m_data(data), m_dim(data.size()) {
  set_dim();
}


Init::Init(State const& state)
    : m_data(state.deme()), m_dim(state.deme()) {
  auto deme = state.deme();
  for (decltype(deme) d = 0; d < deme; ++d) {
    for (decltype(deme) a = 0; a < deme; ++a) {
      m_data[unsign(d)] += state[d + a * deme];
    }
  }

  set_dim();
}


Init::Init(AFS const& afs)
    : m_data(afs.deme()), m_dim(afs.deme()) {
  auto deme = afs.deme();
  for (auto allele: afs) {
    for (decltype(deme) i = 0; i < deme; ++i) {
      m_data[unsign(i)] += (allele.first)[i] * allele.second;
    }
  }

  set_dim();
}


Init::Init(Allele const& allele)
    : m_data(allele.begin(), allele.end()), m_dim(allele.deme()) {
  set_dim();
}


esf_uint_t Init::deme() const {
  return m_data.size();
}


esf_uint_t Init::operator[](esf_uint_t id) const {
  return m_data[id];
}


esf_uint_t Init::dim(esf_uint_t deme) const {
  return m_dim[deme];
}


esf_uint_t Init::dim() const {
  return accumulate(m_dim.begin(), m_dim.end(), 1U, multiplies<esf_uint_t>());
}


void Init::set_dim() {
  auto deme = m_data.size();
  transform(m_data.begin(), m_data.end(), m_dim.begin(),
            [deme](esf_uint_t i) {
              return binomial(i + deme - 1, i);
            });
}


Init::const_iterator Init::begin() const {
  return m_data.begin();
}


Init::const_iterator Init::end() const {
  return m_data.end();
}


bool operator==(Init const& a, Init const& b) {
  assert(a.deme() == b.deme());

  return a.m_data == b.m_data && a.m_dim == b.m_dim;
}


Init operator+(Init const& a, Init const& b) {
  assert(a.deme() == b.deme());

  Init::value_type data{a.m_data};
  transform(data.begin(), data.end(), b.begin(), data.begin(), plus<esf_uint_t>());
  return Init{data};
}


}
