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
#include <functional>
#include <numeric>
#include <vector>

#include "afs.hh"
#include "init.hh"
#include "state.hh"
#include "typedef.hh"
#include "util.hh"

namespace esf {


Init::Init(::std::vector<Index> const& data)
    : m_data(data) {

  set_size();

}


Init::Init(State const& state) {

  auto deme = state.deme();

  m_data.resize(unsign(deme));

  for (decltype(deme) d = 0; d < deme; ++d) {

    for (decltype(deme) a = 0; a < deme; ++a) {

      m_data[unsign(d)] += state[d + a * deme];

    }

  }

  set_size();

}


Init::Init(AFS const& afs) {

  auto deme = afs.deme();

  m_data.resize(unsign(deme));

  for (auto allele: afs) {

    for (decltype(deme) i = 0; i < deme; ++i) {

      m_data[unsign(i)] += (allele.first)[i] * allele.second;

    }

  }

  set_size();

}


Index Init::deme() const {

  return sign(m_data.size());

}


Index Init::operator[](Index id) const {

  return m_data[unsign(id)];

}


Index Init::size(Index deme) const {

  return m_size[unsign(deme)];

}


Index Init::size() const {

  using ::std::accumulate;

  return accumulate(m_size.begin(), m_size.end(), 1, multiplies<Index>());

}


void Init::set_size() {

  auto deme = sign(m_data.size());

  m_size.resize(unsign(deme));

  ::std::transform(m_data.begin(), m_data.end(), m_size.begin(),
                   [deme](Index i)
                   {
                     return binomial<Index>(i + deme - 1, i);
                   }
                   );

}


Init::iterator Init::begin() {

  return m_data.begin();

}


Init::const_iterator Init::begin() const {

  return m_data.begin();

}


Init::iterator Init::end() {

  return m_data.end();

}


Init::const_iterator Init::end() const {

  return m_data.end();

}


bool operator==(Init const& a, Init const& b) {

  return a.m_data == b.m_data && a.m_size == b.m_size;

}


}
