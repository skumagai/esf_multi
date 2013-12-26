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
#include <vector>

#include "typedef.hh"


namespace esf {


Init::Init(::std::vector<Index> const& data)
    : m_data(data), m_deme(data.size()) {

  set_size();

}


Init::Init(State const& state)
    : m_deme(state.deme()) {

  m_data.resize(m_deme);

  for (Index d = 0; d < m_deme; ++d) {

    for (Index a = 0; a < m_deme; ++a) {

      m_data[d] += state[d + a * m_deme];

    }

  }

  set_size();

}


Index Init::deme() const {

  return m_deme;

}


Index Init::operator[](Index id) {

  return m_data[id];

}


Index Init::size(Index deme) const {

  return m_size[deme];

}


Index Init::size() const {

  return accumulate(m_size.begin(), m_size.end(), 1, ::std::multiplies<Index>());

}


void Init::set_size() {

  m_size.resize(m_data.size());

  ::std::transform(m_data.begin(), m_data.end(), m_size.begin(),
                   [m_deme](Index i)
                   {
                     return binomial(i + m_deme - 1, i);
                   }
                   );

}


}
