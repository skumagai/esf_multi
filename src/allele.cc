// -*- mode: c++; coding: utf-8; -*-

// allele.cc - Allele

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


#include <numeric>

#include "allele.hh"


namespace esf {


Allele::Allele(Allele const& allele, Index deme, Mode mode)
    : data(allele.data), total(allele.total) {

  if (mode == Mode::ADD) {

    data[deme] += 1;
    total += 1;

  } else {

    data[deme] -= 1;
    total -= 1;

  }

}


Allele::Allele(Allele&& allele, Index deme, Mode mode)
    : total(allele.total), data(::std::move(allele.data)) {

  if (mode == Mode::ADD) {

    data[deme] += 1;
    total += 1;

  } else {

    data[deme] -= 1;
    total -= 1;

  }

}


Allele::Allele(::std::vector<Index> const& d)
    : data(d), total(::std::accumulate(d.begin(), d.end(), 0)), m_deme(d.size()) {}


Allele Allele::remove(Index deme) const {

  return Allele(*this, deme, Mode::REMOVE);

}


bool Allele::singleton() const {

  return total == 1;

}


Index Allele::size() const {

  return total;

}


Index Allele::deme() const {

  return m_deme;

}


Allele Allele::add(Index deme) const {

  return Allele(*this, deme, Mode::ADD);

}


Allele::iterator Allele::begin() {

  return data.begin();

}


Allele::const_iterator Allele::begin() const {

  return data.begin();

}


Allele::iterator Allele::end() {

  return data.end();

}


Allele::const_iterator Allele::end() const {

  return data.end();

}


bool Allele::operator<(Allele const& allele) const {

  return data < allele.data;

}


Index Allele::operator[](Index deme) const {

  return data[deme];

}


}
