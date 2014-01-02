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
#include <utility>

#include "allele.hh"


namespace esf {


namespace {

using ::std::pair;
using ::std::vector;

vector<Allele> move_genes(Allele, Index, Index);

vector<pair<Allele, vector<Index>>> combine(Allele,
                                            vector<Index>,
                                            vector<vector<Allele>>::iterator,
                                            vector<vector<Allele>>::iterator);

}


Allele::Allele(Allele const& allele, Index deme, Mode mode)
    : data(allele.data), total(allele.total), m_deme(allele.m_deme) {

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


::std::vector<::std::pair<Allele, ::std::vector<Index>>> Allele::reacheable() const {

  using ::std::vector;

  vector<Index> vals(m_deme);

  Allele base(vals);

  vector<vector<Allele>> retvals;

  for (auto i = 0; i < m_deme; ++i) {

    retvals.push_back(move_genes(base, 0, data[i]));

  }

  vector<::std::pair<Allele, vector<Index>>> data;

  for (auto a: retvals[0]) {

    auto pairs = combine(a, a.data, retvals.begin() + 1, retvals.end());

    data.insert(data.end(), pairs.begin(), pairs.end());

  }

  return data;

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


Index& Allele::operator[](Index deme) {

  return data[deme];

}


Index const& Allele::operator[](Index deme) const {

  return data[deme];

}


bool Allele::operator==(Allele const& other) const {

  return data == other.data;

}


namespace {


::std::vector<Allele> move_genes(Allele orig, Index begin, Index remaining) {

  if (begin == orig.deme() || !remaining) {

    return {orig};

  }

  using ::std::vector;

  vector<Allele> data;

  for (auto i = begin; i < orig.deme(); ++i) {

    auto alleles = move_genes(orig.add(i), i, remaining - 1);

    data.insert(data.end(), alleles.begin(), alleles.end());
  }

  return data;

}


vector<pair<Allele, vector<Index>>> combine(Allele allele,
                                            vector<Index> state,
                                            vector<vector<Allele>>::iterator begin,
                                            vector<vector<Allele>>::iterator end) {

  using ::std::make_pair;

  if (begin == end) {

    return {make_pair(allele, state)};

  }

  vector<pair<Allele, vector<Index>>> data;

  for (auto a: *begin) {

    auto s = state;

    s.insert(s.end(), a.begin(), a.end());

    auto b = allele;

    for (auto i = 0; i < allele.deme(); ++i) {

      for (auto j = a[i]; j > 0; --j) {

        b = b.add(i);

      }

    }

    auto retval = combine(b, s, begin + 1, end);

    data.insert(data.end(), retval.begin(), retval.end());

  }

  return data;

}


}


}
