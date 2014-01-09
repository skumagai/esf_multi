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
#include "util.hh"

namespace esf {


namespace {

using ::std::pair;
using ::std::size_t;
using ::std::vector;

vector<Allele> move_genes(Allele, Index, Index);

vector<ExitAllelePair> combine(Allele,
                               vector<Index>,
                               vector<vector<Allele>>::iterator,
                               vector<vector<Allele>>::iterator);

}


Allele::Allele(::std::vector<Index> const& d)
    : m_data(d),
      m_total(::std::accumulate(d.begin(), d.end(), static_cast<Index>(0))) {}


Index Allele::size() const {

  return m_total;

}


Index Allele::deme() const {

  return sign(m_data.size());

}


Allele Allele::remove(Index deme) const {

  value_type data = m_data;

  --data[unsign(deme)];

  return Allele(data);

}


Allele Allele::add(Index deme) const {

  value_type data = m_data;

  ++data[unsign(deme)];

  return Allele(data);

}


bool Allele::singleton() const {

  return m_total == 1;

}


::std::vector<ExitAllelePair> Allele::reacheable() const {

  using ::std::vector;

  auto d = deme();

  vector<Index> vals(unsign(d));

  Allele base(vals);

  vector<vector<Allele>> retvals;

  for (decltype(d) i = 0; i < d; ++i) {

    retvals.push_back(move_genes(base, 0, m_data[unsign(i)]));

  }

  vector<ExitAllelePair> data;

  for (auto a: retvals[0]) {

    auto pairs = combine(a, a.m_data, retvals.begin() + 1, retvals.end());

    data.insert(data.end(), pairs.begin(), pairs.end());

  }

  return data;

}


Allele::iterator Allele::begin() {

  return m_data.begin();

}


Allele::const_iterator Allele::begin() const {

  return m_data.begin();

}


Allele::iterator Allele::end() {

  return m_data.end();

}


Allele::const_iterator Allele::end() const {

  return m_data.end();

}


bool Allele::operator<(Allele const& allele) const {

  return m_data < allele.m_data;

}


Index& Allele::operator[](Index deme) {

  return m_data[unsign(deme)];

}


Index Allele::operator[](Index deme) const {

  return m_data[unsign(deme)];

}




bool operator==(Allele const& a, Allele const& b) {

  return a.m_data == b.m_data;

}


bool operator==(ExitAllelePair const& a, ExitAllelePair const& b) {

  return a.allele == b.allele && a.state == b.state;

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


vector<ExitAllelePair> combine(Allele allele,
                               vector<Index> state,
                               vector<vector<Allele>>::iterator begin,
                               vector<vector<Allele>>::iterator end) {

  if (begin == end) {

    return {ExitAllelePair({allele, state})};

  }

  vector<ExitAllelePair> data;

  for (auto a: *begin) {

    auto s = state;

    s.insert(s.end(), a.begin(), a.end());

    auto b = allele;

    for (decltype(allele.deme()) i = 0; i < allele.deme(); ++i) {

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
