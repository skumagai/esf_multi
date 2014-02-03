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

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <utility>

#include "allele.hh"
#include "util.hh"

namespace esf {


namespace {

using ::std::accumulate;
using ::std::copy;
using ::std::make_pair;
using ::std::ostream;
using ::std::ostream_iterator;
using ::std::pair;
using ::std::plus;
using ::std::size_t;
using ::std::transform;
using ::std::vector;


ExitAlleleData synthesize_genes(ExitAlleleData const&, ExitAlleleData const&);

vector<vector<ExitAlleleData>> move_genes(Allele const&);

vector<ExitAlleleData> assign_genes(Index, Index, ExitAlleleData const&, size_t);

vector<ExitAlleleData> combine_allele_parts(ExitAlleleData const&,
                                            vector<vector<ExitAlleleData>>::const_iterator,
                                            vector<vector<ExitAlleleData>>::const_iterator);

}


Allele::Allele(::std::vector<Index> const& d)
    : m_data(d) {}


Index Allele::size() const {
  using ::std::accumulate;
  return accumulate(m_data.begin(), m_data.end(), static_cast<Index>(0));
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

  return size() == 1;

}


vector<ExitAlleleData> Allele::reacheable() const {

  auto retvals = move_genes(*this);

  auto deme = this->deme();
  ExitAlleleData archetype{vector<Index>(unsign(deme)), vector<Index>(unsign(deme * deme)), 1.0};
  return combine_allele_parts(archetype, retvals.begin(), retvals.end());
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


bool operator==(ExitAlleleData const& a, ExitAlleleData const& b) {

  return a.allele == b.allele && a.state == b.state;

}


ostream& operator<<(ostream& str, Allele const& allele) {

  str << "Allele(";
  copy(allele.begin(), allele.end(), ostream_iterator<Index>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;

}


Allele const operator+(Allele const& a, Allele const&b) {
  vector<Index> data{a.deme};
  copy(a.begin(), a.end(), data.begin());
  transform(data.begin(), data.end(), b.begin(), data.begin(), plus<Index>());
  return Allele{data};
}


namespace {

// generate all assignments of genes within an allelic class.
vector<vector<ExitAlleleData>> move_genes(Allele const& orig) {

  vector<vector<ExitAlleleData>> parts;
  auto deme = orig.deme();
  // for (auto i = 0; i < 1; ++i) {
  for (auto i = 0; i < deme; ++ i) {
    vector<Index> zeros(unsign(deme));
    vector<Index> states(unsign(deme * deme));
    vector<Index> init{unsign(deme)};
    init[i] = orig[i];
    parts.push_back(assign_genes(i, orig[i], {zeros, {init, states}, 1.0}, 0));
  }

  return parts;
}


ExitAlleleData synthesize_genes(ExitAlleleData const& a,
                                ExitAlleleData const& b) {
  return {a.allele + b.allele, a.state + b.state, a.factor * b.factor};
}


vector<ExitAlleleData> combine_allele_parts(ExitAlleleData const& archetype,
                                            vector<vector<ExitAlleleData>>::const_iterator begin,
                                            vector<vector<ExitAlleleData>>::const_iterator end) {
  if (begin == end) {
    ExitAlleleData a(archetype);
    Allele& allele = a.allele;
    vector<Index>& state = a.state;
    auto deme = allele.deme();
    for (Index i = 0; i < deme; ++i) {
      auto num = allele[i];
      for (Index j = 0; j < deme; ++j) {
        auto k = state[unsign(deme * j + i)];
        a.factor *= binomial(num, k);
        num -= k;
      }
    }
    return {a};
  }

  vector<ExitAlleleData> retval;
  for (auto data: *begin) {
    auto ret = synthesize_genes(archetype, data);
    auto allele_list = combine_allele_parts(ret, begin + 1, end);
    retval.insert(retval.end(), allele_list.begin(), allele_list.end());
  }
  return retval;
}


// generate all assignments of genes from a single deme within an
// alleleic class.
vector<ExitAlleleData> assign_genes(Index origin,
                                    Index count,
                                    ExitAlleleData const& archetype,
                                    size_t deme) {

  auto factor = archetype.factor;
  auto dim = archetype.allele.deme();

  // base case: no degree of freedom left
  if (deme == unsign(dim - 1)) {
    Allele base(archetype.allele);
    base[sign(deme)] = count;

    vector<Index> bstate(archetype.state);
    bstate[unsign(origin * dim) + deme] += count;

    return {ExitAlleleData{base, bstate, factor}};
  }

  // regular cases: factor increases by the number of available
  // aassignments.
  vector<ExitAlleleData> holder;
  for (Index i = 0; i <= count; ++i) {
    Allele base(archetype.allele);
    base[sign(deme)] = i;

    vector<Index> bstate(archetype.state);
    bstate[unsign(origin * dim) + deme] += i;

    auto retval = assign_genes(origin,
                               count - i,
                               {base, bstate, factor},
                               deme + 1);
    holder.insert(holder.end(), retval.begin(), retval.end());
  }

  return holder;
}


}


}
