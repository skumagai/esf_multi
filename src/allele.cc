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

using ::std::accumulate;
using ::std::copy;
using ::std::make_pair;
using ::std::ostream;
using ::std::ostream_iterator;
using ::std::pair;
using ::std::plus;
using ::std::size_t;
using ::std::transform;


namespace {

vector<vector<ExitAlleleData>> move_genes(Allele const&);

vector<ExitAlleleData> assign_genes(esf_uint_t,
                                    esf_uint_t,
                                    esf_uint_t,
                                    Init const&,
                                    Allele const,
                                    State::value_type&);

vector<ExitAlleleData> combine_allele_parts(ExitAlleleData const&,
                                            vector<vector<ExitAlleleData>>::const_iterator,
                                            vector<vector<ExitAlleleData>>::const_iterator);

ExitAlleleData synthesize_genes(ExitAlleleData const&, ExitAlleleData const&);

}


Allele::Allele(value_type const& d)
    : m_data(d) {}


esf_uint_t Allele::size() const {
  return accumulate(m_data.begin(), m_data.end(), 0U);
}


esf_uint_t Allele::deme() const {
  return m_data.size();
}

Allele Allele::remove(esf_uint_t deme, esf_uint_t num) const {
  value_type data = m_data;
  if (data[deme] >= num) {
    data[deme] -= num;
  } else {
    data[deme] = 0U;
  }
  return Allele{data};
}


Allele Allele::add(esf_uint_t deme, esf_uint_t num) const {
  value_type data = m_data;
  data[deme] += num;
  return Allele(data);
}


bool Allele::singleton() const {
  return size() == 1;
}


vector<ExitAlleleData> Allele::reacheable() const {
  auto retvals = move_genes(*this);

  auto deme = this->deme();
  Init::value_type vals(deme);
  State::value_type s(deme * deme);
  ExitAlleleData archetype{value_type(deme), {vals, s}, 1.0};
  return combine_allele_parts(archetype, retvals.begin(), retvals.end());
}


Allele::const_iterator Allele::begin() const {
  return m_data.begin();
}


Allele::const_iterator Allele::end() const {
  return m_data.end();
}


bool Allele::operator<(Allele const& allele) const {
  return m_data < allele.m_data;
}


bool Allele::operator<=(Allele const& allele) const {
  return *this < allele or *this == allele;
}


bool Allele::operator>(Allele const& allele) const {
  return not (*this <= allele);
}


bool Allele::operator>=(Allele const& allele) const {
  return not (*this < allele);
}


esf_uint_t Allele::operator[](esf_uint_t deme) const {
  return m_data[deme];
}


bool operator==(Allele const& a, Allele const& b) {
  return a.m_data == b.m_data;
}


bool operator==(ExitAlleleData const& a, ExitAlleleData const& b) {
  return a.allele == b.allele && a.state == b.state;
}


ostream& operator<<(ostream& str, Allele const& allele) {

  str << "Allele(";
  copy(allele.begin(), allele.end(), ostream_iterator<esf_uint_t>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;

}


Allele operator+(Allele const& a, Allele const&b) {
  Allele::value_type data{a.begin(), a.end()};
  transform(data.begin(), data.end(), b.begin(), data.begin(), plus<esf_uint_t>());
  return Allele{data};
}


namespace {

// generate all assignments of genes within an allelic class.
vector<vector<ExitAlleleData>> move_genes(Allele const& orig) {

  vector<vector<ExitAlleleData>> parts;
  auto deme = orig.deme();
  // for (auto i = 0; i < 1; ++i) {
  for (esf_uint_t d = 0; d < deme; ++d) {
    Allele::value_type zeros(deme);
    State::value_type states(deme * deme);
    Init::value_type init(deme);
    init[d] = orig[d];
    parts.push_back(assign_genes(d, orig[d], 0U, init, zeros, states));
  }
  return parts;
}


// generate all assignments of genes from a single deme within an
// alleleic class.
vector<ExitAlleleData> assign_genes(esf_uint_t origin,
                                    esf_uint_t count,
                                    esf_uint_t deme,
                                    Init const& init,
                                    Allele const allele,
                                    State::value_type& svec) {
  auto dim = allele.deme();

  // base case: no degree of freedom left
  if (deme == dim - 1) {
    State::value_type s{svec};
    s[origin * dim + deme] = count;
    return {{allele.add(deme, count), {init, s}, 1.0}};
  }

  // regular cases: factor increases by the number of available
  // aassignments.

  vector<ExitAlleleData> holder;
  for (esf_uint_t i = 0; i <= count; ++i) {
    Allele base = allele.add(deme, i);
    State::value_type s{svec};
    s[origin * dim + deme] += i;
    auto retval = assign_genes(origin, count - i, deme + 1, init, base, s);
    holder.insert(holder.end(), retval.begin(), retval.end());
  }

  return holder;
}


vector<ExitAlleleData> combine_allele_parts(ExitAlleleData const& archetype,
                                            vector<vector<ExitAlleleData>>::const_iterator begin,
                                            vector<vector<ExitAlleleData>>::const_iterator end) {
  if (begin == end) {
    ExitAlleleData a(archetype);
    auto deme = a.allele.deme();
    for (esf_uint_t i = 0; i < deme; ++i) {
      auto num = a.allele[i];
      for (esf_uint_t j = 0; j < deme; ++j) {
        auto k = a.state[deme * j + i];
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


ExitAlleleData synthesize_genes(ExitAlleleData const& a,
                                ExitAlleleData const& b) {
  return {a.allele + b.allele, a.state + b.state, a.factor * b.factor};
}


}


}
