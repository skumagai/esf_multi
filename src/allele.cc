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

using ::std::make_pair;
using ::std::pair;
using ::std::plus;
using ::std::size_t;
using ::std::transform;
using ::std::vector;

struct AlleleData_ {
  vector<Index> allele;
  double factor;
};


pair<AlleleData_, vector<Index>> synthesize_genes(AlleleData_ const&,
                                                  vector<Index> const&,
                                                  AlleleData_ const&);

vector<vector<AlleleData_>> move_genes(Allele const&);

vector<AlleleData_> assign_genes(Index, AlleleData_ const&, size_t);

vector<ExitAlleleData> combine_allele_parts(AlleleData_ const&,
                                            vector<Index> const&,
                                            vector<vector<AlleleData_>>::const_iterator,
                                            vector<vector<AlleleData_>>::const_iterator);

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


::std::vector<ExitAlleleData> Allele::reacheable() const {

  using ::std::vector;

  auto retvals = move_genes(*this);

  AlleleData_ archetype{vector<Index>(unsign(this->deme())), 1.0};
  vector<Index> state;
  return combine_allele_parts(archetype, state, retvals.begin(), retvals.end());
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


::std::ostream& operator<<(::std::ostream& str, Allele const& allele) {

  using ::std::copy;
  using ::std::ostream_iterator;

  str << "Allele(";
  copy(allele.begin(), allele.end(), ostream_iterator<Index>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;

}



namespace {

// generate all assignments of genes within an allelic class.
vector<vector<AlleleData_>> move_genes(Allele const& orig) {

  vector<vector<AlleleData_>> parts;
  for (auto i: orig) {
    vector<Index> zeros(unsign(orig.deme()));
    parts.push_back(assign_genes(i, {zeros, 1.0}, static_cast<size_t>(0)));
  }

  return parts;
}


pair<AlleleData_, vector<Index>> synthesize_genes(AlleleData_ const& archetype,
                                                  vector<Index> const& state,
                                                  AlleleData_ const& new_info) {
  auto begin = new_info.allele.begin();
  auto end = new_info.allele.end();

  AlleleData_ d(archetype);
  transform(d.allele.begin(), d.allele.end(), begin, d.allele.begin(), plus<Index>());

  vector<Index> s(state);
  s.insert(s.end(), begin, end);

  d.factor *= new_info.factor;

  return make_pair(d, s);
}


vector<ExitAlleleData> combine_allele_parts(AlleleData_ const& archetype,
                                            vector<Index> const& state,
                                            vector<vector<AlleleData_>>::const_iterator begin,
                                            vector<vector<AlleleData_>>::const_iterator end) {
  if (begin == end) {
    // for (auto data: *begin) {
    //   auto ret = synthesize_genes(archetype, state, data);
    //   retval.push_back({ret.first.allele, ret.second, ret.first.factor});
    // }
    // return retval;

    return {{archetype.allele, state, archetype.factor}};
  }

  vector<ExitAlleleData> retval;
  for (auto data: *begin) {
    auto ret = synthesize_genes(archetype, state, data);
    auto allele_list = combine_allele_parts(ret.first, ret.second, begin + 1, end);
    retval.insert(retval.end(), allele_list.begin(), allele_list.end());
  }
  return retval;
}


// generate all assignments of genes from a single deme within an
// alleleic class.
vector<AlleleData_> assign_genes(Index count, AlleleData_ const& archetype, size_t deme) {

  auto orig_data = archetype.allele;
  auto factor = archetype.factor;

  // base case: no degree of freedom left
  if (deme == orig_data.size() - 1) {
    decltype(orig_data) base(orig_data);
    base[deme] = count;
    return {AlleleData_{base, factor}};
  }

  // regular cases: factor increases by the number of available
  // aassignments.
  vector<AlleleData_> holder;
  for (Index i = 0; i <= count; ++i) {
    decltype(orig_data) base(orig_data);
    base[deme] = i;
    auto retval = assign_genes(count - i,
                               {base, factor * binomial(count, i)},
                               deme + 1);
    holder.insert(holder.end(), retval.begin(), retval.end());
  }

  return holder;
}


}


}
