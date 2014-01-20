// -*- mode: c++; coding: utf-8; -*-

// afs.cc - Implementation of allele frequency spectrum

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
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "afs.hh"
#include "allele.hh"
#include "init.hh"
#include "util.hh"

namespace esf {


AFS::AFS(::std::vector<Allele> const& avec) {

  for (auto allele: avec) {

    m_data[allele] += 1;

  }

}


AFS::AFS(AFS::data_type const& data)
    : m_data(data) {}


AFS AFS::add(Allele const& allele) const {

  auto data = m_data;

  if (allele.size() > 0) {

    ++data[allele];

  }

  return AFS(data);

}


AFS AFS::remove(Allele const& allele) const {

  auto data = m_data;

  if (allele.size() > 0) {

    if (data[allele] > 1) {

      --data[allele];

    } else {

      data.erase(allele);

    }

  }

  return AFS(data);

}


bool AFS::singleton() const {

  return ::std::any_of(m_data.begin(), m_data.end(),
                       [](value_type p)
                       {
                         return p.first.singleton();
                       }
                       );

}


Index AFS::operator[](const Allele& allele) const {

  try {

    return m_data.at(allele);

  } catch (::std::out_of_range&) {

    return 0;

  }

}


Index AFS::size(Index deme) const {

  Index start = 0;
  return ::std::accumulate(m_data.begin(), m_data.end(), start,
                           [deme](Index a, value_type p)
                           {
                             return a + (p.first)[deme] * p.second;
                           }
                           );

}


Index AFS::size() const {

  Index start = 0;
  return ::std::accumulate(this->begin(), this->end(), start,
                           [](Index a, value_type p)
                           {
                             return a + (p.first).size() * p.second;
                           }
                           );
}


Index AFS::deme() const {

  return m_data.begin()->first.deme();

}


::std::vector<ExitAFSData> AFS::reacheable() const {

  using ::std::vector;

  vector<Allele> alleles;
  vector<Index> state_vec(unsign(deme() * deme()));

  return build(alleles, state_vec, m_data.begin(), m_data.end());

}


::std::vector<ExitAFSData> AFS::build(::std::vector<Allele> alleles,
                                      ::std::vector<Index> states,
                                      data_type::const_iterator begin,
                                      data_type::const_iterator end) const {

  if (begin == end) {

    return {ExitAFSData({AFS(alleles), State(Init(*this), states)})};

  }

  auto reacheables = begin->first.reacheable();

  auto itr_b = reacheables.begin();
  auto itr_e = reacheables.end();

  return sub_build(alleles, states, begin, end, begin->second, itr_b, itr_e);

}


::std::vector<ExitAFSData>
AFS::sub_build(::std::vector<Allele> alleles,
               ::std::vector<Index> states,
               data_type::const_iterator begin,
               data_type::const_iterator end,
               Index count,
               ::std::vector<ExitAlleleData>::const_iterator allele_begin,
               ::std::vector<ExitAlleleData>::const_iterator allele_end) const {

  if (count == 0 || allele_begin == allele_end) {

    auto begin_copy = begin;

    ++begin_copy;

    return build(alleles, states, begin_copy, end);

  }

  using ::std::vector;

  vector<ExitAFSData> retval;

  for (auto a_itr = allele_begin; a_itr != allele_end; ++a_itr) {

    auto a_copy = alleles;
    auto s_copy = states;

    a_copy.push_back(a_itr->allele);

    auto s = a_itr->state;
    for (decltype(s.size()) i = 0; i != s.size(); ++i) {

      s_copy[i] += s[i];

    }

    auto p = sub_build(a_copy, s_copy, begin, end, count - 1, a_itr, allele_end);

    retval.insert(retval.end(), p.begin(), p.end());

  }

  return retval;

}


AFS::iterator AFS::begin() {

  return m_data.begin();

}


AFS::const_iterator AFS::begin() const {

  return m_data.begin();

}


AFS::iterator AFS::end() {

  return m_data.end();

}


AFS::const_iterator AFS::end() const {

  return m_data.end();

}


bool operator==(AFS const& a, AFS const& b) {

  return a.m_data == b.m_data;

}


bool operator<(AFS const& a, AFS const& b) {

  return a.m_data < b.m_data;

}


bool operator==(ExitAFSData const& a, ExitAFSData const& b) {

  return a.afs == b.afs && a.state == b.state;

}


bool operator<(ExitAFSData const& a, ExitAFSData const& b) {

  if (a.afs == b.afs) {

    return a.state < b.state;

  }

  return a.afs < b.afs;

}


::std::ostream& operator<<(::std::ostream& str, ::std::pair<Allele, Index> const& p) {
  using ::std::copy;
  using ::std::ostream_iterator;

  str << p.first << ": " << p.second;

  return str;

}


::std::ostream& operator<<(::std::ostream& str, AFS const& afs) {

  using ::std::copy;
  using ::std::ostream_iterator;
  using ::std::pair;

  str << "AFS(";
  copy(afs.begin(), afs.end(), ostream_iterator<pair<Allele, Index>>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;

}

}
