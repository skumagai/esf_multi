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

using ::std::accumulate;
using ::std::any_of;
using ::std::copy;
using ::std::ostream;
using ::std::ostream_iterator;
using ::std::out_of_range;
using ::std::pair;
using ::std::plus;
using ::std::transform;
using ::std::vector;

AFS::AFS(vector<Allele> const& avec) {

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

  return any_of(m_data.begin(), m_data.end(),
                [](value_type p) {
                  return p.first.singleton();
                });

}


Index AFS::operator[](const Allele& allele) const {

  try {

    return m_data.at(allele);

  } catch (out_of_range&) {

    return 0;

  }

}


Index AFS::size(Index deme) const {

  Index start = 0;
  return accumulate(m_data.begin(), m_data.end(), start,
                    [deme](Index a, value_type p) {
                      return a + (p.first)[deme] * p.second;
                    });

}


Index AFS::size() const {

  Index start = 0;
  return accumulate(this->begin(), this->end(), start,
                    [](Index a, value_type p) {
                      return a + (p.first).size() * p.second;
                    });
}


Index AFS::deme() const {

  return m_data.begin()->first.deme();

}


vector<ExitAFSData> AFS::reacheable() const {

  vector<Allele> alleles;
  vector<Index> state_vec(unsign(deme() * deme()));

  return build(alleles, state_vec, 1.0, m_data.begin(), m_data.end());
}


vector<ExitAFSData> AFS::build(vector<Allele> const& alleles,
                               vector<Index> const& states,
                               double factor,
                               data_type::const_iterator begin,
                               data_type::const_iterator end) const {

  if (begin == end) {
    Init init{*this};
    State state{init, states};
    auto denom = get_denominator(init, state);
    return {ExitAFSData({AFS(alleles), state, factor / denom})};
  }

  // auto reacheables = begin->first.reacheable();

  // auto itr_b = reacheables.begin();
  // auto itr_e = reacheables.end();

  // return sub_build(alleles, states, factor, begin, end,
  //                  begin->second, itr_b, itr_e);
  return sub_build(alleles, states, factor, begin, end, begin->second);
}


vector<ExitAFSData>
AFS::sub_build(vector<Allele> const& alleles,
               vector<Index> const& states,
               double factor,
               data_type::const_iterator begin,
               data_type::const_iterator end,
               Index count) const {

  // if (count == 0 || allele_begin == allele_end) {
  if (count == 0) {
    auto begin_copy = begin;
    ++begin_copy;
    return build(alleles, states, factor, begin_copy, end);
  }


  vector<ExitAFSData> retval;

  for (auto a: (begin->first).reacheable()) {

  // for (auto a_itr = allele_begin; a_itr != allele_end; ++a_itr) {

    auto a_copy (alleles);
    auto s_copy(states);

    // a_copy.push_back(a_itr->allele);
    a_copy.push_back(a.allele);

    // auto s = a_itr->state;
    auto s = a.state;
    for (decltype(s.size()) i = 0; i != s.size(); ++i) {

      s_copy[i] += s[i];

    }

    std::cout << "[DEBUG]: " << a.allele << "\t" << s << "\n";

    // auto p = sub_build(a_copy, s_copy, factor * a_itr->factor,
    //                    begin, end, count - 1, a_itr, allele_end);
    auto p = sub_build(a_copy, s_copy, factor * a.factor, begin, end, count - 1);

    retval.insert(retval.end(), p.begin(), p.end());

  }

  return retval;

}


double AFS::get_denominator(Init const& init, State const& state) const {

  auto deme = this->deme();
  double denom = 1.0;
  for (auto i = 0; i < deme; ++i) {
    Index genes = init[i];
    for (auto j = 0; j < deme; ++j) {
      auto in_the_deme = state[i * deme + j];
      denom *= binomial(genes, in_the_deme);
      genes -= in_the_deme;
    }
  }
  return denom;
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
  } else {
    return a.afs < b.afs;
  }

}


ostream& operator<<(ostream& str, pair<Allele, Index> const& p) {

  str << p.first << ": " << p.second;
  return str;

}


ostream& operator<<(ostream& str, AFS const& afs) {

  str << "AFS(";
  copy(afs.begin(), afs.end(), ostream_iterator<pair<Allele, Index>>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;

}

}
