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
#include <numeric>
#include <stdexcept>
#include <vector>

#include "afs.hh"
#include "allele.hh"
#include "init.hh"

namespace esf {


AFS::AFS(AFS const& afs, Allele const& allele, Mode mode)
    : data(afs.data), m_deme(afs.m_deme) {

  if (allele.size() > 0) {

    if (mode == Mode::ADD) {

      data[allele] += 1;

    } else if (data[allele] > 1) {

      data[allele] -= 1;

    } else {

      data.erase(allele);

    }

  }

}


AFS::AFS(::std::vector<Allele> const& avec) {

  if (avec.size() > 0) {

    m_deme = avec[0].deme();

  }

  for (auto allele: avec) {

    data[allele] += 1;

  }

}


AFS AFS::add(Allele allele) {

  return AFS(*this, allele, Mode::ADD);

}


AFS AFS::remove(Allele allele) {

  return AFS(*this, allele, Mode::REMOVE);

}


bool AFS::singleton() const {

  return ::std::any_of(data.begin(), data.end(),
                       [](value_type p)
                       {
                         return p.first.singleton();
                       }
                       );

}


Index AFS::operator[](const Allele& allele) const {

  try {

    return data.at(allele);

  } catch (::std::out_of_range& e) {

    return 0;

  }

}


Index AFS::size(Index deme) const {

  return ::std::accumulate(data.begin(), data.end(), 0,
                           [deme](Index a, value_type p)
                           {
                             return a + (p.first)[deme] * p.second;
                           }
                           );

}


Index AFS::size() const {

  return ::std::accumulate(this->begin(), this->end(), 0,
                           [](Index a, value_type p)
                           {
                             return a + (p.first).size() * p.second;
                           }
                           );
}


Index AFS::deme() const {

  return m_deme;

}


::std::vector<ExitAFSPair> AFS::reacheable() const {

  ::std::vector<Allele> alleles;
  ::std::vector<Index> state_vec(m_deme * m_deme);

  return build(alleles, state_vec, data.begin(), data.end());

}


::std::vector<ExitAFSPair> AFS::build(::std::vector<Allele> alleles,
                                      ::std::vector<Index> states,
                                      data_type::const_iterator begin,
                                      data_type::const_iterator end) const {

  if (begin == end) {

    return {ExitAFSPair({AFS(alleles), State(Init(*this), states)})};

  }

  auto reacheables = begin->first.reacheable();

  auto itr_b = reacheables.begin();
  auto itr_e = reacheables.end();

  return sub_build(alleles, states, begin, end, begin->second, itr_b, itr_e);

}


::std::vector<ExitAFSPair>
AFS::sub_build(::std::vector<Allele> alleles,
               ::std::vector<Index> states,
               data_type::const_iterator begin,
               data_type::const_iterator end,
               Index count,
               ::std::vector<ExitAllelePair>::const_iterator allele_begin,
               ::std::vector<ExitAllelePair>::const_iterator allele_end) const {

  if (count == 0 || allele_begin == allele_end) {

    auto begin_copy = begin;

    ++begin_copy;

    return build(alleles, states, begin_copy, end);

  }

  ::std::vector<ExitAFSPair> retval;

  for (auto a_itr = allele_begin; a_itr != allele_end; ++a_itr) {

    auto a_copy = alleles;
    auto s_copy = states;

    a_copy.push_back(a_itr->allele);

    auto s = a_itr->state;
    for (auto i = 0; i != s.size(); ++i) {

      s_copy[i] += s[i];

    }

    auto p = sub_build(a_copy, s_copy, begin, end, count - 1, a_itr, allele_end);

    retval.insert(retval.end(), p.begin(), p.end());

  }

  return retval;

}


AFS::iterator AFS::begin() {

  return data.begin();

}


AFS::const_iterator AFS::begin() const {

  return data.begin();

}


AFS::iterator AFS::end() {

  return data.end();

}


AFS::const_iterator AFS::end() const {

  return data.end();

}


bool operator==(AFS const& a, AFS const& b) {

  return a.data == b.data;

}


bool operator<(AFS const& a, AFS const& b) {

  return a.data < b.data;

}


bool operator==(ExitAFSPair const& a, ExitAFSPair const& b) {

  return a.afs == b.afs && a.state == b.state;

}


bool operator<(ExitAFSPair const& a, ExitAFSPair const& b) {

  if (a.afs == b.afs) {

    return a.state < b.state;

  }

  return a.afs < b.afs;

}


};
