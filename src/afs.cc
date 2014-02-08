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
#include <map>
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
using ::std::make_pair;
using ::std::map;
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


esf_uint_t AFS::operator[](const Allele& allele) const {
  try {

    return m_data.at(allele);

  } catch (out_of_range&) {

    return 0U;

  }
}


esf_uint_t AFS::size(esf_uint_t deme) const {
  return accumulate(m_data.begin(), m_data.end(), 0U,
                    [deme](esf_uint_t a, value_type p) {
                      return a + (p.first)[deme] * p.second;
                    });

}


esf_uint_t AFS::size() const {
  return accumulate(this->begin(), this->end(), 0U,
                    [](esf_uint_t a, value_type p) {
                      return a + (p.first).size() * p.second;
                    });
}


esf_uint_t AFS::deme() const {
  return m_data.begin()->first.deme();
}


vector<ExitAFSData> AFS::reacheable() const {
  auto deme = this->deme();
  AFS afs{};
  Init::value_type i(deme);
  State::value_type s(deme * deme);
  State state{i, s};
  vector<pair<Allele,Allele>> match;
  return build(afs, state, match, 1.0, m_data.begin(), m_data.end());
}


vector<ExitAFSData> AFS::build(AFS const& afs,
                               State const& state,
                               vector<pair<Allele, Allele>> const& match,
                               double factor,
                               data_type::const_iterator begin,
                               data_type::const_iterator end) const {
  if (begin == end) {
    auto denom = get_denominator(Init{afs}, state);
    auto duplicity = get_duplicity(afs, match);
    return {ExitAFSData({afs, state, duplicity * factor / denom})};
  } else {
    auto reacheables = begin->first.reacheable();
    auto itr_b = reacheables.begin();
    auto itr_e = reacheables.end();
    return sub_build(afs, state, match, factor, begin, end, begin->second, itr_b, itr_e);
  }
}


vector<ExitAFSData>
AFS::sub_build(AFS const& afs,
               State const& state,
               vector<pair<Allele, Allele>> const& match,
               double factor,
               data_type::const_iterator begin,
               data_type::const_iterator end,
               esf_uint_t count,
               vector<ExitAlleleData>::const_iterator allele_begin,
               vector<ExitAlleleData>::const_iterator allele_end) const {
  if (count == 0) {
    ++begin;
    return build(afs, state, match, factor, begin, end);
  } else {
    vector<ExitAFSData> retval;
    for (auto a_itr = allele_begin; a_itr != allele_end; ++a_itr) {
      auto m(match);
      m.push_back(make_pair(begin->first, a_itr->allele));
      auto p = sub_build(afs.add(a_itr->allele),
                         a_itr->state + state,
                         m,
                         factor * a_itr->factor,
                         begin, end, count - 1, a_itr, allele_end);
      retval.insert(retval.end(), p.begin(), p.end());
    }
    return retval;
  }
}


double AFS::get_denominator(Init const& init, State const& state) const {

  auto deme = this->deme();
  double denom = 1.0;
  for (esf_uint_t i = 0; i < deme; ++i) {
    esf_uint_t genes = init[i];
    for (esf_uint_t j = 0; j < deme; ++j) {
      esf_uint_t in_the_deme = state[j * deme + i];
      denom *= binomial(genes, in_the_deme);
      genes -= in_the_deme;
    }
  }
  return denom;
}


double AFS::get_duplicity(AFS const& past, vector<pair<Allele, Allele>> const& a) const {

  double val = 1.0;

  map<Allele, esf_uint_t> count_now;
  for (auto& allele: m_data) {
    count_now[allele.first] = allele.second;
  }

  map<Allele, esf_uint_t> count_prev;
  for (auto& allele: past) {
    count_prev[allele.first] = allele.second;
  }

  AFS tmp{{Allele{{0, 2}}, Allele{{1, 1}}, Allele{{1,1}}, Allele{{2, 0}}}};

  if (past == tmp) {
    std::cout << "ENTER\n";
  }

  for (auto& p: a) {
    // The first element is an allele in the current state, and the
    // second element is an allele immediately before an event.
    auto& now = p.first;
    auto& prev = p.second;

    esf_uint_t choice = count_now[now];
    if (choice > count_prev[prev]) {
      choice = count_prev[prev];
    }

    auto v = static_cast<double>(binomial(count_prev[prev], choice));
    if (v == 0.0) {
      std::cout << count_prev[prev] << " " << choice << "\n";
    }

    val *= static_cast<double>(binomial(count_prev[prev], choice));

    if (past == tmp) {
      std::cout << now << " " << prev << " " << count_prev[prev] << " " << choice << " " << val << "\n";
    }
    count_prev[prev] -= choice;
    count_prev[now] -= 1;
  }

  if (past == tmp) {
    std::cout << "EXIT\n\n";
  }

  return val;
}


AFS::const_iterator AFS::begin() const {
  return m_data.begin();
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


ostream& operator<<(ostream& str, pair<Allele, esf_uint_t> const& p) {
  str << p.first << ": " << p.second;
  return str;
}


ostream& operator<<(ostream& str, AFS const& afs) {
  auto end = afs.end();
  --end;
  str << "AFS(";
  copy(afs.begin(), end, ostream_iterator<pair<Allele, esf_uint_t>>(str, ", "));
  str << *end;
  str << ")";

  return str;
}

}
