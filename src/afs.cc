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
using ::std::make_pair;
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
                [](value_type p) { return p.first.singleton(); });
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
                    [](esf_uint_t a, value_type p) { return a + (p.first).size() * p.second; });
}


esf_uint_t AFS::deme() const {
  return m_data.begin()->first.deme();
}


vector<ExitAFSData> AFS::reacheable() const {
  esf_uint_t deme = this->deme();
  esf_uint_t size = m_data.size();

  vector<esf_uint_t> count;
  count.reserve(size);

  vector<vector<ExitAlleleData>> a_store;
  a_store.reserve(size);

  vector<esf_uint_t> dim;
  dim.reserve(size);

  vector<Allele> orig;
  orig.reserve(size);

  for (auto& p: m_data) {
    count.push_back(p.second);
    orig.push_back(p.first);

    a_store.push_back(p.first.reacheable());

    esf_uint_t k = a_store.back().size();
    dim.push_back(binomial(count.back() + k - 1, count.back()));
  }

  vector<vector<vector<esf_uint_t>>> indecies(size);
  for (esf_uint_t i = 0; i < size; ++i) {
    indecies[i].reserve(dim[i]);
    vector<esf_uint_t> tmp(count[i]);
    vector<esf_uint_t> r(a_store[i].size());
    indecies[i].push_back(tmp);
    generate_n(r.begin(), r.size(), range<esf_uint_t>());
    while (next_k_selection(r.begin(), r.end(), tmp.begin(), tmp.end())) {
      indecies[i].push_back(tmp);
    }
  }

  esf_uint_t total_dim = accumulate(dim.begin(), dim.end(), 1UL, multiplies<esf_uint_t>());

  map<pair<AFS, State>, double> retmap;

  vector<esf_uint_t> jndex(size);
  for (esf_uint_t i = 0; i < total_dim; ++i) {
    index_1_to_n(dim.begin(), dim.end(), jndex.begin(), i);

    AFS afs;
    State state{Init{vector<esf_uint_t>(deme)}, vector<esf_uint_t>(deme * deme)};

    double factor = 1.0;
    map<Allele, map<pair<Allele, State>, esf_uint_t>> repmap;
    for (esf_uint_t j = 0; j < size; ++j) {
      auto& index = indecies[j][jndex[j]];
      for (auto k: index) {
        auto& a = a_store[j][k];
        afs = afs.add(a.allele);
        state = state + a.state;
        factor *= a.factor;
        ++repmap[a.allele][make_pair(orig[j], a.state)];
      }
    }
    for (auto& p: repmap) {
      auto m = p.second;
      vector<esf_uint_t> repvec(m.size());
      transform(m.begin(), m.end(), repvec.begin(),
                [](decltype(m)::value_type& q) { return q.second; });
      factor *= multinomial(repvec.begin(), repvec.end());
    }
    factor /= get_denominator(Init{afs}, state);
    retmap[make_pair(afs, state)] += factor;
  }

  vector<ExitAFSData> retvec(retmap.size());
  transform(retmap.begin(), retmap.end(), retvec.begin(),
            [](decltype(retmap)::value_type& p) {
              return ExitAFSData{p.first.first, p.first.second, p.second};
            });

  return retvec;
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
  str << "AFS(";
  copy(afs.begin(), afs.end(), ostream_iterator<pair<Allele, esf_uint_t>>(str, ", "));
  str << "\b\b";
  str << ")";

  return str;
}

}
