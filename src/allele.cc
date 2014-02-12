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
  esf_uint_t deme = this->deme();

  vector<esf_uint_t> dim(deme);
  dimensions(m_data.begin(), m_data.end(), dim.begin());

  vector<esf_uint_t> offset(deme);
  offsets(dim.begin(), dim.end(), offset.begin());

  esf_uint_t nelem = accumulate(dim.begin(), dim.end(), 1UL, multiplies<esf_uint_t>());
  vector<vector<esf_uint_t>> alleles(nelem);
  for (auto& a: alleles) {
    a.resize(deme);
  }
  vector<vector<esf_uint_t>> dist(nelem);
  for (auto& d: dist) {
    d.resize(deme * deme);
  }

  vector<esf_uint_t> tmp(deme);
  for (esf_uint_t d = 0; d < deme; ++ d) {
    esf_uint_t k = m_data[d];
    for (esf_uint_t j = 0; j < dim[d]; ++j) {
      esf_uint_t i = j * offset[d];
      distribute_n(tmp.begin(), tmp.end(), j, k);
      while (i < nelem) {
        for (esf_uint_t dummy = 0; dummy < offset[d]; ++i, ++dummy) {
          transform(tmp.begin(), tmp.end(),
                    alleles[i].begin(), alleles[i].begin(), plus<esf_uint_t>());
          copy(tmp.begin(), tmp.end(), dist[i].begin() + sign(d * deme));
        }
        i += (dim[d] - 1) * offset[d];
      }
    }
  }

  vector<ExitAlleleData> retval;
  retval.reserve(nelem);
  for (esf_uint_t i = 0; i < nelem; ++i) {
    double factor = 1.0;
    for (esf_uint_t j = 0; j < deme; ++j) {
      esf_uint_t c = alleles[i][j];
      for (esf_uint_t k = 0; k < deme; ++k) {
        esf_uint_t choice = dist[i][j + k * deme];
        factor *= binomial(c, choice);
        c -= choice;
      }
    }
    retval.push_back({Allele{alleles[i]}, State{m_data, dist[i]}, factor});
  }

  return retval;
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


}
