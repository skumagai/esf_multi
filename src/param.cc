// -*- mode: c++; coding: utf-8; -*-

// param.cc - demographic parameter object

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

#include "param.hh"

namespace esf {


double Param::mig_rate(esf_uint_t i, esf_uint_t j) const {
  return m_mig[i + j * m_pop.size()];
}


double& Param::mig_rate(esf_uint_t i, esf_uint_t j) {
  return m_mig[i + j * m_pop.size()];
}


double Param::mut_rate(esf_uint_t i) const {
  return m_mut[i];
}


double& Param::mut_rate(esf_uint_t i) {
  return m_mut[i];
}


double Param::pop_size(esf_uint_t i) const {
  return m_pop[i];
}


double& Param::pop_size(esf_uint_t i) {
  return m_pop[i];
}


}
