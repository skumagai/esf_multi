// -*- mode: c++; coding: utf-8; -*-

// param.hh - demographic parameter object

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

#ifndef ESF_MULTI_PARAM_HH
#define ESF_MULTI_PARAM_HH


#include <vector>


#include "typedef.hh"


namespace esf {


// This class stores demographic parameters.  Those parameters at
// present are migration rates, population size relative to
// (idealized) ancestral population, and mutation rates.  All
// parameters including mutation rates are subpopulation-specific.
class Param {

 public:

  typedef typename ::std::vector<double> value_type;

 private:

  value_type m_mig;

  value_type m_pop;

  value_type m_mut;

 public:

  Param() = default;

  Param(Param const&) = default;
  Param(Param&&) = default;

  Param& operator=(Param const&) = default;
  Param& operator=(Param&&) = default;

  ~Param() = default;

  // This constructor takes migration rate matrix, relative population
  // size vector, and mutation rate vector.  If there are n
  // subpopulations, respective size is n^2, n, and n.  Note that the
  // migration rate matrix needs to be flattened and that the matrix
  // has to be specified in column-major order (FORTRAN's storate mode
  // rather than C's).
  Param(value_type mi, value_type p, value_type mu)
      : m_mig(mi), m_pop(p), m_mut(mu) {}

  // This function takes two integers corresponding to source and
  // target demes, and it returns the corresponding migration rate.
  // Note that direction of migration is defined in backward in
  // time.  Deme index starts from zero.
  double mig_rate(esf_uint_t, esf_uint_t) const;

  double& mig_rate(esf_uint_t, esf_uint_t);

  // This function returns a mutation rate in a deme.
  double mut_rate(esf_uint_t) const;

  double& mut_rate(esf_uint_t);

  // This function returns a relative population size of a deme.
  double pop_size(esf_uint_t) const;

  double& pop_size(esf_uint_t);

};


}


#endif // ESF_MULTI_PARAM_HH
