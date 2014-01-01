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


class Param {

 public:

  typedef typename ::std::vector<double> value_type;

 private:

  value_type m_pop;

  value_type m_mut;

  value_type m_mig;

 public:

  Param(value_type mi, value_type p, value_type mu)
      : m_mig(mi), m_pop(p), m_mut(mu) {}

  // migration rate is stored in column-major order
  double& mig_rate(Index, Index);

  double& mut_rate(Index);

  double& pop_size(Index);

};


}


#endif // ESF_MULTI_PARAM_HH
