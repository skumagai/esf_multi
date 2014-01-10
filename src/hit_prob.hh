// -*- mode: c++; coding: utf-8; -*-

// hit_prob.hh - Computing hitting probabilities

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

#ifndef ESF_MULTI_HIT_PROB_HH
#define ESF_MULTI_HIT_PROB_HH

#include <vector>

#include "init.hh"
#include "param.hh"
#include "typedef.hh"

namespace esf {

using ::std::vector;


// Forward declaration
class Init;
class State;


// This class computes the probabilities of coalescence given initial
// locations of genes and demographic parameters. Genotype of genes is
// not considered. The probabilities are computed for location of
// genes immediately prior to coalescence and the location of
// coalescent event.
class HitProb {

 public:

  typedef typename ::std::vector<double> value_type;

 private:

  Init m_init;

  Param m_param;

  vector<double> m_prob;

  void compute();

  double compute_u(State const&, State const&) const;

 public:

  // The hittng probabilities are computated upon construction of
  // HitProb object.
  HitProb(Init const&, Param const&);

  // Returns the hitting probability of i-th state and coalescence in
  // j-th deme. States contain information of initial and current
  // placement of genes, but they do not contain information on
  // coalescence.
  double get(Index, Index) const;

  // Returns the hitting probability of specified state and deme.
  double get(State const&, Index) const;

  HitProb update(Param const&) const;

};


}


#endif // ESF_MULTI_HIT_PROB_HH
