// -*- mode: c++; coding: utf-8; -*-

// esf_prob.hh - Compute probability of structured ESF

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

#ifndef ESF_MULTI_ESF_PROB_HH
#define ESF_MULTI_ESF_PROB_HH

#include "typedef.hh"
#include "afs.hh"
#include "hit_prob.hh"
#include "param.hh"
#include "typedef.hh"

namespace esf {


template <typename KEY, typename VALUE> class Cache;


// An instance of this class computes probabilities of an allele
// frequency spectrum according to population-structured extension to
// Ewen's sampling formula.
class ESFProb {

 private:

  AFS const m_afs;

  Init const m_init;

  Param const m_param;

  Cache<AFS, double>*  m_esf_prob_cache;

  Cache<Init, HitProb>* m_hit_prob_cache;

  double compute_with_singleton();

  double compute_without_singleton();

  double compute_coal_probs(ExitAFSPair const&, HitProb const&);

 public:

  // ESFProb() = delete;

  // ESFProb(ESFProb const&);

  ESFProb& operator=(ESFProb const&);

  ~ESFProb();

  // This constructor takes an allele frequency spectrum (of AFS
  // class) and demographic paramegers of (Param class).  The acutual
  // computation is deferred until compute method is explicitly invoked.
  ESFProb(AFS const&, Param const&);

  ESFProb(AFS const&, Param const&, Cache<AFS, double>*, Cache<Init, HitProb>*);

  // This function implements actual computation of
  // population-structured ESP, and it's return value is a probability
  // according to the sample state (allele frequency spectrum) and
  // demographic parameters.  The computation is performed recursively.
  double compute();

  friend void swap(ESFProb&, ESFProb&);

};


}


#endif // ESF_MULTI_ESF_PROB_HH
