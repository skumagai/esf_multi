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
#include "hitprob.hh"

namespace esef {


class ESFProb {

 private:

  AFS afs;

  Init init;

  Index ndeme;

  Params params;

  HitProb hitprob;

  Value compute_with_singleton();

  Value compute_without_singleton();

 public:

  ESFProb(AFS, Params);

  Value compute();

};


};


#endif // ESF_MULTI_ESF_PROB_HH
