// -*- mode: c++; coding: utf-8; -*-

// esf_prob.cc - Compute probability of structured ESF.

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

#include "afs.hh"
#include "allele.hh"
#include "enum.hh"
#include "esf_prob.hh"

namespace esf {


ESFProb::ESFProb(AFS a, Params p)
    : afs(a), init(a.to_init()), ndeme(a.demes()), params(p) {

  hitprob = HitProb(init, params);

}


Value ESFProb::compute() {

  if (afs.singleton()) {

    return compute_with_singleton();

  } else {

    return compute_without_singleton();

  }

}


Value ESFProb::compute_without_singleton() {

  if (afs.size() == 1) {

    return 1.0;

  }

  Value val = 0.0;

  HitProb hp = HitProb(init, params);

  AFS other;
  Allele na;

  for (auto spec: afs.reacheable()) {

    val += compute_coal_probs(spec);

  }

}


Value ESFProb::compute_with_singleton() {

  auto singleton = ::std::find_if(afs.begin(), afs.end(),
                                  [](AFS::value_type p)
                                  {
                                    return p.first.singleton();
                                  });

  Allele na, allele = *singleton;

  Index deme = 0;
  while (allele[deme] == 0) {

    ++deme;

  }

  Index dsize = afs.size(deme);

  AFS base = afs.remove(allele).add(allele.remove(deme)), other;

  Value val = ESFProb(base, params).compute();

  for (auto a: base) {

    na = a.add(deme);

    other = base.remove(a).add(na);

    val -= ESFProb(other, params).compute() * other[na] * na.size(deme) / dsize;

  }

  val *= allele.size(deme) / dsize;

  return val;

}


Value compute_coal_probs(AFS afs) {

  Value val = 0.0;
  State state = afs.to_state();
  AFS other0, other1;
  Allele na;

  for (auto allele: afs) {

    other0 = afs.remove(allele);

    for (Index i = 0; i < ndeme; ++i) {

      if (allele.size(i) > 1) {

        na = allele.remove(i);

        other1 = other0.add(na);

        val += other1.compute() * other1.count(na) * na.size(i) / other1.size(i);

      }

    }

  }

  return val;

}

};
