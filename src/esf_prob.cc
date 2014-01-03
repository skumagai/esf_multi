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
#include "hit_prob.hh"
#include "esf_prob.hh"

namespace esf {


ESFProb::ESFProb(AFS a, Param p)
    : afs(a), init(a), ndeme(a.deme()), param(p) {}


double ESFProb::compute() {

  if (afs.singleton()) {

    return compute_with_singleton();

  } else {

    return compute_without_singleton();

  }

}


double ESFProb::compute_without_singleton() {

  double val = 0.0;

  HitProb hp(init, param);

  for (auto spec: afs.reacheable()) {

    val += compute_coal_probs(spec, hp);

  }

  return val;

}


Value ESFProb::compute_with_singleton() {

  if (afs.size() == 1) {

    return 1.0;

  }

  using ::std::find_if;

  auto singleton = find_if(afs.begin(), afs.end(),
                           [](AFS::value_type p)
                           {
                             return p.first.singleton();
                           });

  Allele allele = singleton->first;

  Index deme = 0;
  while (allele[deme] == 0) {

    ++deme;

  }

  double dsize = static_cast<double>(afs.size(deme));

  AFS base = afs.remove(allele).add(allele.remove(deme));

  // probability of a sample excluding one of singleton alleles.
  double val = ESFProb(base, param).compute();

  for (auto a: base) {

    // add a gene to previously an observed allele.
    Allele na = a.first.add(deme);

    AFS other = base.remove(a.first).add(na);

    val -= ESFProb(other, param).compute() * other[na] * na[deme] / dsize;


  }

  val *= dsize / afs[allele];

  return val;

}


double ESFProb::compute_coal_probs(ExitAFSPair const& pair, HitProb const& hp) {

  double val = 0.0;
  AFS afs = pair.afs;
  State state = pair.state;

  for (auto a: afs) {

    auto allele = a.first;

    AFS other0 = afs.remove(allele);

    for (Index i = 0; i < ndeme; ++i) {

      if (allele[i] > 1) {

        Allele na = allele.remove(i);

        AFS other1 = other0.add(na);

        double tmp = hp.get(state, i) * ESFProb(other1, param).compute();

        tmp *= na.size() * other1[na] / other1.size();

        val += tmp;

      }

    }

  }

  return val;

}


};
