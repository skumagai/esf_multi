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

#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "afs.hh"
#include "allele.hh"
#include "cache.hh"
#include "hit_prob.hh"
#include "esf_prob.hh"
#include "util.hh"

namespace esf {


ESFProb::~ESFProb() {

  if (&(m_esf_prob_cache->root()) == &m_afs) {

    delete m_esf_prob_cache;
    delete m_hit_prob_cache;

  }

}


ESFProb::ESFProb(AFS const& a, Param const& p)
    : m_afs(a), m_init(a), m_param(p),
      m_esf_prob_cache(new Cache<AFS, double>(m_afs)),
      m_hit_prob_cache(new Cache<Init, HitProb>(m_init)) {}


ESFProb::ESFProb(AFS const& a, Param const& p,
                 Cache<AFS, double>* esf_prob_cache,
                 Cache<Init, HitProb>* hit_prob_cache)
    : m_afs(a), m_init(a), m_param(p),
      m_esf_prob_cache(esf_prob_cache),
      m_hit_prob_cache(hit_prob_cache) {}


double ESFProb::compute() {

  try {

    return m_esf_prob_cache->at(m_afs);

  } catch (::std::out_of_range&) {

    double val;

    if (m_afs.singleton()) {

      val = compute_with_singleton();


    } else {

      val = compute_without_singleton();

    }

    (*m_esf_prob_cache)[m_afs] = val;

    return val;

  }

}


double ESFProb::compute_without_singleton() {

  double val = 0.0;

  HitProb hp;

  try {

    hp = m_hit_prob_cache->at(m_init);

  } catch (::std::out_of_range&) {

    hp = HitProb(m_init, m_param);

    (*m_hit_prob_cache)[m_init] = hp;

  }

  for (auto spec: m_afs.reacheable()) {

    val += compute_coal_probs(spec, hp);

  }

  return val;

}


Value ESFProb::compute_with_singleton() {

  if (m_afs.size() == 1) {
    return 1.0;
  }

  using ::std::find_if;

  auto singleton = find_if(m_afs.begin(), m_afs.end(),
                           [](AFS::value_type p)
                           {
                             return p.first.singleton();
                           });

  Allele allele = singleton->first;

  decltype(m_afs.deme()) deme = 0;
  while (allele[deme] == 0) {
    ++deme;
  }

  double dsize = static_cast<double>(m_afs.size(deme));

  AFS base = m_afs.remove(allele).add(allele.remove(deme));

  // probability of a sample excluding one of singleton alleles.
  double val = ESFProb(base, m_param, m_esf_prob_cache, m_hit_prob_cache).compute();

  for (auto a: base) {

    // add a gene to previously an observed allele.
    Allele na = a.first.add(deme);

    AFS other = base.remove(a.first).add(na);

    val -= ESFProb(other, m_param, m_esf_prob_cache, m_hit_prob_cache).compute() * \
        other[na] * na[deme] / dsize;


  }

  val *= dsize / m_afs[allele];

  return val;

}


double ESFProb::compute_coal_probs(ExitAFSData const& data, HitProb const& hp) {

  double val = 0.0;
  AFS afs = data.afs;
  State state = data.state;
  double factor = data.factor;

  auto deme = m_init.deme();

  for (auto a: afs) {

    auto allele = a.first;

    AFS other0 = afs.remove(allele);

    for (Index i = 0; i < deme; ++i) {

      if (allele[i] > 1) {

        Allele na = allele.remove(i);
        AFS other1 = other0.add(na);
        // double tmp = hp.get(state, i) * \
        //     ESFProb(other1, m_param, m_esf_prob_cache, m_hit_prob_cache).compute();
        // val += tmp * factor * na[i] * other1[na] / other1.size(i);

        std::cout << other1 << ": " << state << ": " << (double)na[i] * other1[na] / other1.size(i) << "\n";

      }

    }

  }

  return val;

}


}
