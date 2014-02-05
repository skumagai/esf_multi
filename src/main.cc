// -*- mode: c++; coding: utf-8; -*-

// main.cc - brief description

// Copyright (C) 2014 Seiji Kumagai

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

#include <cstdlib>
#include <iostream>
#include <vector>

#include "afs.hh"
#include "allele.hh"
#include "esf_prob.hh"
#include "param.hh"

#include "cache.hh"
#include "hit_prob.hh"
#include "init.hh"

using ::std::strtod;
using ::std::strtoul;
using ::std::cerr;
using ::std::cout;
using ::std::strtod;
using ::std::vector;

using ::esf::AFS;
using ::esf::Allele;
using ::esf::ESFProb;
using ::esf::Param;

int main(int argc, char *argv[]) {

  if (argc < 6) {
    cerr << "usage: esf ndeme mig_rate pop_size mut_rate afs...\n";
    return 1;
  }

  auto ndeme = strtoul(argv[1], nullptr, 10);
  if (ndeme < 1) {
    cerr << "require at least one deme.\n";
    return 1;
  }

  if ((unsigned long)(argc - 5) % ndeme != 0) {
    cerr << "syntax error in allele frequency spectrum.\n";
    return 1;
  }

  vector<double> mig(ndeme * ndeme, strtod(argv[2], nullptr));
  for (decltype(ndeme) i = 0; i < ndeme; ++i) {
    mig[i * ndeme + i] = 0.0;
  }
  vector<double> pop(ndeme, strtod(argv[3], nullptr));
  vector<double> mut(ndeme, strtod(argv[4], nullptr));

  Param param = {mig, pop, mut};

  long demes = (long) ndeme;
  vector<Allele> alleles;

  for (auto i = 5; i < argc; i += ndeme) {
    vector<::esf::esf_uint_t> avec;
    for (auto j = 0; j < demes; ++j) {
      avec.push_back(strtoul(argv[i + j], nullptr, 10));
    }
    alleles.push_back(Allele{avec});
  }

  // for (auto a: alleles[0].reacheable()) {
  //   std::cout << a.allele << '\t' << a.state << '\t' << a.factor << '\n';
  // }

  AFS afs(alleles);

  ESFProb ep(afs, param);

  cout.precision(9);

  cout << ep.compute() << '\n';

  // for (auto a: afs.reacheable()) {
  //   cout << "factor: " << a.factor << "\t";
  //   cout << a.state << "\t";
  //   cout << a.afs << "\n";
  // }


  return 0;
}
