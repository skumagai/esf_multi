// -*- mode: c++; coding: utf-8; -*-

// util.cc - a collection of (small) utilitiy functions

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

#include <functional>
#include <numeric>
#include <vector>

#include "typedef.hh"


using ::std::vector;


namespace esf {


using std::multiplies;
using std::partial_sum;


namespace {

vector<Index> mult_factors(vector<Index>);

};


Index binomial(Index n, Index k) {
  // compute in double then cast back to long
  Index value = 1;

  if (k > n - k) {
    k = n - k;
  }

  for (Index i = 0; i < k; ++i) {
    value *= (n - i);
    value /= (i + 1);
  }

  return value;
};


Index index_n_to_1(vector<Index> dim, vector<Index> idx) {
  auto nsize = dim.size();
  auto accum = mult_factors(dim);

  Index val = 0;
  for (Index i = 0; i < nsize; ++i) {
    val += idx[i] * accum[i];
  }

  return val;
};


vector<Index> index_1_to_n(vector<Index> dim, Index idx) {
  auto nsize = dim.size();
  auto accum = mult_factors(dim);

  vector<Index> vals(nsize);
  for (Index i = 0; i < nsize; ++i) {
    vals[i] = (idx / accum[i]) % dim[i];
  }

  return vals;
};


namespace {


vector<Index> mult_factors(vector<Index> dim) {
  vector<Index> accum(dim.size());
  accum[0] = 1;
  partial_sum(dim.begin(), dim.end() - 1, accum.begin() + 1, multiplies<Index>());
  return accum;
};


};


};
