// -*- mode: c++; coding: utf-8; -*-

// test.cc - brief description

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

#include "util.hh"
#include "state_space.hh"

int main() {
  using std::cout;
  using std::endl;
  using esf::Index;
  using esf::Binomial;
  using esf::StateSpace;

  Index n = 3;

  // for (Index k = 0; k <= n; ++k) {
  //   cout << Binomial(n, k) << endl;
  // }

  // cout << "state: ";
  // for (auto s: StateSpace::IndexToState({1,2,3}, 8)) {
  //   cout << s << " ";
  // }
  // cout << endl;

  // cout << StateSpace::StateToIndex({1,2,3}, {1, 0, 0, 0, 2, 0, 0,
  // 0, 3}) << endl;

  for (auto s: StateSpace::Neighbors({1,2,3}, 8)) {
    cout << s << " ";
  }
  cout << endl;

  return 0;
};
