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
#include <iterator>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
using std::next;

#include "util.hh"
#include "state_space.hh"

int main() {
  using std::cout;
  using std::endl;
  using esf::Index;
  using namespace esf;
  // using namespace esf::StateSpace;
  using esf::binomial;

  Index n = 3;

  // for (Index k = 0; k <= n; ++k) {
  //   cout << binomial(n, k) << endl;
  // }

  // cout << "state:";
  // for (auto s: StateSpace::index_to_state({1,2,3}, 4)) {
  //   cout << ' ' << s;
  // }
  // cout << endl;

  // for (int i = 0; i < 10; ++i) {
  //   cout << "state(" << i << "):";
  //   for (auto s: StateSpace::index_to_state({1,2,3}, i)) {
  //     cout << ' ' << s;
  //   }
  //   cout << endl;
  // }


  cout << "state(1, 0, 0, 0, 2, 0, 0, 0, 3): ";
  cout << StateSpace::state_to_index({1,2,3}, {1, 0, 0, 0, 2, 0, 0,
  0, 3}) << endl;

  // for (auto s: StateSpace::neighbors({1,2,3}, 8)) {
  //   cout << "neighbors:";
  //   for (auto t: s) {
  //     cout << ' ' << t;
  //   }
  //   cout << endl;
  // }
  // cout << endl;

  // esf::IndexList dim = {2, 3, 4};
  // esf::IndexList idx = {1, 1, 0};
  // cout << esf::index_n_to_1(dim, idx) << endl;

  // for (auto a: esf::index_1_to_n(dim, 3)) {
  //   cout << a << endl;
  // }

  // cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

  // vector<int> a = {0, 1, 2, 3, 4, 5, 6, 7};
  // auto itr = next(a.begin(), 3);

  // cout << *itr << endl;
  // itr += 1;
  // cout << *itr << endl;


  return 0;
};
