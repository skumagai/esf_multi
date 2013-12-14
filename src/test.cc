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

#include <cassert>
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
  using esf::binomial;

  Index n = 3;

  // Test conversion of state to index and vice versa.
  cout << "Check conversion of state to index vice versa: ";
  for (int i = 0; i < 180; ++i) {
    auto state = StateSpace::index_to_state({1,2,3}, i);
    assert (i == StateSpace::state_to_index({1,2,3}, state));
  }
  cout << "passed\n";

  return 0;
};
