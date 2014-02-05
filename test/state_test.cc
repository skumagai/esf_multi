// -*- mode: c++; coding: utf-8; -*-

// state_test.cc - [unit test] state

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
#include <vector>
#include "init.hh"
#include "state.hh"
#include "typedef.hh"
#include "gtest/gtest.h"

namespace {

using ::std::sort;
using ::std::vector;

using ::esf::esf_uint_t;
using ::esf::Init;
using ::esf::State;


class StateTest: public ::testing::Test {

 protected:

  StateTest()
      : orig({2, 3, 2}), state(orig, {1, 0, 1, 0, 2, 1, 0, 1, 1}) {}

  Init orig;
  State state;

};


TEST_F(StateTest, TwoDemeID) {

  Init init({2, 3});

  for (auto i = 0; i < init.dim(); ++i) {

    State s(init, i);

    vector<esf_uint_t> v;

    for (auto j: s) {

      v.push_back(j);

    }

    State t(init, v);

    EXPECT_EQ(i, t.id());

  }

}


TEST_F(StateTest, TwoDemeNeighbors) {

  Init init({2, 3});

  State s(init, 0);

  vector<esf_uint_t> ids;

  for (auto adj: s.neighbors()) {

    ids.push_back(adj.id());

  }

  sort(ids.begin(), ids.end());

  vector<vector<esf_uint_t>> states =
      {
        {1, 1, 0 ,2},
        {0, 2, 1, 2}
      };

  vector<esf_uint_t> test_data;

  for (auto ss: states) {

    State stmp (init, ss);

    test_data.push_back(stmp.id());

  }

  sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

  ids.clear();

  State s3(init,{1, 1, 1, 2});

  for (auto adj: s3.neighbors()) {

    ids.push_back(adj.id());

  }

  states =
      {
        {2, 0, 1, 2},
        {0, 2, 1, 2},
        {1, 1, 0, 3},
        {1, 1, 2, 1}

      };

  test_data.clear();

  for (auto ss: states) {

    State stmp(init, ss);

    test_data.push_back(stmp.id());

  }

  sort(ids.begin(), ids.end());
  sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

}


TEST_F(StateTest, TwoDemeConversionFromInit) {

  Init orig({2, 3});

  State state(orig, {1, 1, 2, 1});

  Init test(state);

  vector<esf_uint_t> data = {3, 2};

  for (auto i = 0; i < test.deme(); ++i) {

    EXPECT_EQ(data[i], test[i]);

  }

}


TEST_F(StateTest, ThreeDemeID) {

  Init init({2, 3, 2});

  for (auto i = 0; i < init.dim(); ++i) {

    State s(init, i);

    vector<esf_uint_t> v;

    for (auto j: s) {

      v.push_back(j);

    }

    State t(init, v);

    EXPECT_EQ(i, t.id());

  }

}


TEST_F(StateTest, ThreeDemeNeighbors) {

  Init init({2, 3, 2});

  State s(init, 0);

  vector<esf_uint_t> ids;

  for (auto adj: s.neighbors()) {

    ids.push_back(adj.id());

  }

  sort(ids.begin(), ids.end());

  vector<vector<esf_uint_t>> states =
      {
        {0, 1, 1, 0, 0, 3, 0, 0, 2},
        {1, 0 ,1, 0, 0, 3, 0, 0, 2},
        {0, 0, 2, 1, 0, 2, 0, 0, 2},
        {0, 0, 2, 0, 1, 2, 0, 0, 2},
        {0, 0, 2, 0, 0, 3, 1, 0, 1},
        {0, 0, 2, 0, 0, 3, 0, 1, 1}
      };

  vector<esf_uint_t> test_data;

  for (auto ss: states) {

    State stmp (init, ss);

    test_data.push_back(stmp.id());

  }

  sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

  ids.clear();

  State s3(init, {1, 1, 0, 1, 1, 1, 1, 0, 1});

  for (auto adj: s3.neighbors()) {

    ids.push_back(adj.id());

  }

  states =
      {
        {0, 2, 0, 1, 1, 1, 1, 0, 1},
        {0, 1, 1, 1, 1, 1, 1, 0, 1},
        {2, 0, 0, 1, 1, 1, 1, 0, 1},
        {1, 0, 1, 1, 1, 1, 1, 0, 1},

        {1, 1, 0, 0, 2, 1, 1, 0, 1},
        {1, 1, 0, 0, 1, 2, 1, 0, 1},
        {1, 1, 0, 2, 0, 1, 1, 0, 1},
        {1, 1, 0, 1, 0, 2, 1, 0, 1},
        {1, 1, 0, 2, 1, 0, 1, 0, 1},
        {1, 1, 0, 1, 2, 0, 1, 0, 1},

        {1, 1, 0, 1, 1, 1, 0, 1, 1},
        {1, 1, 0, 1, 1, 1, 0, 0, 2},
        {1, 1, 0, 1, 1, 1, 2, 0, 0},
        {1, 1, 0, 1, 1, 1, 1, 1, 0}
      };

  test_data.clear();

  for (auto ss: states) {

    State stmp(init, ss);

    test_data.push_back(stmp.id());

  }

  sort(ids.begin(), ids.end());
  sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

}


TEST_F(StateTest, ThreeDemeConversionFromInit) {

  Init test(state);

  vector<esf_uint_t> data = {1, 3, 3};

  for (auto i = 0; i < test.deme(); ++i) {

    EXPECT_EQ(data[i], test[i]);

  }

}


TEST_F(StateTest, ValidComparison) {
  State exp{state};
  EXPECT_EQ(exp, state);

  // Only init is different.
  exp = State{Init{{1, 2, 3}}, {1, 0, 1, 0, 2, 1, 0, 1, 1}};
  EXPECT_FALSE(exp == state);

  // Only state vector is different.
  exp = State{Init{{2, 3, 2}}, {0, 1, 0, 3, 0, 0, 0, 1, 1}};
  EXPECT_FALSE(exp == state);

  // Both init and state vector are different.
  exp = State{Init{{1, 2, 3}}, {1, 0, 0, 2, 0, 0, 3, 0, 0}};
  EXPECT_FALSE(exp == state);
}


TEST_F(StateTest, ValidAddition) {
  State exp{Init{{4, 6, 4}}, {2, 0, 2, 0, 4, 2, 0, 2, 2}};
  EXPECT_EQ(exp, state + state);
}


}
