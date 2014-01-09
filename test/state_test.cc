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


class StateTest: public ::testing::Test {

 protected:

  StateTest() {}

};


TEST_F(StateTest, TwoDemeID) {

  ::esf::Init init({2, 3});

  for (auto i = 0; i < init.size(); ++i) {

    ::esf::State s(init, i);

    ::std::vector<::esf::Index> v;

    for (auto j: s) {

      v.push_back(j);

    }

    ::esf::State t(init, v);

    EXPECT_EQ(i, t.id());

  }

}


TEST_F(StateTest, TwoDemeNeighbors) {

  ::esf::Init init({2, 3});

  ::esf::State s(init, 0);

  ::std::vector<::esf::Index> ids;

  for (auto adj: s.neighbors()) {

    ids.push_back(adj.id());

  }

  ::std::sort(ids.begin(), ids.end());

  ::std::vector<::std::vector<::esf::Index>> states =
      {
        {1, 1, 0 ,2},
        {0, 2, 1, 2}
      };

  ::std::vector<::esf::Index> test_data;

  for (auto ss: states) {

    ::esf::State stmp (init, ss);

    test_data.push_back(stmp.id());

  }

  ::std::sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

  ids.clear();

  ::esf::State s3(init,{1, 1, 1, 2});

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

    ::esf::State stmp(init, ss);

    test_data.push_back(stmp.id());

  }

  ::std::sort(ids.begin(), ids.end());
  ::std::sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

}


TEST_F(StateTest, TwoDemeConversionFromInit) {

  ::esf::Init orig({2, 3});

  ::esf::State state(orig, {1, 1, 2, 1});

  ::esf::Init test(state);

  ::std::vector<::esf::Index> data = {3, 2};

  for (auto i = 0; i < test.deme(); ++i) {

    EXPECT_EQ(data[i], test[i]);

  }

}


TEST_F(StateTest, ThreeDemeID) {

  ::esf::Init init({2, 3, 2});

  for (auto i = 0; i < init.size(); ++i) {

    ::esf::State s(init, i);

    ::std::vector<::esf::Index> v;

    for (auto j: s) {

      v.push_back(j);

    }

    ::esf::State t(init, v);

    EXPECT_EQ(i, t.id());

  }

}


TEST_F(StateTest, ThreeDemeNeighbors) {

  ::esf::Init init({2, 3, 2});

  ::esf::State s(init, 0);

  ::std::vector<::esf::Index> ids;

  for (auto adj: s.neighbors()) {

    ids.push_back(adj.id());

  }

  ::std::sort(ids.begin(), ids.end());

  ::std::vector<::std::vector<::esf::Index>> states =
      {
        {0, 1, 1, 0, 0, 3, 0, 0, 2},
        {1, 0 ,1, 0, 0, 3, 0, 0, 2},
        {0, 0, 2, 1, 0, 2, 0, 0, 2},
        {0, 0, 2, 0, 1, 2, 0, 0, 2},
        {0, 0, 2, 0, 0, 3, 1, 0, 1},
        {0, 0, 2, 0, 0, 3, 0, 1, 1}
      };

  ::std::vector<::esf::Index> test_data;

  for (auto ss: states) {

    ::esf::State stmp (init, ss);

    test_data.push_back(stmp.id());

  }

  ::std::sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

  ids.clear();

  ::esf::State s3(init, {1, 1, 0, 1, 1, 1, 1, 0, 1});

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

    ::esf::State stmp(init, ss);

    test_data.push_back(stmp.id());

  }

  ::std::sort(ids.begin(), ids.end());
  ::std::sort(test_data.begin(), test_data.end());

  EXPECT_EQ(ids.size(), test_data.size());

  for (auto i = 0; i < ids.size(); ++i) {

    EXPECT_EQ(ids[i], test_data[i]);

  }

}


TEST_F(StateTest, ThreeDemeConversionFromInit) {

  ::esf::Init orig({2, 3, 2});

  ::esf::State state(orig, {1, 0, 1, 0, 2, 1, 0, 1, 1});

  ::esf::Init test(state);

  ::std::vector<::esf::Index> data = {1, 3, 3};

  for (auto i = 0; i < test.deme(); ++i) {

    EXPECT_EQ(data[i], test[i]);

  }

}


}
