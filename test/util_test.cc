// -*- mode: c++; coding: utf-8; -*-

// util_test.cc - [unit test] util

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


#include <vector>

#include "typedef.hh"
#include "util.hh"
#include "gtest/gtest.h"

namespace {

using ::std::vector;

using ::esf::binomial;
using ::esf::esf_uint_t;
using ::esf::index_n_to_1;
using ::esf::index_1_to_n;


class UtilTest: public ::testing::Test {

 protected:

  UtilTest() {}

};


TEST_F(UtilTest, BinomiaHandleZeroZero) {

  EXPECT_EQ(1, binomial(0, 0));

}


TEST_F(UtilTest, BinomialHandleNonZeroZero) {

  EXPECT_EQ(1, binomial(1, 0));
  EXPECT_EQ(1, binomial(2, 0));
  EXPECT_EQ(1, binomial(3, 0));

}


TEST_F(UtilTest, BinomialHandleZeroNonZero) {

  EXPECT_EQ(0, binomial(0, 1));
  EXPECT_EQ(0, binomial(0, 2));
  EXPECT_EQ(0, binomial(0, 3));

}


TEST_F(UtilTest, BinomialHandleNormalCases) {

  EXPECT_EQ(3, binomial(3, 1));
  EXPECT_EQ(3, binomial(3, 2));
  EXPECT_EQ(1, binomial(3, 3));

}


TEST_F(UtilTest, Convertesf_uint_tNtoOne) {

  vector<esf_uint_t> i{2, 3};
  vector<vector<esf_uint_t>> o = {{0, 0},
                                  {1, 0},
                                  {0, 1},
                                  {1, 1},
                                  {0, 2},
                                  {1, 2}};

  for (esf_uint_t j = 0; j < 6; ++j) {
    EXPECT_EQ(j, index_n_to_1(i.begin(), i.end(), o[j].begin()));
  }

  i = {2, 3, 2};
  o = {{0, 0, 0},
       {1, 0, 0},
       {0, 1, 0},
       {1, 1, 0},
       {0, 2, 0},
       {1, 2, 0},
       {0, 0, 1},
       {1, 0, 1},
       {0, 1, 1},
       {1, 1, 1},
       {0, 2, 1},
       {1, 2, 1}};
  for (esf_uint_t j = 0; j < 12; ++j) {
    EXPECT_EQ(j, index_n_to_1(i.begin(), i.end(), o[j].begin()));
  }
}


TEST_F(UtilTest, Convertesf_uint_tOneToN) {

  vector<esf_uint_t> i1 = {2, 3};
  vector<vector<esf_uint_t>> o = {{0, 0},
                                  {1, 0},
                                  {0, 1},
                                  {1, 1},
                                  {0, 2},
                                  {1, 2}};
  vector<esf_uint_t> i2(2);
  for (esf_uint_t j = 0; j < 6; ++j) {
    index_1_to_n(i1.begin(), i1.end(), i2.begin(), j);
    EXPECT_EQ(o[j], i2);
  }

  i1 = {2, 3, 2};
  i2.resize(3);
  o = {{0, 0, 0},
       {1, 0, 0},
       {0, 1, 0},
       {1, 1, 0},
       {0, 2, 0},
       {1, 2, 0},
       {0, 0, 1},
       {1, 0, 1},
       {0, 1, 1},
       {1, 1, 1},
       {0, 2, 1},
       {1, 2, 1}};

  for (esf_uint_t j = 0; j < 6; ++j) {
    index_1_to_n(i1.begin(), i1.end(), i2.begin(), j);
    EXPECT_EQ(o[j], i2);
  }
}


}
