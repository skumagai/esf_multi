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


class UtilTest: public ::testing::Test {

 protected:

  UtilTest() {}

};


TEST_F(UtilTest, BinomiaHandleZeroZero) {

  EXPECT_EQ(1, ::esf::binomial(0, 0));

}


TEST_F(UtilTest, BinomialHandleNonZeroZero) {

  EXPECT_EQ(1, ::esf::binomial(1, 0));
  EXPECT_EQ(1, ::esf::binomial(2, 0));
  EXPECT_EQ(1, ::esf::binomial(3, 0));

}


TEST_F(UtilTest, BinomialHandleZeroNonZero) {

  EXPECT_EQ(0, ::esf::binomial(0, 1));
  EXPECT_EQ(0, ::esf::binomial(0, 2));
  EXPECT_EQ(0, ::esf::binomial(0, 3));

}


TEST_F(UtilTest, BinomialHandleNormalCases) {

  EXPECT_EQ(3, ::esf::binomial(3, 1));
  EXPECT_EQ(3, ::esf::binomial(3, 2));
  EXPECT_EQ(1, ::esf::binomial(3, 3));

}


TEST_F(UtilTest, ConvertIndexNtoOne) {

  EXPECT_EQ(0, ::esf::index_n_to_1({2, 3}, {0, 0}));
  EXPECT_EQ(1, ::esf::index_n_to_1({2, 3}, {1, 0}));
  EXPECT_EQ(2, ::esf::index_n_to_1({2, 3}, {0, 1}));
  EXPECT_EQ(3, ::esf::index_n_to_1({2, 3}, {1, 1}));
  EXPECT_EQ(4, ::esf::index_n_to_1({2, 3}, {0, 2}));
  EXPECT_EQ(5, ::esf::index_n_to_1({2, 3}, {1, 2}));

  EXPECT_EQ(0, ::esf::index_n_to_1({2, 3, 2}, {0, 0, 0}));
  EXPECT_EQ(1, ::esf::index_n_to_1({2, 3, 2}, {1, 0, 0}));
  EXPECT_EQ(2, ::esf::index_n_to_1({2, 3, 2}, {0, 1, 0}));
  EXPECT_EQ(3, ::esf::index_n_to_1({2, 3, 2}, {1, 1, 0}));
  EXPECT_EQ(4, ::esf::index_n_to_1({2, 3, 2}, {0, 2, 0}));
  EXPECT_EQ(5, ::esf::index_n_to_1({2, 3, 2}, {1, 2, 0}));
  EXPECT_EQ(6, ::esf::index_n_to_1({2, 3, 2}, {0, 0, 1}));
  EXPECT_EQ(7, ::esf::index_n_to_1({2, 3, 2}, {1, 0, 1}));
  EXPECT_EQ(8, ::esf::index_n_to_1({2, 3, 2}, {0, 1, 1}));
  EXPECT_EQ(9, ::esf::index_n_to_1({2, 3, 2}, {1, 1, 1}));
  EXPECT_EQ(10, ::esf::index_n_to_1({2, 3, 2}, {0, 2, 1}));
  EXPECT_EQ(11, ::esf::index_n_to_1({2, 3, 2}, {1, 2, 1}));

}


TEST_F(UtilTest, ConvertIndexOneToN) {

  EXPECT_EQ(::std::vector<::esf::Index>({0, 0}), ::esf::index_1_to_n({2, 3}, 0));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 0}), ::esf::index_1_to_n({2, 3}, 1));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 1}), ::esf::index_1_to_n({2, 3}, 2));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 1}), ::esf::index_1_to_n({2, 3}, 3));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 2}), ::esf::index_1_to_n({2, 3}, 4));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 2}), ::esf::index_1_to_n({2, 3}, 5));

  EXPECT_EQ(::std::vector<::esf::Index>({0, 0, 0}), ::esf::index_1_to_n({2, 3, 2}, 0));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 0, 0}), ::esf::index_1_to_n({2, 3, 2}, 1));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 1, 0}), ::esf::index_1_to_n({2, 3, 2}, 2));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 1, 0}), ::esf::index_1_to_n({2, 3, 2}, 3));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 2, 0}), ::esf::index_1_to_n({2, 3, 2}, 4));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 2, 0}), ::esf::index_1_to_n({2, 3, 2}, 5));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 0, 1}), ::esf::index_1_to_n({2, 3, 2}, 6));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 0, 1}), ::esf::index_1_to_n({2, 3, 2}, 7));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 1, 1}), ::esf::index_1_to_n({2, 3, 2}, 8));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 1, 1}), ::esf::index_1_to_n({2, 3, 2}, 9));
  EXPECT_EQ(::std::vector<::esf::Index>({0, 2, 1}), ::esf::index_1_to_n({2, 3, 2}, 10));
  EXPECT_EQ(::std::vector<::esf::Index>({1, 2, 1}), ::esf::index_1_to_n({2, 3, 2}, 11));

}


}
