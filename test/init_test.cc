// -*- mode: c++; coding: utf-8; -*-

// init_test.cc - [unit test] init

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

#include "afs.hh"
#include "allele.hh"
#include "init.hh"
#include "state.hh"
#include "gtest/gtest.h"

namespace {

using ::esf::AFS;
using ::esf::Allele;
using ::esf::Init;
using ::esf::State;


class InitTest: public ::testing::Test {

 protected:

  InitTest()
      : i2({2, 3}), i21({2, 3}), i22({3, 2}), i3({2, 3, 2}) {}

  Init i2, i21, i22;  // 2-deme
  Init i3;  // 3-deme

};


TEST_F(InitTest, TwoDemePopDim) {

  EXPECT_EQ(2, i2[0]);
  EXPECT_EQ(3, i2[1]);

}


TEST_F(InitTest, TwoDemeDimPerDeme) {

  EXPECT_EQ(3, i2.dim(0));
  EXPECT_EQ(4, i2.dim(1));

}


TEST_F(InitTest, TwoDemeTotalDim) {

  EXPECT_EQ(12, i2.dim());

}


TEST_F(InitTest, TwoDemeNumber) {

  EXPECT_EQ(2, i2.deme());

}


TEST_F(InitTest, ThreeDemePopDim) {

  EXPECT_EQ(2, i3[0]);
  EXPECT_EQ(3, i3[1]);
  EXPECT_EQ(2, i3[2]);

}


TEST_F(InitTest, ThreeDemeDimPerDeme) {

  EXPECT_EQ(6, i3.dim(0));
  EXPECT_EQ(10, i3.dim(1));
  EXPECT_EQ(6, i3.dim(2));

}


TEST_F(InitTest, ThreeDemeTotalDim) {

  EXPECT_EQ(360, i3.dim());

}


TEST_F(InitTest, ThreeDemeNumber) {

  EXPECT_EQ(3, i3.deme());

}


TEST_F(InitTest, ConversionFromState) {

  State state(Init({2, 5}), {1, 1, 1, 4});
  Init init(state);

  EXPECT_EQ(init[0], 2);
  EXPECT_EQ(init[1], 5);

}


TEST_F(InitTest, ConversionFromAFS) {

  Allele a0({0,2}), a1({2,1}), a2({0,2});
  AFS afs({a0, a1, a2});
  Init init(afs);

  EXPECT_EQ(init[0], 2);
  EXPECT_EQ(init[1], 5);

}


TEST_F(InitTest, ConversionFromAllele) {
  Allele a{{3,2}};
  Init init{a};

  EXPECT_EQ(init[0], 3);
  EXPECT_EQ(init[1], 2);
}


TEST_F(InitTest, CompatibleComparison) {
  EXPECT_TRUE(i2 == i21);
  EXPECT_FALSE(i2 == i22);
}


TEST_F(InitTest, ImcompatibleComparision) {
  EXPECT_DEATH({bool cond = (i2 == i3);}, ".*");
}


TEST_F(InitTest, CompatibleAddition) {
  Init exp{{5, 5}};
  EXPECT_EQ(exp, i2 + i22);
}


TEST_F(InitTest, IncompatibleAddition) {
  EXPECT_DEATH({auto cond = i2 + i3;}, ".*");
}


}
