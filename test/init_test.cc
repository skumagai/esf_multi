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


class InitTest: public ::testing::Test {

 protected:

  InitTest()
      : i2({2, 3}), i3({2, 3, 2}) {}

  ::esf::Init i2;  // 2-deme
  ::esf::Init i3;  // 3-deme

};


TEST_F(InitTest, TwoDemePopSize) {

  EXPECT_EQ(2, i2[0]);
  EXPECT_EQ(3, i2[1]);

}


TEST_F(InitTest, TwoDemeSizePerDeme) {

  EXPECT_EQ(3, i2.size(0));
  EXPECT_EQ(4, i2.size(1));

}


TEST_F(InitTest, TwoDemeTotalSize) {

  EXPECT_EQ(12, i2.size());

}


TEST_F(InitTest, TwoDemeNumber) {

  EXPECT_EQ(2, i2.deme());

}


TEST_F(InitTest, ThreeDemePopSize) {

  EXPECT_EQ(2, i3[0]);
  EXPECT_EQ(3, i3[1]);
  EXPECT_EQ(2, i3[2]);

}


TEST_F(InitTest, ThreeDemeSizePerDeme) {

  EXPECT_EQ(6, i3.size(0));
  EXPECT_EQ(10, i3.size(1));
  EXPECT_EQ(6, i3.size(2));

}


TEST_F(InitTest, ThreeDemeTotalSize) {

  EXPECT_EQ(360, i3.size());

}


TEST_F(InitTest, ThreeDemeNumber) {

  EXPECT_EQ(3, i3.deme());

}


TEST_F(InitTest, ConversionFromState) {

  ::esf::State state(::esf::Init({2, 5}), {1, 1, 1, 4});
  ::esf::Init init(state);

  EXPECT_EQ(init[0], 2);
  EXPECT_EQ(init[1], 5);

}


TEST_F(InitTest, ConversionFromAFS) {

  ::esf::Allele a0({0,2}), a1({2,1}), a2({0,2});
  ::esf::AFS afs({a0, a1, a2});
  ::esf::Init init(afs);

  EXPECT_EQ(init[0], 2);
  EXPECT_EQ(init[1], 5);

}


}
