// -*- mode: c++; coding: utf-8; -*-

// param_test.cc - [unit test] param

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

#include "param.hh"
#include "gtest/gtest.h"

namespace {

class ParamTest: public ::testing::Test {

 protected:

  ParamTest()
      : p2d({0.0, 1.0, 2.0, 0.0}, {1.0, 1.0}, {1., 1.}),
        p3d({0.0, 0.25, 0.5, 0.75, 0.0, 1.0, 1.25, 1.5, 0.0}, {1., 1., 1.}, {1., 1., 1.}) {}

  ::esf::Param p2d;  // 2-deme
  ::esf::Param p3d;  // 3-deme

};


TEST_F(ParamTest, TwoDemeMigration) {

  EXPECT_DOUBLE_EQ(0.0, p2d.mig_rate(0, 0));
  EXPECT_DOUBLE_EQ(1.0, p2d.mig_rate(1, 0));
  EXPECT_DOUBLE_EQ(2.0, p2d.mig_rate(0, 1));
  EXPECT_DOUBLE_EQ(0.0, p2d.mig_rate(1, 1));

}

TEST_F(ParamTest, ThreeDemeMigration) {

  EXPECT_DOUBLE_EQ(0.0, p3d.mig_rate(0, 0));
  EXPECT_DOUBLE_EQ(0.25, p3d.mig_rate(1, 0));
  EXPECT_DOUBLE_EQ(0.5, p3d.mig_rate(2, 0));
  EXPECT_DOUBLE_EQ(0.75, p3d.mig_rate(0, 1));
  EXPECT_DOUBLE_EQ(0.0, p3d.mig_rate(1, 1));
  EXPECT_DOUBLE_EQ(1.0, p3d.mig_rate(2, 1));
  EXPECT_DOUBLE_EQ(1.25, p3d.mig_rate(0, 2));
  EXPECT_DOUBLE_EQ(1.5, p3d.mig_rate(1, 2));
  EXPECT_DOUBLE_EQ(0.0, p3d.mig_rate(2, 2));

}


}
