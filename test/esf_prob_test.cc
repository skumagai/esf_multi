// -*- mode: c++; coding: utf-8; -*-

// esf_prob_test.cc - [unit test] ESFProb

// Copyright (C) 2014 Seiji Kumagai

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

#include "afs.hh"
#include "allele.hh"
#include "esf_prob.hh"
#include "param.hh"
#include "gtest/gtest.h"

namespace {

using ::esf::Allele;
using ::esf::AFS;
using ::esf::Param;
using ::esf::ESFProb;
using vector = ::std::vector<Allele>;


class ESFProbTest: public ::testing::Test {

 protected:

  ESFProbTest()
      : s20(vector({Allele({1, 0}), Allele({0, 1})})),
        s21(vector({Allele({1, 0}), Allele({1, 0})})),
        s30(vector({Allele({1, 0}), Allele({0, 1}), Allele({1, 0})})),
        s31(vector({Allele({1, 0}), Allele({1, 0}), Allele({1, 0})})),
        s32(vector({Allele({2, 0}), Allele({1, 0})})),
        s33(vector({Allele({2, 0}), Allele({0, 1})})),
        ns2(vector({Allele({2, 0})})),
        ns3(vector({Allele({3, 0})})),
        p({0.0, 1.0, 0.5, 0.0}, {1.0, 1.5}, {0.2, 0.4}) {}

  AFS s20, s21, s30, s31, s32, s33, ns2, ns3;
  Param p;

};


TEST_F(ESFProbTest, BaseCase) {

  ESFProb prob(AFS({Allele({1,0})}), p);

  EXPECT_DOUBLE_EQ(1., prob.compute());

}


TEST_F(ESFProbTest, WithSingleton) {

  ESFProb p0(s20, p);
  EXPECT_NEAR(0.4815596672047684, p0.compute(), 1e-6);

  ESFProb p1(s21, p);
  EXPECT_NEAR(0.25928225506022595, p1.compute(), 1e-6);

  ESFProb p2(s30, p);
  EXPECT_NEAR(0.09009334137786906, p2.compute(), 1e-6);

  ESFProb p3(s31, p);
  EXPECT_NEAR(0.03589004507443738, p3.compute(), 1e-6);

  ESFProb p4(s32, p);
  EXPECT_NEAR(0.33508831497868286, p4.compute(), 1e-6);

  ESFProb p5(s33, p);
  EXPECT_NEAR(0.3068718689857209, p5.compute(), 1e-6);

}


TEST_F(ESFProbTest, WithoutSingleton) {

  ESFProb p0(ns2, p);
  EXPECT_NEAR(0.740717744939774, p0.compute(), 1e-6);

  ESFProb p1(ns3, p);
  EXPECT_NEAR(0.6290216399468798, p1.compute(), 1e-6);

}


}
