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

class ESFProbTest: public ::testing::Test {


 protected:

  ESFProbTest()
      : allele1({1,0}),
        afs1({allele1}),
        param({0,1,1,0},{1,1},{1,1}) {}

  Allele allele1;
  AFS afs1;
  Param param;

};


TEST_F(ESFProbTest, Base) {

  double exp = 1.0;
  auto val = ESFProb(afs1, param).compute();

  EXPECT_DOUBLE_EQ(exp, val);

}


}
