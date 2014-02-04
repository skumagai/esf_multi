// -*- mode: c++; coding: utf-8; -*-

// allele_test.cc - [unit test] allele

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

#include "allele.hh"
#include "init.hh"
#include "state.hh"
#include "typedef.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <vector>

namespace {

using ::std::find;
using ::std::vector;

using ::esf::Allele;
using ::esf::esf_uint_t;
using ::esf::ExitAlleleData;
using ::esf::Init;
using ::esf::State;

class AlleleTest: public ::testing::Test {

 protected:

  AlleleTest(): s(Allele({1, 0, 0})),
                ns(Allele({1, 2, 0})) {}

  Allele s;  // singleton
  Allele ns; // non-singleton

};

TEST_F(AlleleTest, Singleton) {

  EXPECT_EQ(1, s.size());

}

TEST_F(AlleleTest, SingletonAdd) {

  auto a = s.add(1);
  auto b = s.add(0).add(0);

  EXPECT_EQ(1, a[0]);
  EXPECT_EQ(1, a[1]);
  EXPECT_EQ(0, a[2]);
  EXPECT_EQ(3, b[0]);
  EXPECT_EQ(0, b[1]);
  EXPECT_EQ(0, b[2]);

  EXPECT_EQ(2, a.size());
  EXPECT_EQ(3, b.size());

}


TEST_F(AlleleTest, SingletonRemove) {

  auto a = s.remove(0);

  for (auto i: a) {

    EXPECT_EQ(0, i);

  }

}


TEST_F(AlleleTest, NonSingleton) {

  EXPECT_EQ(3, ns.size());

}


TEST_F(AlleleTest, NonSingletonADD) {

  auto a = ns.add(2);

  EXPECT_EQ(1, a[0]);
  EXPECT_EQ(2, a[1]);
  EXPECT_EQ(1, a[2]);


}


TEST_F(AlleleTest, CheckSingleton) {

  EXPECT_EQ(true, s.singleton());
  EXPECT_NE(true, ns.singleton());

}


TEST_F(AlleleTest, AlleleCompare) {

  bool less = s < ns;

  EXPECT_EQ(true, less);
  // == and > are not implemented.

}


// TEST_F(AlleleTest, ReacheableAlleles) {

//   vector<ExitAlleleData> exp =
//       {
//         ExitAlleleData{Allele{{3, 0}}, State{Init{{2, 1}}, {2, 0, 1, 0}}, 1.0},
//         ExitAlleleData{Allele{{2, 1}}, State{Init{{2, 1}}, {2, 0, 0, 1}}, 1.0},
//         ExitAlleleData{Allele{{2, 1}}, State{Init{{2, 1}}, {1, 1, 1, 0}}, 1.0},
//         ExitAlleleData{Allele{{1, 2}}, State{Init{{2, 1}}, {1, 1, 0, 1}}, 1.0},
//         ExitAlleleData{Allele{{1, 2}}, State{Init{{2, 1}}, {0, 2, 1, 0}}, 1.0},
//         ExitAlleleData{Allele{{0, 3}}, State{Init{{2, 1}}, {0, 2, 0, 1}}, 1.0}
//       };


//   Allele allele({2, 1});

//   auto data = allele.reacheable();

//   auto crit = exp.end();

//   EXPECT_EQ(exp.size(), data.size());

//   for (auto datum: data) {

//     EXPECT_NE(crit, find(exp.begin(), exp.end(), datum));

//   }

// }


}  // namespace
