// -*- mode: c++; coding: utf-8; -*-

// afs_test.cc - Unit tests for AFS class

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
#include "gtest/gtest.h"

#include <iostream>

namespace {

class AFSTest: public ::testing::Test {

 protected:

  AFSTest()
      : a100({1, 0, 0}),
        a010({0, 1, 0}),
        a200({2, 0, 0}),
        a020({0, 2, 0}),
        a002({0, 0, 2}),
        s({::esf::Allele({1,0,0}),
              ::esf::Allele({1,0,0}),
              ::esf::Allele({0,1,0}),
              ::esf::Allele({0,0,2})}),
        ns({::esf::Allele({2,0,0}),
                ::esf::Allele({0,2,0}),
                ::esf::Allele({0,0,2})}),
        c0({::esf::Allele({1,0,0}),
                ::esf::Allele({1,0,0})}),
        c1({::esf::Allele({1,0,0}),
                ::esf::Allele({1,0,0})}) {}


  ::esf::Allele a100, a010, a200, a020, a002;
  ::esf::AFS s;  // singleton only
  ::esf::AFS ns; // non singleton present
  ::esf::AFS c0, c1;  // for comparison

};


TEST_F(AFSTest, Singleton) {

  EXPECT_EQ(5, s.size());

  EXPECT_EQ(2, s.size(0));
  EXPECT_EQ(1, s.size(1));
  EXPECT_EQ(2, s.size(2));


  EXPECT_EQ(2, s[a100]);
  EXPECT_EQ(1, s[a010]);
  EXPECT_EQ(1, s[a002]);

}


TEST_F(AFSTest, SingletonAdd) {

  auto a0 = s.add(::esf::Allele({1, 0, 0}), 0);

  EXPECT_EQ(0, a0[::esf::Allele({1, 0, 0})]);
  EXPECT_EQ(1, a0[::esf::Allele({2, 0, 0})]);

}


TEST_F(AFSTest, SingletonRemove) {

  auto a0 = s.remove(::esf::Allele({1, 0, 0}), 0);

  EXPECT_EQ(0, a0[::esf::Allele({1, 0, 0})]);

}


TEST_F(AFSTest, NonSingleton) {

  EXPECT_EQ(6, ns.size());

  EXPECT_EQ(2, ns.size(0));
  EXPECT_EQ(2, ns.size(1));
  EXPECT_EQ(2, ns.size(2));

  EXPECT_EQ(1, ns[a200]);
  EXPECT_EQ(1, ns[a020]);
  EXPECT_EQ(1, ns[a002]);

}


TEST_F(AFSTest, NonSingletonAdd) {

  auto a0 = ns.add(::esf::Allele({1, 0, 0}), 0);

  EXPECT_EQ(0, a0[::esf::Allele({1, 0, 0})]);
  EXPECT_EQ(2, a0[::esf::Allele({2, 0, 0})]);

}


TEST_F(AFSTest, NonSingletonRemove) {

  auto a0 = ns.remove(::esf::Allele({2, 0, 0}), 0);

  EXPECT_EQ(1, a0[::esf::Allele({1, 0, 0})]);
  EXPECT_EQ(0, a0[::esf::Allele({2, 0, 0})]);

}


TEST_F(AFSTest, CheckSingleton) {

  EXPECT_EQ(true, s.singleton());
  EXPECT_NE(true, ns.singleton());

}



}  // namespace


int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();

}
