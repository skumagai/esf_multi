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
#include "init.hh"
#include "state.hh"
#include "gtest/gtest.h"

#include <algorithm>
#include <vector>

namespace {

class AFSTest: public ::testing::Test {

 protected:

  AFSTest()
      : a100({1, 0, 0}),
        a010({0, 1, 0}),
        a200({2, 0, 0}),
        a020({0, 2, 0}),
        a002({0, 0, 2}),
        vs({a100, a100, a010, a002}),
        vns({a200, a020, a002}),
        vc0({a100, a100}),
        vc1({a100, a100}),
        s(vs),
        ns(vns),
        c0(vc0),
        c1(vc1) {}

  ::esf::Allele a100, a010, a200, a020, a002;
  ::std::vector<::esf::Allele> vs, vns, vc0, vc1;
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

  auto afs0 = s.add(a100);

  EXPECT_EQ(3, afs0[a100]);
  EXPECT_EQ(0, afs0[a200]);

}


TEST_F(AFSTest, SingletonRemove) {

  auto afs0 = s.remove(a100);
  auto afs1 = s.remove(a100).remove(a100);

  EXPECT_EQ(2, s[a100]);
  EXPECT_EQ(1, afs0[a100]);
  EXPECT_EQ(0, afs1[a100]);

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

  auto a0 = ns.add(a200);

  EXPECT_EQ(0, a0[a100]);
  EXPECT_EQ(2, a0[a200]);

}


TEST_F(AFSTest, NonSingletonRemove) {

  auto afs0 = ns.remove(a100);
  auto afs1 = ns.remove(a200);

  EXPECT_EQ(0, afs0[a100]);
  EXPECT_EQ(1, afs0[a200]);
  EXPECT_EQ(0, afs1[a100]);
  EXPECT_EQ(0, afs1[a200]);

}


TEST_F(AFSTest, CheckSingleton) {

  EXPECT_EQ(true, s.singleton());
  EXPECT_NE(true, ns.singleton());

}


TEST_F(AFSTest, AllReacheable) {

  using ::std::find;
  using ::std::vector;
  using ::esf::Allele;
  using ::esf::AFS;
  using ::esf::Init;
  using ::esf::State;
  using ::esf::ExitAFSPair;

  vector<Allele> v0 = {Allele({1, 0}), Allele({1, 0})};
  vector<Allele> v1 = {Allele({1, 0}), Allele({0, 1})};
  vector<Allele> v2 = {Allele({0, 1}), Allele({0, 1})};

  auto test = AFS(v1).reacheable();

  using ::std::vector;

  Init init({1, 1});

  vector<ExitAFSPair> exp =
      {
        ExitAFSPair({AFS(v0), State(init, {1, 0, 1, 0})}),
        ExitAFSPair({AFS(v1), State(init, {1, 0, 0, 1})}),
        ExitAFSPair({AFS(v1), State(init, {0, 1, 1, 0})}),
        ExitAFSPair({AFS(v2), State(init, {0, 1, 0, 1})})
      };

  EXPECT_EQ(exp.size(), test.size());

  auto end = exp.end();

  for (auto val: test) {

    EXPECT_NE(end, find(exp.begin(), exp.end(), val));

  }

  init = Init({2, 0});

  exp =
      {
        ExitAFSPair({AFS(v0), State(init, {2, 0, 0, 0})}),
        ExitAFSPair({AFS(v1), State(init, {1, 1, 0, 0})}),
        ExitAFSPair({AFS(v2), State(init, {0, 2, 0, 0})})
      };

  test = AFS(v0).reacheable();

  EXPECT_EQ(exp.size(), test.size());

  end = exp.end();

  for (auto val: test) {

    EXPECT_NE(end, find(exp.begin(), exp.end(), val));

  }

}


}  // namespace
