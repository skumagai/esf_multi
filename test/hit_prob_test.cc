// -*- mode: c++; coding: utf-8; -*-

// hit_prob_test.cc - [unit test] hit_prob

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
using ::std::vector;

#include "hit_prob.hh"
#include "init.hh"
#include "param.hh"
#include "state.hh"
#include "gtest/gtest.h"

namespace {

class HitProbTest: public ::testing::Test {

 protected:

  HitProbTest()
      : init2({2, 3}),
        init3({1, 2, 1}),
        param2({0.0, 1.0, 0.5, 0.0}, {1.0, 1.5}, {0.2, 0.4}),
        param3({0.0, 1.0, 0.5, 1.5, 0.0, 2.0, 2.5, 3.0, 0.0},
               {1.0, 1.5, 2.0}, {0.2, 0.4, 0.6}),
        hp2(init2, param2),
        hp3(init3, param3) {}

  ::esf::Init init2, init3;
  ::esf::Param param2, param3;
  ::esf::HitProb hp2, hp3;

};


TEST_F(HitProbTest, TwoDeme) {

  auto nstate = init2.dim();
  auto ndeme = init2.deme();

  vector<double> exp =
      {
        0., 0.00103878, 0., 0.0455599, 0.12159, 0.547156, 0., 0.000561967,
        0.00256614, 0.0115476, 0.080003, 0.0400015, 0.0000264864,
        0.000119189, 0.00185621, 0.000928104, 0.0200632, 0., 0.0000108476,
        5.4238e-6, 0.00028529, 0., 0.00144315, 0.
      };

  for (auto i = 0; i < nstate; ++i) {

    ::esf::State state(init2, i);

    for (auto j = 0; j < ndeme; ++j) {

      EXPECT_NEAR(exp[i * ndeme + j], hp2.get(i, j), 1.0e-6);
      EXPECT_NEAR(exp[i * ndeme + j], hp2.get(state, j), 1.0e-6);

    }

  }

}


TEST_F(HitProbTest, ThreeDeme) {


  auto nstate = init3.dim();
  auto ndeme = init3.deme();

  vector<double> exp =
      {
        0., 0., 0.019806, 0., 0., 0.0100659, 0., 0., 0.0359177, 0., 0.,
        0.0611395, 0., 0.0112255, 0.0149673, 0., 0., 0.0848855, 0.,
        0.0238654, 0.0318206, 0., 0.0403591, 0., 0., 0.169849, 0., 0., 0.,
        0.0127268, 0., 0., 0.00390996, 0.00532553, 0., 0.0106511, 0., 0.,
        0.0104442, 0., 0.00500293, 0., 0.0155173, 0., 0., 0.000520704, 0.,
        0.00104141, 0.000404654, 0., 0., 0.00273678, 0., 0., 0., 0.,
        0.00614411, 0., 0.00121538, 0.00162051, 0., 0., 0.00492505, 0.,
        0.00535094, 0.00713459, 0., 0.00908355, 0., 0., 0.0134821, 0., 0.,
        0.0127084, 0., 0., 0.0111906, 0., 0., 0.0437093, 0., 0., 0.,
        0.00243712, 0., 0.0012199, 0., 0.00190634, 0., 0., 0., 0.00246842,
        0., 0., 0.00344907, 0., 0.00281571, 0.00422357, 0., 0.000249283, 0.,
        0., 0.000132574, 0.000198862, 0., 0.000776759, 0., 0., 0., 0.,
        0.00272357, 0., 0., 0.000852867, 0.000847104, 0., 0.00169421, 0., 0.,
        0.00345647, 0., 0.00164827, 0., 0.00330404, 0., 0., 0., 0.00201258,
        0., 0., 0.00276916, 0., 0.00340535, 0.00510803, 0., 0.000457663, 0.,
        0.000915326, 0.000357125, 0., 0., 0.00171498, 0., 0., 0.000687756,
        0., 0., 0.000352559, 0.000528838, 0., 0.00275034, 0., 0.,
        0.000241668, 0., 0., 0.00014761, 0., 0., 0.000418154, 0., 0.
      };

  for (auto i = 0; i < nstate; ++i) {

    ::esf::State state(init3, i);

    for (auto j = 0; j < ndeme; ++j) {

      EXPECT_NEAR(exp[i * ndeme + j], hp3.get(i, j), 1.0e-6);
      EXPECT_NEAR(exp[i * ndeme + j], hp3.get(state, j), 1.0e-6);

    }

  }

}


}
