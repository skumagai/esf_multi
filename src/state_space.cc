// -*- mode: c++; coding: utf-8; -*-

// state_space.cc - State space

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

#include <iostream>
#include <algorithm>
#include <ctgmath>
#include <functional>
#include <numeric>

#include <boost/range/irange.hpp>

#include "state_space.hh"
#include "util.hh"

namespace esf {

using std::accumulate;
using std::copy;
using std::multiplies;
using std::partial_sum;
using std::sort;
using std::transform;
using boost::irange;

namespace {
State GenerateState(Index idx, Index ndeme, Index ngene);
Index GenerateIndex(State::iterator begin, State::iterator end);
};


Index StateSpace::StateToIndex(Init init, State state) {
  auto ndeme = init.size();
  Index sum = 0, mult = 1;
  State::iterator begin = state.begin();
  for (Index i = 0; i < ndeme; ++i, begin += ndeme) {
    sum += mult * GenerateIndex(begin, begin + ndeme);
    mult *= Binomial(i + ndeme, i + 1);
  }
  return sum;
}

State StateSpace::IndexToState(Init init, Index idx) {
  auto ndeme = init.size();
  IndexList dim;
  dim.resize(ndeme);
  transform(init.begin(), init.end() - 1, dim.begin(),
            [ndeme](Index i) {return Binomial(i + ndeme - 1, i);});
  *(dim.end() - 1) = 1;
  IndexList dim2;
  dim2.resize(ndeme);
  dim2[0] = 1;
  partial_sum(dim.begin(), dim.end() - 1, dim2.begin() + 1, multiplies<Index>());
  State state;
  State substate;
  for (Index i = 0; i < ndeme; ++i) {
    substate = GenerateState((idx / dim2[i]) % dim[i], ndeme, init[i]);
    state.insert(state.end(), substate.begin(), substate.end());
  };
  return state;
}

StateList StateSpace::Neighbors(Init init, Index idx) {
  // Based on incoming states rather than outgoing states.
  auto ndeme = init.size();
  auto deme = 0;
  State state = IndexToState(init, idx);
  State new_state(state.size());
  StateList neighbors;
  for (Index tar = 0; tar < state.size(); ++tar) {
    deme = tar / ndeme;
    for (Index src = ndeme * deme; src < ndeme * (deme + 1); ++src) {
      if (state[src] > 0 && src != tar) {
        copy(state.begin(), state.end(), new_state.begin());
        new_state[src] -= 1;
        new_state[tar] += 1;
        neighbors.push_back(new_state);
      }
    }
  }
  return neighbors;
}

Init StateSpace::StateToInit(State state) {
  auto ndeme = static_cast<Index>(sqrt(state.size()));
  Init init(ndeme);
  for (Index d = 0; d < ndeme; ++d) {
    for (Index a = 0; a < ndeme; ++a) {
      init[d] += state[d + a * ndeme];
    }
  }
  return init;
}

State StateSpace::InitToState(Init init) {
  auto ndeme = init.size();
  State state(ndeme * ndeme);
  for (Index d = 0; d < ndeme; ++d) {
    state[d * (1 + ndeme)] = init[d];
  }
  return state;
}

namespace {

State GenerateState(Index idx, Index ndeme, Index ngene) {
  if (ndeme == 1) {
    return {ngene};
  }

  IndexList offsets = {0};
  auto range = irange<Index>(0, ngene);
  IndexList tmp;
  tmp.resize(ngene);
  transform(range.begin(),
            range.end(),
            tmp.begin(),
            [ndeme, ngene](Index i) {return Binomial(ngene - i + ndeme - 2, ngene - i);});
  offsets.insert(offsets.end(), tmp.begin(), tmp.end());
  IndexList::iterator oiter = offsets.begin(), end = offsets.end();
  Index offset = 0;
  while (oiter < end && idx >= *oiter) {
    idx -= *oiter;
    offset += 1;
    oiter += 1;
  }
  if (offset > 0) {
    offset -= 1;
  }

  State state = {offset};
  auto rec = GenerateState(idx, ndeme - 1, ngene - offset);
  state.insert(state.end(), rec.begin(), rec.end());
  return state;
}

Index GenerateIndex(State::iterator begin, State::iterator end) {
  auto size = end - begin;
  if (size == 1) {
    return 0;
  }
  size -= 1;
  auto total = accumulate(begin, end, 0);
  auto range = irange<Index>(0, *begin);
  State tmp;
  tmp.resize(range.size());
  transform(range.begin(),
            range.end(),
            tmp.begin(),
            [total, size](Index i) {return Binomial(total - i + size - 1, total - i);});
  return accumulate(tmp.begin(), tmp.end(), 0) + GenerateIndex(begin + 1, end);
}

};

};
