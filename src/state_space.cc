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

State generate_state(Index, Index, Index);

Index generate_index(State&, Index, Index);

};


Index StateSpace::state_to_index(Init init, State state) {
  auto ndeme = init.size();
  auto range = irange<Index>(0, ndeme * ndeme, ndeme);
  IndexList accum(ndeme);

  transform(range.begin(), range.end(), accum.begin(),
            [ndeme, &state](Index i)
            {
              return generate_index(state, i, i + ndeme);
            });

  return index_n_to_1(state_dim(init), accum);
}


State StateSpace::index_to_state(Init init, Index idx) {
  auto ndeme = init.size();

  auto idx_list = index_1_to_n(state_dim(init), idx);

  State state;
  State substate;
  for (Index i = 0; i < ndeme; ++i) {
    substate = generate_state(idx_list[i], ndeme, init[i]);
    state.insert(state.end(), substate.begin(), substate.end());
  };

  return state;
}


StateList StateSpace::neighbors(Init init, Index idx) {
  // List of nodes accessible from idx by a single emigration.
  auto ndeme = init.size();
  auto deme = 0;

  State state = index_to_state(init, idx);
  State new_state(state.size());

  StateList neighbors;
  for (Index src = 0; src < state.size(); ++src) {
    if (state[src] != 0) {
      deme = src / ndeme;

      for (Index tar = ndeme * deme; tar < ndeme * (deme + 1); ++tar) {
        if (src != tar) {
          copy(state.begin(), state.end(), new_state.begin());
          new_state[src] -= 1;
          new_state[tar] += 1;

          neighbors.push_back(new_state);
        }
      }
    }
  }

  return neighbors;
}


Init StateSpace::state_to_init(State state) {
  auto ndeme = static_cast<Index>(sqrt(state.size()));
  Init init(ndeme);

  for (Index d = 0; d < ndeme; ++d) {
    for (Index a = 0; a < ndeme; ++a) {
      init[d] += state[d + a * ndeme];
    }
  }

  return init;
}


State StateSpace::init_to_state(Init init) {
  auto ndeme = init.size();

  State state(ndeme * ndeme);
  for (Index d = 0; d < ndeme; ++d) {
    state[d * (1 + ndeme)] = init[d];
  }

  return state;
}


IndexList StateSpace::state_dim(Init init) {
  auto ndeme = init.size();
  IndexList dim(ndeme);

  transform(init.begin(), init.end(), dim.begin(),
            [ndeme](Index i)
            {
              return binomial(i + ndeme - 1, i);
            });

  return dim;
}


Index StateSpace::total_state(Init init) {
  auto dim = state_dim(init);
  return accumulate(dim.begin(), dim.end(), 1, multiplies<Index>());
};


namespace {


State generate_state(Index idx, Index ndeme, Index ngene) {

  Index j, offset = 0;
  State state(ndeme);

  for (Index i = 0; i < ndeme - 1; ++i) {

    j = 0;

    while (idx > 0 &&
           (offset = binomial(ngene - j + ndeme - i - 2, ngene - j)) <= idx) {

      idx -= offset;
      ++j;

    }

    ngene -= j;
    state[i] = j;

  }

  state[ndeme - 1] = ngene;

  return state;
}


Index generate_index(State& state, Index b, Index e) {

  Index idx = 0, curr = 0, ndeme = e - b;
  Index size = accumulate(state.begin() + b, state.begin() + e, 0);

  for (Index i = 0; i < ndeme - 1; ++i) {

    curr = state[b++];

    for (Index j = 0; j < curr; ++j) {

      idx += binomial(size - j + ndeme - i - 2, size - j);

    }

    size -= curr;

  }

  return idx;
}


};

};
