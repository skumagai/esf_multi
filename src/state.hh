// -*- mode: c++; coding: utf-8; -*-

// state.hh - State in structured ESF

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

#ifndef ESF_MULTI_STATE_HH
#define ESF_MULTI_STATE_HH


#include <vector>

#include "typedef.hh"
#include "init.hh"


namespace esf {

using ::std::ostream;
using ::std::vector;


// This class represents a state. States contain information of inital
// location of genes and the current location.
class State {

 public:

  typedef vector<esf_uint_t> value_type;

  typedef typename value_type::const_iterator const_iterator;

  private:

  Init m_init;

  value_type m_data;

  value_type compute_state(esf_uint_t);

  esf_uint_t compute_id() const;

  value_type expand_init() const;

 public:

  State() = default;

  State(State const&) = default;

  State(State&&) = default;

  State& operator=(State const&) = default;

  State& operator=(State&&) = default;

  ~State() = default;

  // Create a state from initial condition. This conversion assumes
  // that all genes remain in their initial locations.
  State(Init const&);

  // Createa i-th state associated with an initial condition.
  State(Init const&, esf_uint_t);

  // Create a state associated with an initial condition. The state is
  // specified as a list of number of genes in each deme grouped by
  // their initial location. That said, the list has length of n^2 for
  // n the number of demes.  Order of elements is {g_00, g_01,...g_10,
  // g_11,..., a_nn} for g_ij number of genes where i the inital
  // location and j the current location.
  State(Init const&, value_type const&);

  State move(esf_uint_t, esf_uint_t) const;

  // Computes a list of states, which are diferent from the current
  // state by one migration event.  This excludes the currentstate itself.
  vector<State> neighbors() const;

  // Returns an ID associated with the current state. The order of
  // states is stable, but forward and backward comptatibilities are
  // not guaranteed.
  esf_uint_t id() const;

  // Returns the number of demes.
  esf_uint_t deme() const;

  esf_uint_t operator[](esf_uint_t) const;

  const_iterator begin() const;

  const_iterator end() const;

  // Comparison operators.  Some of STL containers require object to
  // implement these.
  friend bool operator==(State const&, State const&);

  friend bool operator<(State const&, State const&);

  friend State operator+(State const&, State const&);

};


ostream& operator<<(ostream&, State const&);

}


#endif // ESF_MULTI_STATE_HH
