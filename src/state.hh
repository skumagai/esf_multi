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

namespace esf {


class Init;


class State {

 public:

  typedef typename ::std::vector<Index> value_type;

  typedef typename value_type::iterator iterator;

  typedef typename value_type::const_iterator const_iterator;

  private:

  ::std::vector<Index> m_data;

  Index m_id;

  Index m_deme;

  Init *m_init;

  void compute_state();

  void compute_id();

  void expand_init();

 public:

  State();

  State(State const&);

  State(State&&) = default;

  State& operator=(State const&);

  State& operator=(State&&);

  ~State();

  State(Init const&, ::std::vector<Index> const&);

  State(Init const&);

  State(Init const&, Index);

  State(Init const&, ::std::vector<Index> const&);

  ::std::vector<State> neighbors();

  Index id() const;

  Index deme() const;

};


};


#endif // ESF_MULTI_STATE_HH
