// -*- mode: c++; coding: utf-8; -*-

// init.hh - Initial state

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

#ifndef ESF_MULTI_INIT_HH
#define ESF_MULTI_INIT_HH


#include <vector>

#include "typedef.hh"

namespace esf {


class AFS;

class State;


class Init {

 public:

  typedef typename ::std::vector<Index> value_type;

  typedef typename value_type::iterator iterator;

  typedef typename value_type::const_iterator const_iterator;

 private:

  value_type m_data;

  Index m_deme;

  ::std::vector<Index> m_size;

  void set_size();

 public:

  Init() = default;

  Init(Init const&) = default;

  Init(Init&&) = default;

  Init& operator=(Init const&) = default;

  Init& operator=(Init&&) = default;

  Init(::std::vector<Index> const&);

  Init(State const&);

  Init(AFS const&);

  Index deme() const;

  Index size(Index) const;

  Index size() const;

  Index operator[](Index) const;

  friend bool operator==(Init const&, Init const&);

  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

};




}


#endif // ESF_MULTI_INIT_HH
