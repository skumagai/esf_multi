// -*- mode: c++; coding: utf-8; -*-

// allele.hh - Distribution of an allele across demes

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

#ifndef ESF_MULTI_ALLELE_HH
#define ESF_MULTI_ALLELE_HH

#include <vector>

#include "typedef.hh"
#include "enum.hh"

namespace esf {


struct ExitAllelePair;


class Allele {

 private:

  typedef ::std::vector<Index> value_type;

 public:

  typedef value_type::iterator iterator;

  typedef value_type::const_iterator const_iterator;

 private:

  value_type data;

  Index total;

  Index m_deme;

  Allele(Allele const&, Index, Mode);

  Allele(Allele&&, Index, Mode);

 public:

  Allele(::std::vector<Index> const&);

  Allele(Allele const&) = default;

  Allele(Allele&&) = default;

  Allele& operator=(Allele const&) = default;

  Allele& operator=(Allele&&) = default;

  Index size() const;

  Index deme() const;

  Allele remove(Index) const;

  Allele add(Index) const;

  bool singleton() const;

  ::std::vector<ExitAllelePair> reacheable() const;

  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

  bool operator<(Allele const&) const;

  Index& operator[](Index);

  Index const& operator[](Index) const;

  friend bool operator==(Allele const&, Allele const&);

};



struct ExitAllelePair {

  Allele allele;

  ::std::vector<Index> state;

};


bool operator==(ExitAllelePair const&, ExitAllelePair const&);


};


#endif // ESF_MULTI_ALLELE_HH
