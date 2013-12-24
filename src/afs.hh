// -*- mode: c++; coding: utf-8; -*-

// afs.hh - Allele frequency spectrum

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

#ifndef ESF_MULTI_AFS_HH
#define ESF_MULTI_AFS_HH

#include <map>

#include "typedef.hh"
#include "allele.hh"
#include "enum.hh"

namespace esf {


class AFS {

 private:

  typedef typename ::std::map<Allele, Index> data_type;

 public:

  typedef typename data_type::key_type key_type;

  typedef typename data_type::mapped_type mapped_type;

  typedef typename data_type::value_type value_type;

  typedef typename data_type::iterator iterator;

  typedef typename data_type::const_iterator const_iterator;

 private:

  data_type data;

  AFS(AFS const&, Allele const&, Mode);

 public:

  AFS(::std::vector<Allele> const&);

  AFS(AFS const&) = default;

  AFS& operator=(AFS const&) = default;

  ~AFS() {};

  AFS add(Allele);

  AFS remove(Allele);

  bool singleton() const;

  Index operator[](const Allele&) const;

  Index size(Index) const;

  Index size() const;

  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

};


};


#endif // ESF_MULTI_AFS_HH
