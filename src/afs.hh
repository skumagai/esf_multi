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
#include "state.hh"

namespace esf {


// Forward declarations

class Init;

// Internally used structs for bookkeeping purpose
struct ExitAllelePair;

struct ExitAFSPair;


// This class provides a representation of allele frequency spectrum
// (AFS).  It provides convenient access to underlying allele as well
// as creation of new AFS by modifying the current AFS.  This object
// is immutable.
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

  Index m_deme;

  AFS(AFS const&, Allele const&, Mode);

  ::std::vector<ExitAFSPair> build(::std::vector<Allele>,
                                   ::std::vector<Index>,
                                   data_type::const_iterator,
                                   data_type::const_iterator) const;

  ::std::vector<ExitAFSPair> sub_build(::std::vector<Allele>,
                                       ::std::vector<Index>,
                                       data_type::const_iterator,
                                       data_type::const_iterator,
                                       Index,
                                       ::std::vector<ExitAllelePair>::const_iterator,
                                       ::std::vector<ExitAllelePair>::const_iterator) const;

 public:

  AFS(::std::vector<Allele> const&);

  AFS(AFS const&) = default;

  AFS& operator=(AFS const&) = default;

  ~AFS() = default;

  // Create a new AFS by adding another allele to the current AFS.
  AFS add(Allele);

  // Create a new AFS by removing a preexisting allele from the
  // current AFS.
  AFS remove(Allele);

  // Test if the current AFS is singleton.  AFS is singleton if it
  // contains a single allele and if the allele is singleton.
  bool singleton() const;

  // Returns the multiplicity of the specified allele in AFS.
  // Different allele might have identical placement of genes.  THen
  // the multiplicty of allele is said to be more than one.
  Index operator[](Allele const&) const;

  // Return the number of genes in the specified deme.
  Index size(Index) const;

  // Return the total number of genes in AFS.
  Index size() const;

  // Returns the number of demes.
  Index deme() const;

  // Returns a list of AFS reacheable by one or more migrations.
  // Because origin and destination of genes need to be tracked,
  // return value is a list of pairs, whose first element is AFS and
  // the second element is its corresponding states.
  ::std::vector<ExitAFSPair> reacheable() const;

  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

  // Equality and less than operator are implemented to satisfy
  // requirement for storing in some of STL containers.
  friend bool operator==(AFS const&, AFS const&);

  friend bool operator<(AFS const&, AFS const&);

};


// This struct keeps AFS and its corresponding state together.
struct ExitAFSPair {

  AFS afs;

  State state;

};


bool operator==(ExitAFSPair const&, ExitAFSPair const&);

bool operator<(ExitAFSPair const&, ExitAFSPair const&);


};


#endif // ESF_MULTI_AFS_HH
