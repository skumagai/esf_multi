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


// Forward declarations

// Internally used struct for bookkeeping allele and its corresponding state.
struct ExitAllelePair;


// This class represents an allele, a collection of genes with an
// identical genotype, and it keeps track of the number of genes
// present in each demes.
class Allele {

 private:

  typedef ::std::vector<Index> value_type;

 public:

  typedef value_type::iterator iterator;

  typedef value_type::const_iterator const_iterator;

 private:

  value_type m_data;

  Index m_total;

  Index m_deme;

  // This constructure creates a new object of Allele class by adding
  // or removing one gene to an already existing allele.  Because of
  // immutability of Allele class, new object has to be created every
  // time a gene is added to ore removed from pre-existing Allele class.
  Allele(Allele const&, Index, Mode);

  Allele(Allele&&, Index, Mode);

 public:

  // This constructor is designed to be invoked with data, which is
  // the locations and numbers of genes in present-day samples.
  Allele(::std::vector<Index> const&);

  Allele(Allele const&) = default;

  Allele(Allele&&) = default;

  Allele& operator=(Allele const&) = default;

  Allele& operator=(Allele&&) = default;

  // This function returns total number of genes in the allele.
  Index size() const;

  // Returns numbers of demes.
  Index deme() const;

  // Create a new allele by removing one gene at a specified deme in
  // the current allele.
  Allele remove(Index) const;

  // Create a new allele by adding one gene at a specified deme in the
  // current allele.
  Allele add(Index) const;

  // Returns true if this allele is singleton.  Singleton allele
  // contains only one gene.
  bool singleton() const;

  // Return a list of alleles reacheable from this allele by single
  // migration event.  One gene is taken from the present location,
  // and the same gene is placed to a new deme.  This function ignores
  // migrations where the source and target demes are the same.
  ::std::vector<ExitAllelePair> reacheable() const;

  // Exposes the iterator of underlying container.
  iterator begin();

  const_iterator begin() const;

  iterator end();

  const_iterator end() const;

  // Implements less-than operator for establishing strict weak
  // ordering.  This is required for storing objects of Allele in some
  // types of STL containers.
  bool operator<(Allele const&) const;

  // Returns number of genes in a specified deme.
  Index& operator[](Index);

  Index operator[](Index) const;

  // Implements equality test for the same reason as less-than operator.
  friend bool operator==(Allele const&, Allele const&);

};



// This struct serves as a temporaly storage holding together Allele object
// with its corresponding state, which keeps track of origin and
// current location of genes.
struct ExitAllelePair {

  Allele allele;

  ::std::vector<Index> state;

};


bool operator==(ExitAllelePair const&, ExitAllelePair const&);


}


#endif // ESF_MULTI_ALLELE_HH
