// -*- mode: c++; coding: utf-8; -*-

// afs.cc - Implementation of allele frequency spectrum

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
#include <numeric>
#include <stdexcept>
#include <vector>

#include "afs.hh"
#include "allele.hh"

namespace esf {


AFS::AFS(AFS const& afs, Allele const& allele, Mode mode)
    : data(afs.data) {

  if (allele.size() > 0) {

    if (mode == Mode::ADD) {

      data[allele] += 1;

    } else if (data[allele] > 1) {

      data[allele] -= 1;

    } else {

      data.erase(allele);

    }

  }

}


AFS::AFS(::std::vector<Allele> const& avec) {

  for (auto allele: avec) {

    data[allele] += 1;

  }

}


AFS AFS::add(Allele allele) {

  return AFS(*this, allele, Mode::ADD);

}


AFS AFS::remove(Allele allele) {

  return AFS(*this, allele, Mode::REMOVE);

}


bool AFS::singleton() const {

  return ::std::any_of(data.begin(), data.end(),
                       [](value_type p)
                       {
                         return p.first.singleton();
                       }
                       );

}


Index AFS::operator[](const Allele& allele) const {

  try {

    return data.at(allele);

  } catch (::std::out_of_range& e) {

    return 0;

  }

}


Index AFS::size(Index deme) const {

  return ::std::accumulate(data.begin(), data.end(), 0,
                           [deme](Index a, value_type p)
                           {
                             return a + (p.first)[deme] * p.second;
                           }
                           );

}


Index AFS::size() const {

  return ::std::accumulate(this->begin(), this->end(), 0,
                           [](Index a, value_type p)
                           {
                             return a + (p.first).size() * p.second;
                           }
                           );
}


AFS::iterator AFS::begin() {

  return data.begin();

}


AFS::const_iterator AFS::begin() const {

  return data.begin();

}


AFS::iterator AFS::end() {

  return data.end();

}


AFS::const_iterator AFS::end() const {

  return data.end();

}


};
