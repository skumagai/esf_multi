// -*- mode: c++; coding: utf-8; -*-

// cache.hh - Implementation of cache

// Copyright (C) 2014 Seiji Kumagai

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

#ifndef ESF_MULTI_CACHE_HH
#define ESF_MULTI_CACHE_HH

#include <unordered_map>


namespace esf {


template <typename KEY, typename VALUE>
class Cache {

 public:

  typedef KEY key_type;
  typedef VALUE value_type;

 private:

  ::std::unordered_map<KEY, VALUE, std::hash<KEY>> m_cache;

  KEY const& m_root;

 public:

  Cache(KEY const&);

  KEY const& root() const;

  VALUE& at(KEY const&);
  VALUE at(KEY const&) const;

  VALUE& operator[](KEY const&);
  VALUE operator[](KEY const&) const;

};


template <typename KEY, typename VALUE>
Cache<KEY, VALUE>::Cache(KEY const& root)
    : m_root(root) {}


template <typename KEY, typename VALUE>
KEY const& Cache<KEY, VALUE>::root() const {

  return m_root;

}


template <typename KEY, typename VALUE>
VALUE& Cache<KEY, VALUE>::at(KEY const& key) {

  return m_cache.at(key);

}


template <typename KEY, typename VALUE>
VALUE Cache<KEY, VALUE>::at(KEY const& key) const {

  return m_cache.at(key);

}


template <typename KEY, typename VALUE>
VALUE& Cache<KEY, VALUE>::operator[](KEY const& key) {

  return m_cache[key];

}


template <typename KEY, typename VALUE>
VALUE Cache<KEY, VALUE>::operator[](KEY const& key) const {

  return m_cache[key];

}


}


#endif // ESF_MULTI_CACHE_HH
