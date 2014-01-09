// -*- mode: c++; coding: utf-8; -*-

// util.hh - a collection of small utility functions

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

#ifndef ESF_MULTI_UTIL_HH
#define ESF_MULTI_UTIL_HH

#include <functional>
#include <numeric>
#include <type_traits>
#include <vector>


namespace esf {


using ::std::enable_if;
using ::std::is_integral;
using ::std::is_unsigned;
using ::std::is_signed;
using ::std::make_signed;
using ::std::make_unsigned;
using ::std::multiplies;
using ::std::partial_sum;
using ::std::vector;


// Convert unsigned integer types to corresponding signed types.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
typename make_signed<T>::type sign(T val);

// Convert singed integer types to corresponding unsinged types.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
typename make_unsigned<T>::type unsign(T val);

// Compute binomial coefficient.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
T binomial(T n, T k);

// Converts multidimensional index to 1 dimensional index. The most
// rapidly chaning element is the first element, and the most slowly
// chaning element is the last element.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
T index_n_to_1(vector<T> dim, vector<T> idx);


// Convert one dimensional index to multidimensional index.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
vector<T> index_1_to_n(vector<T> dim, T idx);


template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
vector<T> mult_factors(vector<T> dim);


// template function definitions

template <typename T, class>
typename make_signed<T>::type sign(T val) {

  typedef typename make_signed<T>::type sT;

  return static_cast<sT>(val);

}


template <typename T, class>
typename make_unsigned<T>::type unsign(T val) {

  typedef typename make_unsigned<T>::type uT;

  return static_cast<uT>(val);

}


template <typename T, class>
vector<T> mult_factors(vector<T> dim) {

  vector<T> accum(dim.size());

  accum[0] = 1;

  partial_sum(dim.begin(), dim.end() - 1, accum.begin() + 1, multiplies<T>());

  return accum;

}


template <typename T, class>
T binomial(T n, T k) {

  T value = 1;

  if (k > n - k && n - k >= 0) {

    k = n - k;

  }

  for (T i = 0; i < k; ++i) {

    value *= (n - i);
    value /= (i + 1);

  }

  return value;

}


// Converts multidimensional index to 1 dimensional index. The most
// rapidly chaning element is the first element, and the most slowly
// chaning element is the last element.
template <typename T, class>
T index_n_to_1(vector<T> dim, vector<T> idx) {

  typename vector<T>::size_type size = dim.size();

  auto accum = mult_factors(dim);

  T value = 0;

  for (decltype(size) i = 0; i < size; ++i) {

    value += idx[i] * accum[i];

  }

  return value;

}


// Convert one dimensional index to multidimensional index.
template <typename T, class>
vector<T> index_1_to_n(vector<T> dim, T idx) {

  typename vector<T>::size_type size = dim.size();

  auto accum = mult_factors(dim);

  vector<T> values(size);

  for (decltype(size) i = 0; i < size; ++i) {

    values[i] = (idx / accum[i]) % dim[i];

  }

  return values;

}


}


#endif // ESF_MULTI_UTIL_HH
