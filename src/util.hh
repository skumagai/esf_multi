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


using ::std::distance;
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
T index_n_to_1(vector<T> const& dim, vector<T> const& idx);


// Convert one dimensional index to multidimensional index.
template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
vector<T> index_1_to_n(vector<T> const& dim, T idx);


template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
vector<T> mult_factors(vector<T> const& dim);

template <typename InputIterator, typename OutputIterator>
void offsets(InputIterator b, InputIterator e, OutputIterator o);

template <typename T,
          typename OutputIterator,
          class = typename enable_if<is_integral<T>::value>::type>
void distribute_n(OutputIterator b, OutputIterator e, T idx, T k);

template <typename InputIterator, typename OutputIterator>
void dimensions(InputIterator b, InputIterator e, OutputIterator o);

template <typename T, typename InputIterator>
T convert_id(InputIterator b, InputIterator e);




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
vector<T> mult_factors(vector<T> const& dim) {

  vector<T> accum(dim.size());

  accum[0] = 1;

  partial_sum(dim.begin(), dim.end() - 1, accum.begin() + 1, multiplies<T>());

  return accum;

}


template <typename InputIterator, typename OutputIterator>
void offsets(InputIterator b, InputIterator e, OutputIterator o) {
  *o = 1;
  auto o0 = *o;
  partial_sum(b, e - 1, o + 1, multiplies<decltype(o0)>());
}


template <typename T, class>
T binomial(T n, T k) {

  T value = 1;

  if (k > n - k && n >= k) {

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
T index_n_to_1(vector<T> const& dim, vector<T> const& idx) {

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
vector<T> index_1_to_n(vector<T> const& dim, T idx) {

  typename vector<T>::size_type size = dim.size();

  auto accum = mult_factors(dim);

  vector<T> values(size);

  for (decltype(size) i = 0; i < size; ++i) {

    values[i] = (idx / accum[i]) % dim[i];

  }

  return values;

}


template <typename T, typename OutputIterator, class>
void distribute_n(OutputIterator b, OutputIterator e, T idx, T k) {
  auto len = unsign(distance(b, e));
  T offset = 0;
  for (T i = 0; i < len - 1; ++i) {
    T j = 0, idx_old = idx;
    while (idx <= idx_old &&
           (offset = binomial(k - j + len - i - 2, k - j)) <= idx) {
      idx_old = idx;
      idx -= offset;
      ++j;
    }
    k -= j;
    *b = j;
    ++b;
  }
  *b = k;
}


template <typename InputIterator, typename OutputIterator>
void dimensions(InputIterator b, InputIterator e, OutputIterator o) {
  auto o0 = *o;
  auto n = unsign(distance(b, e));
  transform(b, e, o, [n](decltype(o0) k) { return binomial(k + n - 1U, k); });
}


template <typename T, typename InputIterator>
T convert_id(InputIterator b, InputIterator e) {
  auto deme = unsign(distance(b, e));
  auto size = accumulate(b, e, 0U);

  T idx = 0;
  for (T i = 0; i < deme - 1; ++i, ++b) {
    for (T j = 0; j < *b; ++j) {
      idx += binomial(size - j + deme - i - 2, size - j);
    }
    size -= *b;
  }

  return idx;
}


}


#endif // ESF_MULTI_UTIL_HH
