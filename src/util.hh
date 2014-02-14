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

template <typename InputIterator,
          class = typename enable_if<
            is_integral<typename InputIterator::value_type>::value>::type>
auto multinomial(InputIterator b, InputIterator e)
    -> typename InputIterator::value_type;

// Converts multidimensional index to 1 dimensional index. The most
// rapidly chaning element is the first element, and the most slowly
// chaning element is the last element.
template <typename InputIterator1,
          typename InputIterator2,
          class = typename enable_if<
            is_integral<typename InputIterator1::value_type>::value>::type>
auto index_n_to_1(InputIterator1 b1, InputIterator1 e1, InputIterator2 b2)
    -> typename InputIterator1::value_type;

// Convert one dimensional index to multidimensional index.
template <typename T,
          typename InputIterator,
          typename OutputIterator,
          class = typename enable_if<is_integral<T>::value>::type>
void index_1_to_n(InputIterator b, InputIterator e, OutputIterator o, T idx);

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

template <typename InputIterator1, typename InputIterator2>
bool next_k_selection(InputIterator1 b1, InputIterator1 e1,
                      InputIterator2 b2, InputIterator2 e2);

template <typename T,
          class = typename enable_if<is_integral<T>::value>::type>
class range {
 private:
  T m_data;
  T m_step;
 public:
  range(T init = 0, T step = 1) : m_data(init), m_step(step) {}
  T operator()() { T retval = m_data; m_data += m_step; return retval; }
};





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

template <typename InputIterator, class>
auto multinomial(InputIterator b, InputIterator e)
    -> typename InputIterator::value_type {
  typedef typename InputIterator::value_type T;
  T total = accumulate(b, e, static_cast<T>(0));
  T retval = 1;
  for (; b != e; ++b) {
    retval *= binomial(total, *b);
    total -= *b;
  }
  return retval;
}


template <typename InputIterator1, typename InputIterator2, class>
auto index_n_to_1(InputIterator1 b1, InputIterator1 e1, InputIterator2 b2)
    -> typename InputIterator1::value_type {
  typedef typename InputIterator1::value_type T;
  T len = unsign(distance(b1, e1));
  vector<T> accum(len);
  offsets(b1, e1, accum.begin());
  T value = 0;
  for (T i = 0; i < len; ++i, ++b2) {
    value += *b2 * accum[i];
  }
  return value;
}


template <typename T, typename InputIterator, typename OutputIterator, class>
void index_1_to_n(InputIterator b, InputIterator e, OutputIterator o, T idx) {
  T len = unsign(distance(b, e));
  vector<T> accum(len);
  offsets(b, e, accum.begin());
  for (T i = 0; i < len; ++i, ++b, ++o) {
    *o = (idx / accum[i]) % *b;
  }
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


template <typename InputIterator1, typename InputIterator2>
bool next_k_selection(InputIterator1 b1, InputIterator1 e1,
                      InputIterator2 b2, InputIterator2 e2) {

  InputIterator2 pos = e2 - 1;
  while (true) {
    auto iterr = find(b1, e1, *pos);
    auto iterl = find(b1, e1, *(pos - 1));
    if (pos == b2 and iterr == e1 - 1) {
      return false;
    }

    if (pos == b2 or iterr < iterl) {
      *pos = *(iterr + 1);
      for (auto b = pos + 1; b != e2; ++b) {
        *b = *b1;
      }
      return true;
    }

    if (iterr == iterl) {
      --pos;
    }
  }
}


}


#endif // ESF_MULTI_UTIL_HH
