// -*- mode: c++; coding: utf-8; -*-

// mat_skel.hh - Skelton of matrix to speed up repeated computation of
// hitting probability

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

#ifndef ESF_MULTI_MAT_SKEL_HH
#define ESF_MULTI_MAT_SKEL_HH

#include "typedef.hh"

namespace esf {


class TransMatSkel {
 private:
  Index const row_idx;

  Index const col_idx;

  Value const coeff;

  Value* const val_ptr;

 public:
  TransMatSkel(Index const r, Index const c, Value const v, Value* const vp)
      : row_idx(r), col_idx(c), coeff(v), val_ptr(vp) {};

  Index column() const;

  Index row() const;

  Value value() const;
 };


class AbsMatSkel {
 private:
  Index const idx;

  Init const init;

  ValueList* const pop_size;

  ValueList* const mut_rate;

 public:
  AbsMatSkel(Index const i, Init const j, ValueList* const p, ValueList* m)
      : idx(i), init(j), pop_size(p), mut_rate(m) {};

  Index column() const;

  Index row() const;

  Value value() const;
};


};


#endif // ESF_MULTI_MAT_SKEL_HH
