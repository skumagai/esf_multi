# -*- mode: python; coding: utf-8; -*-

# afs2.py - brief description

# Copyright (C) 2014 Seiji Kumagai

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice (including the next
# paragraph) shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

def comb(n, k):

    value = 1

    if k > n - k and n - k >= 0:
        k = n - k

    for i in range(k):
        value *= (n - i)
        value /= (i + 1)

    return value


def assign(count, template, idx):

    if idx == len(template) - 1:
        base = list(template[0])
        base[idx] = count

        return [(base, template[1])]

    holder = []

    for i in range(count + 1):

        base = list(template[0])
        base[idx] = i

        holder.extend(assign(count - i, (base, template[1] * comb(count, i)), idx + 1))


    return holder


# def assign2(count, templates, idx):

#     holder = []

#     if idx == len(templates[0]) - 1:

#         for t, c in zip(templates, count):
#             base = list(t)
#             base[idx] = c

#             holder.append(base)

#         return holder

#     counts = []

#     for t, c in zip(templates, count):
#         for i in range(c + 1):

#             base = list(t)
#             base[idx] = i

#             holder.append(base)
#             counts.append(c - i)

#     return assign2(counts, holder, idx + 1)

def total_assignment(deme, genes):

    return deme**genes


if __name__ == '__main__':

    base = ([0, 0], 1)

    print assign(6, base, 0)

    print total_assignment(2, 3)

    print total_assignment(3, 5)

    # val1 = assign2([5], [base], 0)

    # print val0 == val1
