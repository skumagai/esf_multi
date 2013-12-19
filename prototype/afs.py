# -*- mode: python; coding: utf-8; -*-

# afs.py - Prototype for allele frequency spectrum

# Copyright (C) 2013 Seiji Kumagai

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

def combine(seqs, ndeme, vec):

    if len(seqs) == 0:
        return [vec]

    head, rest = seqs[0], seqs[1:]
    data = []
    for v in relocate_rec(head, [0] * ndeme):
        w = list(vec)
        w.extend(v)
        for x in combine(rest, ndeme, w):
            data.append(x)

    return [list(s) for s in set(tuple(i) for i in data)]


def relocate_rec(n, data):
    if n == 0:
        return [data]

    ndeme = len(data)

    ret = []
    for i in range(ndeme):
        vec = [0] * ndeme
        vec[i] = 1
        val = add(vec, data)
        for v in relocate_rec(n - 1, val):
            ret.append(v)

    return ret


def add(*lists):
    return [sum(vals) for vals in zip(*lists)]


def all_states(afs):

    ndeme = len(afs[0])
    data = [combine(a, ndeme, []) for a in afs]
    return subf(data)

def subf(data):
    if len(data) == 1:
        return [[d] for d in data[0]]

    head, rest = data[0], data[1:]

    ret = []
    for h in head:
        for i in subf(rest):
            ret.append([h] + i)

    return ret


if __name__ == '__main__':

    allele0 = [1, 3]
    allele1 = [2, 1]

    # afs = [allele0, allele1]
    allele2 = [1, 1]
    afs = [allele0, allele1, allele2]

    # print combine(allele0, len(allele0), [])
    # print add([1,2,3], [3,4,5])

    print all_states(afs)
