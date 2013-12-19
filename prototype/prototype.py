# -*- mode: python; coding: utf-8; -*-

# test.py - brief description

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

from operator import mul

def binom(n, k):
    return int(reduce(mul, [float(n - (k - i)) / i for i in range(1, k + 1)], 1))


def state2idx(init, state):
    ndeme = len(init)
    dim = [1] + [binom(i + ndeme - 1, i) for i in init[:-1]]
    dim = [reduce(mul, dim[:i+1], 1) for i in range(ndeme)]
    idx = [genidx(state[i * ndeme:(i + 1) * ndeme]) for i in range(ndeme)]
    return sum([i * j for i, j in zip(idx, dim)])

def idx2state(init, idx):
    ndeme = len(init)
    dim = [binom(i + ndeme - 1, i) for i in init[:-1]]
    dim2 = [reduce(mul, dim[:i+1], 1) for i in range(ndeme - 1)]
    state = []
    [state.extend(genstate(ix, ndeme, gene))
     for ix, gene
     in zip([i % j for i, j in zip([idx / k for k in [1] + dim2], dim + [1])],
            init)]
    return state


def neighbors(init, idx):
    ndeme = len(init)
    state = idx2state(init, idx)
    ids = []
    src = [i for i in range(len(state)) if state[i] > 0]
    for s in src:
        deme = s / ndeme
        tar = [i for i in range(ndeme * deme, ndeme * (deme + 1)) if i != s]
        for t in tar:
            new = list(state)
            new[s] -= 1
            new[t] += 1
            ids.append(state2idx(init, new))
    return sorted(ids)

def genidx(state):
    if len(state) == 1:
        return 0

    head, tail = state[0], state[1:]
    total = sum(state)
    return sum([binom(total - i + len(tail) - 1, total - i) for i in range(0, head)]) \
        + genidx(tail)


def genstate(idx, ndeme, ngenes):
    if ndeme == 1:
        return [ngenes]
    offsets = [0] + [binom(ngenes - i + ndeme - 2, ngenes - i) for i in range(0, ngenes)]
    offset = -1
    while len(offsets) > 0 and idx >= offsets[0]:
        idx -= offsets[0]
        offset += 1
        offsets = offsets[1:]
    return [offset] + genstate(idx, ndeme - 1, ngenes - offset)


if __name__ == '__main__':

    # print state2idx([1, 2, 3], [1, 0, 0, 0, 2, 0, 0, 0, 3])
    # print idx2state([1, 2, 3], 8)
    # print neighbors([1, 2, 3], 8)
    # print genidx([0, 0, 3])
    # print genidx([1, 0, 0])
    # print genidx([0, 1, 0])
    # print genidx([0, 0, 1])
    # print genidx([2, 0, 0])
    # print genidx([1, 1, 0])
    # print genidx([1, 0, 1])
    # print genidx([0, 2, 0])
    # print genidx([0, 1, 1])
    # print genidx([0, 0, 2])
    # print genidx([0, 0, 1])
    # print genidx([1, 0])
    # print genidx([0, 1])
    # print genidx([2, 0])
    # print genidx([1, 1])
    # print genidx([0, 2])

    # print genstate(5, 3, 2)
    # print genstate(4, 3, 2)
    # print genstate(3, 3, 2)
    # print genstate(2, 3, 2)
    # print genstate(1, 3, 2)
    # print genstate(0, 3, 2)

    # print genidx([1, 0, 0])
    # print genidx([0, 2, 0])
    # print genidx([0, 0, 3])

    # print genstate(2, 3, 1)
    # print genstate(2, 3, 2)
    # print genstate(0, 3, 3)
    pass
