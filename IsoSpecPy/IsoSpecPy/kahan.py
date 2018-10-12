# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
#
#   This file is part of IsoSpec.
#
#   IsoSpec is free software: you can redistribute it and/or modify
#   it under the terms of the Simplified ("2-clause") BSD licence.
#
#   IsoSpec is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#
#   You should have received a copy of the Simplified BSD Licence
#   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
# 



class SSummator:
    def __init__(self):
        self.partials = []  # sorted, non-overlapping partial sums

    def add(self, x):
        i = 0
        for y in self.partials:
            if abs(x) < abs(y):
                x, y = y, x
            hi = x + y
            lo = y - (hi - x)
            if lo:
                self.partials[i] = lo
                i += 1
            x = hi
        self.partials[i:] = [x]
    def get(self):
        return sum(self.partials)


class Summator:
    def __init__(self, keep_partials = False):
        self.sum = 0.0
        self.c = 0.0
        self.keep_partials = keep_partials
        self.partials = []

    def add(self, what):
        y = what - self.c
        t = self.sum + y
        c = (t - self.sum) - y
        self.sum = t
        if self.keep_partials:
            self.partials.append(self.sum)

    def get(self):
        return self.sum

    def get_partials(self):
        return self.partials

class PosNegSummator:
    def __init__(self):
        self.pos = Summator()
        self.neg = Summator()
    def add(self, what):
        if what >= 0:
            self.pos.add(what)
        else:
            self.neg.add(-what)
    def get(self):
        return self.pos.get() - self.neg.get()



