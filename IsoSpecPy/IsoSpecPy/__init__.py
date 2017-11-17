# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015-2016 Mateusz Łącki and Michał Startek.
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

from isoFFI import isoFFI
import re
import types
import PeriodicTbl



regex_symbols = re.compile("\D+")
regex_atom_cnts = re.compile("\d+")

def IsoParamsFromFormula(formula):
    global regex_symbols, regex_atom_cnts
    symbols = regex_symbols.findall(formula)
    atomCounts = [int(x) for x in regex_atom_cnts.findall(formula)]

    if not len(symbols) == len(atomCounts):
        raise ValueError("Invalid formula")

    try:
        masses = tuple(PeriodicTbl.symbol_to_masses[s] for s in symbols)
        probs  = tuple(PeriodicTbl.symbol_to_probs[s]  for s in symbols)
        isotopeNumbers = tuple(len(PeriodicTbl.symbol_to_probs[s]) for s in symbols)
    except KeyError:
        raise ValueError("Invalid formula")

    return (len(atomCounts), isotopeNumbers, atomCounts, masses, probs)



class Iso(object):
    def __init__(self, formula=None,
                 get_confs=False,
                 dimNumber=None,
                 isotopeNumbers=None,
                 atomCounts=None,
                 isotopeMasses=None,
                 isotopeProbabilities=None):
        """Initialize Iso. TODO write it."""

        if formula is None and not all([dimNumber, isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities]):
            raise Exception("Either formula or ALL of: dimNumber, isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities must not be None")

        if formula is not None:
            self.dimNumber, self.isotopeNumbers, self.atomCounts, \
            self.isotopeMasses, self.isotopeProbabilities = IsoParamsFromFormula(formula)

        if dimNumber is not None:
            self.dimNumber = dimNumber

        if atomCounts is not None:
            self.atomCounts = atomCounts

        if isotopeNumbers is not None:
            self.isotopeNumbers = isotopeNumbers

        if isotopeMasses is not None:
            self.isotopeMasses = isotopeMasses

        if isotopeProbabilities is not None:
            self.isotopeProbabilities = isotopeProbabilities

        self.get_confs = get_confs
        self.ffi = isoFFI.clib

        self.iso = self.ffi.setupIso(self.dimNumber, self.isotopeNumbers,
                                     self.atomCounts,
                                     [i for s in self.isotopeMasses for i in s],
                                     [i for s in self.isotopeProbabilities for i in s])

    def __iter__(self):
        if self.get_confs:
            raise NotImplementedError() # TODO: implement
        else:
            for i in xrange(self.size):
                yield (self.masses[i], self.probs[i])
                
    def __del__(self):
        self.ffi.deleteIso(self.iso)
        


class IsoThreshold(Iso):
    def __init__(self, threshold=.0001, absolute=False, get_confs = False, **kwargs):
        super(IsoThreshold, self).__init__(get_confs = get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute

        self.generator = self.ffi.setupIsoThresholdGenerator(self.iso, threshold, absolute, 1000, 1000)
        self.tabulator = self.ffi.setupThresholdTabulator(self.generator, True, True, True, get_confs)

        self.size = self.ffi.confs_noThresholdTabulator(self.tabulator)

        def c(typename, what, mult = 1):
            return isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what)

        self.masses = c("double", self.ffi.massesThresholdTabulator(self.tabulator))
        self.lprobs = c("double", self.ffi.lprobsThresholdTabulator(self.tabulator))
        self.probs  = c("double", self.ffi.probsThresholdTabulator(self.tabulator))

        if get_confs:
            self.confs = c("int", self.ffi.confsThresholdTabulator(self.tabulator), mult = self.dimNumber)

    def __len__(self):
        return self.size

    def __del__(self):
        pass



'''

class IsoLayered(Iso):
    def __init__(self, delta=-10.0, **kwargs):
        super(IsoLayered, self).__init__(**kwargs)
        self.delta = delta

    def __iter__(self):
        return IsoLayeredGenerator(src=self, get_confs=self.get_confs)

    def __del__(self):
        pass

'''

class IsoGenerator(Iso):
    def __init__(self, get_confs=False, **kwargs):
        super(IsoGenerator, self).__init__(get_confs = get_confs, **kwargs)

    def __iter__(self):
        cgen = self.cgen
        if self.get_confs:
            raise NotImplemented()
        else:
            while self.advancer(cgen):
                yield (self.mass_getter(cgen), self.lprob_getter(cgen))

        

class IsoThresholdGenerator(IsoGenerator):
    def __init__(self, threshold=0.0001, absolute=False, get_confs=False, **kwargs):
        super(IsoThresholdGenerator, self).__init__(get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute
        self.cgen = self.ffi.setupIsoThresholdGenerator(self.iso,
                                                        threshold,
                                                        absolute,
                                                        1000,
                                                        1000)
        self.advancer = self.ffi.advanceToNextConfigurationIsoThresholdGenerator
        self.lprob_getter = self.ffi.lprobIsoThresholdGenerator
        self.mass_getter = self.ffi.massIsoThresholdGenerator

    def __del__(self):
        self.ffi.deleteIsoThresholdGenerator(self.cgen)


class IsoLayeredGenerator(IsoGenerator):
    def __init__(self, delta = -10.0, get_confs=False, **kwargs):
        assert delta < 0.0
        super(IsoLayeredGenerator, self).__init__(get_confs, **kwargs)
        self.delta = delta
        self.cgen = self.ffi.setupIsoLayeredGenerator(self.iso,
                                                      self.delta,
                                                      1000,
                                                      1000)
        self.advancer = self.ffi.advanceToNextConfigurationIsoLayeredGenerator
        self.lprob_getter = self.ffi.lprobIsoLayeredGenerator
        self.mass_getter = self.ffi.massIsoLayeredGenerator

    def __del__(self):
        self.ffi.deleteIsoLayeredGenerator(self.cgen)

class IsoOrderedGenerator(IsoGenerator):
    def __init__(self, get_confs=False, **kwargs):
        super(IsoOrderedGenerator, self).__init__(get_confs, **kwargs)
        self.cgen = self.ffi.setupIsoOrderedGenerator(self.iso,
                                                      1000,
                                                      1000)
        self.advancer = self.ffi.advanceToNextConfigurationIsoOrderedGenerator
        self.lprob_getter = self.ffi.lprobIsoOrderedGenerator
        self.mass_getter = self.ffi.massIsoOrderedGenerator

    def __del__(self):
        self.ffi.deleteIsoLayeredGenerator(self.cgen)


