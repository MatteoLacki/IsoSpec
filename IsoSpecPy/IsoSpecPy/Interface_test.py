import IsoSpecPy
from isoFFI import isoFFI
import re
import types
import PeriodicTbl

def IsoSpec(method="layered", **args):
    pass


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

class IsoGenerator(object):
    def __init__(self, src, get_confs=False):
        self.ffi = isoFFI.clib
        self.cgen = src.get_generator()
        self.src = src
        if get_confs:
            self.next = types.MethodType(self.next_confs, self)

    def next(self):
        cgen = self.cgen
        if self.advancer(cgen):
            return (self.mass_getter(cgen), self.lprob_getter(cgen))
        raise StopIteration()

    def next_confs(self):
        raise NotImplementedException()
        

class IsoThresholdGenerator(IsoGenerator):
    def __init__(self, src, get_confs=False):
        super(IsoThresholdGenerator, self).__init__(src, get_confs)
        self.advancer = self.ffi.advanceToNextConfigurationIsoThresholdGenerator
        self.lprob_getter = self.ffi.lprobIsoThresholdGenerator
        self.mass_getter = self.ffi.massIsoThresholdGenerator

    def __del__(self):
        self.ffi.deleteIsoThresholdGenerator(self.cgen)


class IsoLayeredGenerator(IsoGenerator):
    def __init__(self, src, get_confs=False):
        super(IsoLayeredGenerator, self).__init__(src, get_confs)
        self.advancer = self.ffi.advanceToNextConfigurationIsoLayeredGenerator
        self.lprob_getter = self.ffi.lprobIsoLayeredGenerator
        self.mass_getter = self.ffi.massIsoLayeredGenerator

    def __del__(self):
        self.ffi.deleteIsoLayeredGenerator(self.cgen)

class IsoOrderedGenerator(IsoGenerator):
    def __init__(self, src, get_confs=False):
        super(IsoOrderedGenerator, self).__init__(src, get_confs)
        self.advancer = self.ffi.advanceToNextConfigurationIsoOrderedGenerator
        self.lprob_getter = self.ffi.lprobIsoOrderedGenerator
        self.mass_getter = self.ffi.massIsoOrderedGenerator

    def __del__(self):
        self.ffi.deleteIsoLayeredGenerator(self.cgen)


class Iso(object):
    def __init__(self, formula=None,
                 get_confs=False,
                 dimNumber=None,
                 isotopeNumbers=None,
                 atomCounts=None,
                 isotopeMasses=None,
                 isotopeProbabilities=None):
        """Initialize Iso. TODO write it."""

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
        raise NotImplementedError()

    def __del__(self):
        self.ffi.deleteIso(self.iso)


class IsoThreshold(Iso):
    def __init__(self, threshold=.0001, absolute=False, **kwargs):
        super(IsoThreshold, self).__init__(**kwargs)
        self.threshold = threshold
        self.absolute = absolute

    def get_generator(self):
        return self.ffi.setupIsoThresholdGenerator(self.iso,
                                                   self.threshold,
                                                   self.absolute,
                                                   1000,
                                                   1000)

    def __iter__(self):
        return IsoThresholdGenerator(src=self, get_confs=self.get_confs)

    def __del__(self):
        pass



class IsoLayered(Iso):
    def __init__(self, delta=-10.0, **kwargs):
        super(IsoLayered, self).__init__(**kwargs)
        self.delta = delta

    def get_generator(self):
        return self.ffi.setupIsoLayeredGenerator(self.iso,
                                                self.delta,
                                                1000,
                                                1000)

    def __iter__(self):
        return IsoLayeredGenerator(src=self, get_confs=self.get_confs)

    def __del__(self):
        pass


class IsoOrdered(Iso):
    def __init__(self, **kwargs):
        super(IsoOrdered, self).__init__(**kwargs)

    def get_generator(self):
        return self.ffi.setupIsoOrderedGenerator(self.iso,
                                                 1000,
                                                 1000)

    def __iter__(self):
        return IsoOrderedGenerator(src=self, get_confs=self.get_confs)

    def __del__(self):
        pass



