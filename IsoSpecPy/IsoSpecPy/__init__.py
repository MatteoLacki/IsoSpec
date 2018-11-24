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

from .isoFFI import isoFFI
import re
import types
from . import PeriodicTbl
from .confs_passthrough import ConfsPassthrough


try:
    xrange
except NameError:
    xrange = range

regex_pattern = re.compile('([A-Z][a-z]?)([0-9]*)')

def IsoParamsFromFormula(formula):
    """Parse a chemical formula.

    Arguments
    ---------
    formula : str
        A string with a chemical formula.

    Examples
    --------
    >>> IsoParamsFromFormula("C100H200")
    ([100, 200],
     ((12.0, 13.0033548352), (1.00782503227, 2.01410177819)),
     ((0.9892119418504669, 0.010788058149533084),
      (0.9998842901643079, 0.00011570983569203331)))
    """
    global regex_pattern

    symbols = []
    atomCounts = []
    for elem, cnt in re.findall(regex_pattern, formula):
        symbols.append(elem)
        atomCounts.append(int(cnt) if cnt is not '' else 1)
    try:
        masses = tuple(PeriodicTbl.symbol_to_masses[s] for s in symbols)
        probs  = tuple(PeriodicTbl.symbol_to_probs[s]  for s in symbols)
    except KeyError:
        raise ValueError("Invalid formula")

    return (atomCounts, masses, probs)



class Iso(object):
    def __init__(self, formula=None,
                 get_confs=False,
                 atomCounts=None,
                 isotopeMasses=None,
                 isotopeProbabilities=None):
        """Initialize Iso."""

        self.iso = None

        if formula is None and not all([atomCounts, isotopeMasses, isotopeProbabilities]):
            raise Exception("Either formula or ALL of: atomCounts, isotopeMasses, isotopeProbabilities must not be None")

        if formula is not None:
            self.atomCounts, self.isotopeMasses, self.isotopeProbabilities = IsoParamsFromFormula(formula)

        if atomCounts is not None:
            self.atomCounts = atomCounts


        if isotopeMasses is not None:
            self.isotopeMasses = isotopeMasses

        if isotopeProbabilities is not None:
            self.isotopeProbabilities = isotopeProbabilities

        for sublist in self.isotopeProbabilities:
            for prob in sublist:
                if not (0.0 < prob <= 1.0):
                    raise ValueError("All isotope probabilities p must fulfill: 0.0 < p <= 1.0")

        self.isotopeNumbers = tuple(map(len, self.isotopeMasses))
        assert self.isotopeNumbers == tuple(map(len, self.isotopeProbabilities))

        self.dimNumber = len(self.isotopeNumbers)

        self.get_confs = get_confs
        self.ffi = isoFFI.clib

        offsets = []

        if get_confs:
            i = 0
            for j in xrange(self.dimNumber):
                newl = []
                for k in xrange(self.isotopeNumbers[j]):
                    newl.append(i)
                    i += 1
                offsets.append(tuple(newl))
            self.offsets = tuple(offsets)

        self.iso = self.ffi.setupIso(self.dimNumber, self.isotopeNumbers,
                                     self.atomCounts,
                                     [i for s in self.isotopeMasses for i in s],
                                     [i for s in self.isotopeProbabilities for i in s])

    def __iter__(self):
        if self.get_confs:
            for i in xrange(self.size):
                yield(self.masses[i], self.probs[i], self.confs[i])
        else:
            for i in xrange(self.size):
                yield (self.masses[i], self.probs[i])

    def __del__(self):
        if self.iso is not None:
            self.ffi.deleteIso(self.iso)


    def parse_conf(self, cptr, starting_with = 0):
        return tuple(tuple(cptr[i+starting_with] for i in o) for o in self.offsets)



class IsoThreshold(Iso):
    """The thresholded version of the algorithm.

    Arguments
    ---------
    threshold : float
        If absolute is False, the percentage of the most abundant peak above which 
        we report all the isotopologues. If absolute is True, the probability that
        bounds the reported isotopologues from below.
    """
    def __init__(self, threshold, absolute=False, get_confs = False, **kwargs):
        super(IsoThreshold, self).__init__(get_confs = get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute

        generator = self.ffi.setupIsoThresholdGenerator(self.iso, threshold, absolute, 1000, 1000)
        tabulator = self.ffi.setupThresholdTabulator(generator, True, True, True, get_confs)

        self.size = self.ffi.confs_noThresholdTabulator(tabulator)

        def c(typename, what, mult = 1):
            return isoFFI.ffi.gc(isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what), self.ffi.freeReleasedArray)

        self.masses = c("double", self.ffi.massesThresholdTabulator(tabulator))
        self.lprobs = c("double", self.ffi.lprobsThresholdTabulator(tabulator))
        self.probs  = c("double", self.ffi.probsThresholdTabulator(tabulator))

        if get_confs:
            self.sum_isotope_numbers = sum(self.isotopeNumbers)
            self.raw_confs = c("int", self.ffi.confsThresholdTabulator(tabulator), mult = self.sum_isotope_numbers)
            self.confs = ConfsPassthrough(lambda idx: self._get_conf(idx), self.size)

        self.ffi.deleteThresholdTabulator(tabulator)
        self.ffi.deleteIsoThresholdGenerator(generator)


    def _get_conf(self, idx):
        return self.parse_conf(self.raw_confs, starting_with = self.sum_isotope_numbers * idx)

    def __len__(self):
        return self.size




class IsoLayered(Iso):
    def __init__(self, prob_to_cover, get_confs = False, **kwargs):
        super(IsoLayered, self).__init__(get_confs = get_confs, **kwargs)
        self.prob_to_cover = prob_to_cover

        generator = self.ffi.setupIsoLayeredGenerator(self.iso, prob_to_cover, 0.3, 1000, 1000, True)
        tabulator = self.ffi.setupLayeredTabulator(generator, True, True, True, get_confs)

        self.size = self.ffi.confs_noLayeredTabulator(tabulator)

        def c(typename, what, mult = 1):
            return isoFFI.ffi.gc(isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what), self.ffi.freeReleasedArray)

        self.masses = c("double", self.ffi.massesLayeredTabulator(tabulator))
        self.lprobs = c("double", self.ffi.lprobsLayeredTabulator(tabulator))
        self.probs  = c("double", self.ffi.probsLayeredTabulator(tabulator))

        if get_confs:
            self.sum_isotope_numbers = sum(self.isotopeNumbers)
            self.raw_confs = c("int", self.ffi.confsLayeredTabulator(tabulator), mult = self.sum_isotope_numbers)
            self.confs = ConfsPassthrough(lambda idx: self._get_conf(idx), self.size)

        self.ffi.deleteLayeredTabulator(tabulator)
        self.ffi.deleteIsoLayeredGenerator(generator)

    def _get_conf(self, idx):
        return self.parse_conf(self.raw_confs, starting_with = self.sum_isotope_numbers * idx)

    def __len__(self):
        return self.size



class IsoGenerator(Iso):
    def __init__(self, get_confs=False, **kwargs):
        self.cgen = None
        super(IsoGenerator, self).__init__(get_confs = get_confs, **kwargs)
        self.conf_space = isoFFI.ffi.new("int[" + str(sum(self.isotopeNumbers)) + "]")
        self.firstuse = True

    def __iter__(self):
        if not self.firstuse:
            raise NotImplementedError("Multiple iterations through the same IsoGenerator object are not supported. Either create a new (identical) generator for a second loop-through, or use one of the non-generator classes, which do support being re-used.")
        self.firstuse = False
        cgen = self.cgen
        if self.get_confs:
            while self.advancer(cgen):
                self.conf_getter(cgen, self.conf_space)
                yield (self.mass_getter(cgen), self.xprob_getter(cgen), self.parse_conf(self.conf_space))
        else:
            while self.advancer(cgen):
                yield (self.mass_getter(cgen), self.xprob_getter(cgen))



class IsoThresholdGenerator(IsoGenerator):
    def __init__(self, threshold, absolute=False, get_confs=False, use_lprobs=False, **kwargs):
        super(IsoThresholdGenerator, self).__init__(get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute
        self.cgen = self.ffi.setupIsoThresholdGenerator(self.iso,
                                                        threshold,
                                                        absolute,
                                                        1000,
                                                        1000)
        self.advancer     = self.ffi.advanceToNextConfigurationIsoThresholdGenerator
        self.xprob_getter = self.ffi.lprobIsoThresholdGenerator if use_lprobs else self.ffi.probIsoThresholdGenerator
        self.mass_getter  = self.ffi.massIsoThresholdGenerator
        self.conf_getter  = self.ffi.get_conf_signatureIsoThresholdGenerator

    def __del__(self):
        if self.cgen is not None:
            self.ffi.deleteIsoThresholdGenerator(self.cgen)


class IsoLayeredGenerator(IsoGenerator):
    def __init__(self, prob_to_cover = 0.99, get_confs=False, do_trim = False, use_lprobs=False, **kwargs):
        super(IsoLayeredGenerator, self).__init__(get_confs, **kwargs)
        self.delta = prob_to_cover
        self.cgen  = self.ffi.setupIsoLayeredGenerator(self.iso,
                                                      self.delta,
                                                      0.3,
                                                      1000,
                                                      1000,
                                                      do_trim)
        self.advancer     = self.ffi.advanceToNextConfigurationIsoLayeredGenerator
        self.xprob_getter = self.ffi.lprobIsoLayeredGenerator if use_lprobs else self.ffi.probIsoLayeredGenerator
        self.mass_getter  = self.ffi.massIsoLayeredGenerator
        self.conf_getter  = self.ffi.get_conf_signatureIsoLayeredGenerator

    def __del__(self):
        if self.cgen is not None:
            self.ffi.deleteIsoLayeredGenerator(self.cgen)

class IsoOrderedGenerator(IsoGenerator):
    def __init__(self, get_confs=False, use_lprobs=False, **kwargs):
        super(IsoOrderedGenerator, self).__init__(get_confs, **kwargs)
        self.cgen = self.ffi.setupIsoOrderedGenerator(self.iso,
                                                      1000,
                                                      1000)
        self.advancer = self.ffi.advanceToNextConfigurationIsoOrderedGenerator
        self.xprob_getter = self.ffi.lprobIsoOrderedGenerator if use_lprobs else self.ffi.probIsoOrderedGenerator
        self.mass_getter = self.ffi.massIsoOrderedGenerator
        self.conf_getter = self.ffi.get_conf_signatureIsoOrderedGenerator

    def __del__(self):
        if self.cgen is not None:
            self.ffi.deleteIsoLayeredGenerator(self.cgen)


def isospecify(formula,
               P,
               algorithm = "layered_gen",
               get_counts = False):
    """Run IsoSpec.

    Arguments
    ---------
    formula : str
        A chemical formula.
    P : float
        Either the joint probability of reported isotopologues or minimal probability
        in that set.
    algorithm : str
        Either:
        "layered",
            for O(n) bulk calculation of isotopologues together more probable than P§,
        "threshold",
            for O(n) bulk calculation of isotopologues each more probable than P,
        "threshold_gen", 
            for O(n) generation of isotopologues each more probable than P,
        "layered_gen", 
            for O(n) generation of isotopologues together more probable than P,
        "ordered_gen",
            for O(nlog(n)) generation of isotopologues with decreasing probability,
            but YOU HAVE TO TAKE CARE OF STOPPING HERE.
        where n is the total number of isotopologues.
    get_counts : boolean
        Should the result include counts of isotopologues? By default no.

    Returns
    -------
    An iterable of tuples (mass, probability) if get_counts=False or
    (mass, probability, counts), if get_counts=True.
    """
    iso = None
    if algorithm == "layered":
        iso = IsoLayered(formula        = formula,
                         prob_to_cover  = P,
                         get_confs      = get_counts)
    if algorithm == "threshold":
        iso = IsoThreshold(formula      = formula,
                           threshold    = P,
                           get_confs    = get_counts)
    if algorithm == "threshold_gen":
        iso = IsoThresholdGenerator(formula      = formula,
                                    threshold    = P,
                                    get_confs    = get_counts)
    if algorithm == "layered_gen":
        iso = IsoLayeredGenerator(formula        = formula,
                                  prob_to_cover  = P,
                                  get_confs      = get_counts)
    if algorithm == "ordered_gen":
        iso = IsoLayeredGenerator(formula        = formula,
                                  get_confs      = get_counts)
    if iso is None:
        raise KeyError("Algorithm {} is not implemented".format(algorithm))
    return iso


__version__ = "1.9.1"


# For compatibility with 1.0.X

class CompatIsoWrapper(object):
    def __init__(self):
        from .IsoSpecPy import IsoSpec, IsoSpecify, IsoPlot
        self.IsoSpec = IsoSpec
        self.IsoSpecify = IsoSpecify
        self.IsoPlot = IsoPlot


IsoSpecPy = CompatIsoWrapper()

