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
from collections import namedtuple

try:
    xrange
except NameError:
    xrange = range

regex_pattern = re.compile('([A-Z][a-z]?)([0-9]*)')
ParsedFormula = namedtuple('ParsedFormula', 'atomCount mass prob')


def IsoParamsFromFormula(formula):
    global regex_pattern

    symbols = []
    atomCounts = []
    for elem, cnt in re.findall(regex_pattern, formula):
        symbols.append(elem)
        atomCounts.append(int(cnt) if cnt is not '' else 1)
    try:
        masses = [PeriodicTbl.symbol_to_masses[s] for s in symbols]
        probs  = [PeriodicTbl.symbol_to_probs[s]  for s in symbols]
    except KeyError:
        raise ValueError("Invalid formula")

    return ParsedFormula(atomCounts, masses, probs)



class Iso(object):
    def __init__(self, formula="",
                 get_confs=False,
                 atomCounts=[],
                 isotopeMasses=[],
                 isotopeProbabilities=[]):
        """Initialize Iso."""

        self.iso = None

        if not formula and not all([atomCounts, isotopeMasses, isotopeProbabilities]):
            raise Exception("Either formula or ALL of: atomCounts, isotopeMasses, isotopeProbabilities must not be None")

        if formula:
            if isinstance(formula, dict):
                formula = ''.join(key + str(val) for key, val in formula.items() if val > 0)
            self.atomCounts, self.isotopeMasses, self.isotopeProbabilities = IsoParamsFromFormula(formula)
        else:
            self.atomCounts, self.isotopeMasses, self.isotopeProbabilities = [], [], []

        if atomCounts:
            self.atomCounts.extend(atomCounts)

        if isotopeMasses:
            self.isotopeMasses.extend(isotopeMasses)

        if isotopeProbabilities:
            self.isotopeProbabilities.extend(isotopeProbabilities)

        for sublist in self.isotopeProbabilities:
            for prob in sublist:
                if not (0.0 < prob <= 1.0):
                    raise ValueError("All isotope probabilities p must fulfill: 0.0 < p <= 1.0")

        self.isotopeNumbers = tuple(map(len, self.isotopeMasses))
        assert self.isotopeNumbers == tuple(map(len, self.isotopeProbabilities))
        assert len(self.atomCounts) == len(self.isotopeNumbers) == len(self.isotopeProbabilities)

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
        try:
            if self.iso is not None:
                self.ffi.deleteIso(self.iso)
                self.iso = None
        except AttributeError:
            pass

    def getLightestPeakMass(self):
        return self.ffi.getLightestPeakMassIso(self.iso)

    def getHeaviestPeakMass(self):
        return self.ffi.getHeaviestPeakMassIso(self.iso)

    def getMonoisotopicPeakMass(self):
        return self.ffi.getMonoisotopicPeakMassIso(self.iso)

    def getModeLProb(self):
        return self.ffi.getModeLProbIso(self.iso)

    def getModeMass(self):
        return self.ffi.getModeMassIso(self.iso)

    def getTheoreticalAverageMass(self):
        return self.ffi.getTheoreticalAverageMassIso(self.iso)

    def parse_conf(self, cptr, starting_with = 0):
        return tuple(tuple(cptr[i+starting_with] for i in o) for o in self.offsets)



class IsoThreshold(Iso):
    def __init__(self, threshold, absolute=False, get_confs = False, **kwargs):
        super(IsoThreshold, self).__init__(get_confs = get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute

        self.tabulator = self.ffi.setupThresholdFixedEnvelope(self.iso, threshold, absolute, get_confs, True, True, True)

        self.size = self.ffi.confs_noThresholdFixedEnvelope(self.tabulator)

        def c(typename, what, mult = 1):
            return isoFFI.ffi.gc(isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what), self.ffi.freeReleasedArray)

        self.masses = c("double", self.ffi.massesThresholdFixedEnvelope(self.tabulator))
        self.lprobs = c("double", self.ffi.lprobsThresholdFixedEnvelope(self.tabulator))
        self.probs  = c("double", self.ffi.probsThresholdFixedEnvelope(self.tabulator))

        if get_confs:
            self.sum_isotope_numbers = sum(self.isotopeNumbers)
            self.raw_confs = c("int", self.ffi.confsThresholdFixedEnvelope(self.tabulator), mult = self.sum_isotope_numbers)
            self.confs = ConfsPassthrough(lambda idx: self._get_conf(idx), self.size)

    def __del__(self):
        try:
            if self.tabulator != None:
                self.ffi.deleteThresholdFixedEnvelope(self.tabulator)
                self.tabulator = None
        except AttributeError:
            pass
        super(IsoThreshold, self).__del__()

    def _get_conf(self, idx):
        return self.parse_conf(self.raw_confs, starting_with = self.sum_isotope_numbers * idx)

    def __len__(self):
        return self.size




class IsoTotalProb(Iso):
    def __init__(self, prob_to_cover, get_minimal_pset = True, get_confs = False, **kwargs):
        super(IsoTotalProb, self).__init__(get_confs = get_confs, **kwargs)
        self.prob_to_cover = prob_to_cover

        self.tabulator = self.ffi.setupTotalProbFixedEnvelope(self.iso, prob_to_cover, get_minimal_pset, get_confs, True, True, True)

        self.size = self.ffi.confs_noTotalProbFixedEnvelope(self.tabulator)

        def c(typename, what, mult = 1):
            return isoFFI.ffi.gc(isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what), self.ffi.freeReleasedArray)

        self.masses = c("double", self.ffi.massesTotalProbFixedEnvelope(self.tabulator))
        self.lprobs = c("double", self.ffi.lprobsTotalProbFixedEnvelope(self.tabulator))
        self.probs  = c("double", self.ffi.probsTotalProbFixedEnvelope(self.tabulator))

        if get_confs:
            self.sum_isotope_numbers = sum(self.isotopeNumbers)
            self.raw_confs = c("int", self.ffi.confsTotalProbFixedEnvelope(self.tabulator), mult = self.sum_isotope_numbers)
            self.confs = ConfsPassthrough(lambda idx: self._get_conf(idx), self.size)

    def _get_conf(self, idx):
        return self.parse_conf(self.raw_confs, starting_with = self.sum_isotope_numbers * idx)

    def __len__(self):
        return self.size

    def __del__(self):
        try:
            if self.tabulator != None:
                self.ffi.deleteTotalProbFixedEnvelope(self.tabulator)
                self.tabulator = None
        except AttributeError:
            pass
        super(IsoTotalProb, self).__del__()


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

        
    def __del__(self):
        super(IsoGenerator, self).__del__()


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
        self.advancer = self.ffi.advanceToNextConfigurationIsoThresholdGenerator
        self.xprob_getter = self.ffi.lprobIsoThresholdGenerator if use_lprobs else self.ffi.probIsoThresholdGenerator
        self.mass_getter = self.ffi.massIsoThresholdGenerator
        self.conf_getter = self.ffi.get_conf_signatureIsoThresholdGenerator

    def __del__(self):
        try:
            if self.cgen is not None:
                self.ffi.deleteIsoThresholdGenerator(self.cgen)
                self.cgen = None
        except AttributeError:
            pass
        super(IsoThresholdGenerator, self).__del__()


class IsoLayeredGenerator(IsoGenerator):
    def __init__(self, get_confs=False, use_lprobs=False, **kwargs):
        super(IsoLayeredGenerator, self).__init__(get_confs, **kwargs)
        self.cgen = self.ffi.setupIsoLayeredGenerator(self.iso,
                                                      1000,
                                                      1000,
                                                      )
        self.advancer = self.ffi.advanceToNextConfigurationIsoLayeredGenerator
        self.xprob_getter = self.ffi.lprobIsoLayeredGenerator if use_lprobs else self.ffi.probIsoLayeredGenerator
        self.mass_getter = self.ffi.massIsoLayeredGenerator
        self.conf_getter = self.ffi.get_conf_signatureIsoLayeredGenerator

    def __del__(self):
        try:
            if self.cgen is not None:
                self.ffi.deleteIsoLayeredGenerator(self.cgen)
                self.cgen = None
        except AttributeError:
            pass
        super(IsoLayeredGenerator, self).__del__()


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
        try:
            if self.cgen is not None:
                self.ffi.deleteIsoLayeredGenerator(self.cgen)
                self.cgen = None
        except AttributeError:
            pass
        super(IsoOrderedGenerator, self).__del__()
