# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
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
from collections import namedtuple, OrderedDict
import math

try:
    xrange
except NameError:
    xrange = range

regex_pattern = re.compile('([A-Z][a-z]?)(-?[0-9]*)')
ParsedFormula = namedtuple('ParsedFormula', 'atomCounts masses probs elems')




def ParseFormula(formula):
    """Parse a chemical formula.

    Args:
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O"
    
    Returns:
        A tuple containing element symbols and atomCounts of elements in the parsed formula.
    """
    global regex_pattern

    ret = OrderedDict()

    last = 0
    for match in re.finditer(regex_pattern, formula):
        elem, cnt = match.groups()
        if elem in ret:
            raise ValueError("""Invalid formula: {} (repeating element: "{}")""".format(formula, elem))
        ret[elem] = int(cnt) if cnt != '' else 1
        if last!=match.start():
            raise ValueError("""Invalid formula: {}  (garbage inside: "{}")""".format(formula, formula[last:match.start()]))
        if elem not in PeriodicTbl.symbol_to_masses:
            raise ValueError("""Invalid formula: {} (unknown element symbol: "{}")""".format(formula, elem))
        last = match.end()

    if len(formula) != last:
        raise ValueError('''Invalid formula: {}  (trailing garbage: "{}")'''.format(formula, formula[last:]))

    if len(ret) == 0:
        raise ValueError("Invalid formula (empty)")

    return ret

fasta_parsing_space = isoFFI.ffi.new("int[6]")

def ParseFASTA(fasta):
    if isinstance(fasta, str):
        fasta = fasta.encode("ascii")
    isoFFI.clib.parse_fasta_c(fasta, fasta_parsing_space)
    elements = list("CHNOS")
    if fasta_parsing_space[5] > 0:
        elements.append("Se")
    od = OrderedDict()
    for i in range(len(elements)):
        od[elements[i]] = fasta_parsing_space[i]
    return od


def IsoParamsFromDict(formula, use_nominal_masses = False):
    """Produces a set of IsoSpec parameters from a chemical formula.

    Args:
        formula (dict): a parsed chemical formula, e.g. {"C": 2, "H": 6, "O": 1}
        use_nominal_masses (boolean): use masses of elements rounded to integer numbers (nominal masses)

    Returns:
        ParsedFormula: a tuple containing atomCounts, masses and marginal probabilities of elements in the parsed formula.
    """

    symbols, atomCounts = [], []
    for symbol, atomCount in formula.items():
        symbols.append(symbol)
        atomCounts.append(atomCount)

    try:
        if use_nominal_masses:
            masses = [PeriodicTbl.symbol_to_massNo[s] for s in symbols]
        else:
            masses = [PeriodicTbl.symbol_to_masses[s] for s in symbols]
        probs  = [PeriodicTbl.symbol_to_probs[s]  for s in symbols]
    except KeyError:
        raise ValueError("Invalid formula")

    return ParsedFormula(atomCounts, masses, probs, symbols)



class Iso(object):
    """Virtual class representing an isotopic distribution."""
    def __init__(self, formula="",
                 get_confs=False,
                 atomCounts=None,
                 isotopeMasses=None,
                 isotopeProbabilities=None,
                 use_nominal_masses = False,
                 fasta = "",
                 charge = 1.0):
        """Initialize Iso.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            get_confs (boolean): should we report counts of isotopologues?
            atomCounts (list): a list of atom counts (alternative to 'formula').
            isotopeMasses (list): a list of lists of masses of elements with counts in 'atomCounts'.
            isotopeProbabilities (list): a list of lists of probabilities of elements with counts in 'atomCounts'.
            use_nominal_masses (boolean): should the masses be rounded to the closest integer values.
            charge (float): charge state of the molecule: all masses will be divided by this value to obtain the m/z values.
        """

        self.iso = None

        if len(fasta) > 0:
            molecule = ParseFASTA(fasta)
        else:
            molecule = OrderedDict()

        if len(formula) > 0:
            if isinstance(formula, dict):
                df = formula
            else:
                df = ParseFormula(formula)
            for symbol, count in df.items():
                molecule[symbol] = molecule.get(symbol, 0) + count

        for sym, cnt in molecule.items():
            if cnt < 0:
                raise Exception("Negative count of element " + sym + ": " + str(cnt))

        if len(molecule) == 0 and not all([atomCounts, isotopeMasses, isotopeProbabilities]):
            raise Exception("Either formula, fasta or ALL of: atomCounts, isotopeMasses, isotopeProbabilities must not be None")

        if len(molecule) > 0:
            self.atomCounts, self.isotopeMasses, self.isotopeProbabilities, _ = IsoParamsFromDict(molecule, use_nominal_masses = use_nominal_masses)
        else:
            self.atomCounts, self.isotopeMasses, self.isotopeProbabilities = [], [], []

        if not (atomCounts is None):
            self.atomCounts.extend(atomCounts)

        if not (isotopeMasses is None):
            self.isotopeMasses.extend(isotopeMasses)

        if not (isotopeProbabilities is None):
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
                                     [i/charge for s in self.isotopeMasses for i in s],
                                     [i for s in self.isotopeProbabilities for i in s])

    def __del__(self):
        try:
            if self.iso is not None:
                self.ffi.deleteIso(self.iso)
                self.iso = None
        except AttributeError:
            pass

    def getLightestPeakMass(self):
        """Get the lightest peak in the isotopic distribution."""
        return self.ffi.getLightestPeakMassIso(self.iso)

    def getHeaviestPeakMass(self):
        """Get the heaviest peak in the isotopic distribution."""
        return self.ffi.getHeaviestPeakMassIso(self.iso)

    def getMonoisotopicPeakMass(self):
        """Get the monoisotopic mass of the peak."""
        return self.ffi.getMonoisotopicPeakMassIso(self.iso)

    def getModeLProb(self):
        """Get the log probability of the most probable peak(s) in the isotopic distribution."""
        return self.ffi.getModeLProbIso(self.iso)

    def getModeMass(self):
        """Get the mass of the most probable peak.

        If there are more, return only the mass of one of them."""
        return self.ffi.getModeMassIso(self.iso)

    def getTheoreticalAverageMass(self):
        return self.ffi.getTheoreticalAverageMassIso(self.iso)

    def variance(self):
        return self.ffi.getIsoVariance(self.iso)

    def stddev(self):
        return self.ffi.getIsoStddev(self.iso)

    def getMarginalLogSizeEstimates(self, prob):
        cbuf = isoFFI.clib.getMarginalLogSizeEstimates(self.iso, prob)
        ret = list(isoFFI.ffi.cast('double[' + str(self.dimNumber) + ']', cbuf))
        isoFFI.clib.freeReleasedArray(cbuf)
        return ret

    def parse_conf(self, cptr, starting_with = 0):
        return tuple(tuple(cptr[i+starting_with] for i in o) for o in self.offsets)

    def _get_parse_conf_fun(self):
        # Can't just use the above function as lambda, as the entire class instance will be in closure
        offsets = self.offsets
        def pc(cptr, starting_with = 0):
            return tuple(tuple(cptr[i+starting_with] for i in o) for o in offsets)
        return pc


class IsoDistribution(object):
    """Isotopoic distribution with precomputed vector of masses and probabilities."""
    def np_masses(self):
        """Return computed masses as a numpy array."""
        try:
            import numpy as np
        except ImportError:
            raise Exception(e.msg + "\nThis requires numpy to be installed.")
        return np.frombuffer(isoFFI.ffi.buffer(self.masses))

    def np_probs(self):
        """Return computed probabilities as a numpy array."""
        try:
            import numpy as np
        except ImportError as e:
            raise Exception(e.msg + "\nThis requires numpy to be installed.")
        return np.frombuffer(isoFFI.ffi.buffer(self.probs))

    def __iter__(self):
        if hasattr(self, "confs") and self.confs is not None:
            for i in xrange(self.size):
                yield(self.masses[i], self.probs[i], self.confs[i])
        else:
            for i in xrange(self.size):
                yield (self.masses[i], self.probs[i])

    def __getitem__(self, idx):
        try:
            return (self.masses[idx], self.probs[idx], self.confs[idx])
        except (AttributeError, TypeError):
            return (self.masses[idx], self.probs[idx])

    def __len__(self):
        """Get the number of calculated peaks."""
        return self.size

    def _get_conf(self, idx):
        return self.parse_conf(self.raw_confs, starting_with = self.sum_isotope_numbers * idx)

    def __del__(self):
        pass

    def __init__(self, cobject = None, probs = None, masses = None, get_confs = False, iso = None):
        self.mass_sorted = False
        self.prob_sorted = False
        self._total_prob = float('nan')

        if cobject is not None:
            self.size = isoFFI.clib.confs_noFixedEnvelope(cobject)

            def wrap(typename, what, attrname, mult = 1):
                if what is not None:
                    x = isoFFI.ffi.gc(isoFFI.ffi.cast(typename + '[' + str(self.size*mult) + ']', what), isoFFI.clib.freeReleasedArray)
                    setattr(self, attrname, x)

            wrap("double", isoFFI.clib.massesFixedEnvelope(cobject), "masses")
            wrap("double", isoFFI.clib.probsFixedEnvelope(cobject), "probs")

            if get_confs:
                # Must also be a subclass of Iso...
                self.sum_isotope_numbers = sum(iso.isotopeNumbers)
                wrap("int", isoFFI.clib.confsFixedEnvelope(cobject), "raw_confs", mult = self.sum_isotope_numbers)
                self.confs = ConfsPassthrough(lambda idx: self._get_conf(idx), self.size)
                self.parse_conf = iso._get_parse_conf_fun()

        elif probs is not None or masses is not None:
            assert probs is not None and masses is not None
            assert len(probs) == len(masses)
            self.size = len(probs)
            type_str = "double["+str(self.size)+"]"
            self.probs  = isoFFI.ffi.new(type_str, probs)
            self.masses = isoFFI.ffi.new(type_str, masses)
        elif cobject == probs == masses == get_confs == iso == None:
            self.size = 0
            type_str = "double["+str(self.size)+"]"
            self.probs  = isoFFI.ffi.new(type_str, [])
            self.masses = isoFFI.ffi.new(type_str, [])
            self._total_prob = 0.0
            self.mass_sorted = True
            self.prob_sorted = True
        else:
            raise RuntimeError("Invalid arguments for IsoDistribution constructor")

    def _get_cobject(self):
        return isoFFI.clib.setupFixedEnvelope(self.masses, self.probs, len(self.masses), self.mass_sorted, self.prob_sorted, self._total_prob)

    def copy(self):
        x = self._get_cobject()
        c = isoFFI.clib.copyFixedEnvelope(x)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        ret = IsoDistribution(cobject = c)
        ret._total_prob = self._total_prob
        ret.mass_sorted = self.mass_sorted
        ret.prob_sorted = self.prob_sorted
        isoFFI.clib.deleteFixedEnvelope(c, False)
        return ret

    def __add__(self, other):
        x = self._get_cobject()
        y = other._get_cobject()
        cobject = isoFFI.clib.addEnvelopes(x, y)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        ret = IsoDistribution(cobject = cobject)
        isoFFI.clib.deleteFixedEnvelope(cobject, False)
        return ret

    def __mul__(self, other):
        x = self._get_cobject()
        y = other._get_cobject()
        cobject = isoFFI.clib.convolveEnvelopes(x, y)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        ret = IsoDistribution(cobject = cobject)
        isoFFI.clib.deleteFixedEnvelope(cobject, False)
        return ret

    def total_prob(self):
        if math.isnan(self._total_prob):
            co = self._get_cobject()
            self._total_prob = isoFFI.clib.getTotalProbOfEnvelope(co)
            isoFFI.clib.deleteFixedEnvelope(co, True)
        return self._total_prob

    def normalize(self):
        co = self._get_cobject()
        isoFFI.clib.normalizeEnvelope(co)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        self._total_prob = 1.0

    def normalized(self):
        ret = self.copy()
        ret.normalize()
        return ret

    def add_mass(self, d_mass):
        isoFFI.clib.array_add(self.masses, self.size, d_mass)

    def mul_mass(self, d_mass):
        isoFFI.clib.array_mul(self.masses, self.size, d_mass)

    def add_mul_mass(self, add, mul):
        isoFFI.clib.array_fma(self.masses, self.size, mul, add*mul)

    def mul_add_mass(self, mul, add):
        isoFFI.clib.array_fma(self.masses, self.size, mul, add)

    def scale(self, factor):
        '''Multiplies the pribabilities of spectrum by factor. Works in place.'''
        co = self._get_cobject()
        isoFFI.clib.scaleEnvelope(co, factor)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        self._total_prob *= factor

    def scaled(self, factor):
        '''Returns a copy of the spectrum where each probability was multiplied by factor.'''
        ret = self.copy()
        ret.scale(factor)
        return ret

    def resample(self, ionic_current, beta_bias=1.0):
        co = self._get_cobject()
        isoFFI.clib.resampleEnvelope(co, ionic_current, beta_bias)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        self._total_prob = float(ionic_current)
        self.prob_sorted = False

    def sort_by_prob(self):
        if not self.prob_sorted:
            co = self._get_cobject()
            isoFFI.clib.sortEnvelopeByProb(co)
            isoFFI.clib.deleteFixedEnvelope(co, True)
            self.mass_sorted = False
            self.prob_sorted = True

    def sort_by_mass(self):
        if not self.mass_sorted:
            co = self._get_cobject()
            isoFFI.clib.sortEnvelopeByMass(co)
            isoFFI.clib.deleteFixedEnvelope(co, True)
            self.mass_sorted = True
            self.prob_sorted = False

    def _recalculate_everything(self):
        self._total_prob = float('nan')
        self.mass_sorted = False
        self.prob_sorted = False

    def empiric_average_mass(self):
        co = self._get_cobject()
        ret = isoFFI.clib.empiricAverageMass(co)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        return ret

    def empiric_variance(self):
        co = self._get_cobject()
        ret = isoFFI.clib.empiricVariance(co)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        return ret

    def empiric_stddev(self):
        co = self._get_cobject()
        ret = isoFFI.clib.empiricStddev(co)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        return ret

    def wassersteinDistance(self, other):
        x = self._get_cobject()
        y = other._get_cobject()
        ret = isoFFI.clib.wassersteinDistance(x, y)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        if math.isnan(ret):
            raise ValueError("Both spectra must be normalized before Wasserstein distance can be computed.")
        return ret

    def orientedWassersteinDistance(self, other):
        x = self._get_cobject()
        y = other._get_cobject()
        ret = isoFFI.clib.orientedWassersteinDistance(x, y)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        if math.isnan(ret):
            raise ValueError("Both spectra must be normalized before Wasserstein distance can be computed.")
        return ret

    def abyssalWassersteinDistance(self, other, abyss_depth, other_scale = 1.0):
        x = self._get_cobject()
        y = other._get_cobject()
        ret = isoFFI.clib.abyssalWassersteinDistance(x, y, abyss_depth, other_scale)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        return ret

    def abyssalWassersteinDistanceGrad(self, others, scales, grad, abyss_depth_exp, abyss_depth_the):
        assert len(others) + 1 == len(scales) == len(grad)
        cobjs = [self._get_cobject()]
        cobjs.extend(other._get_cobject() for other in others)
        ret = isoFFI.clib.abyssalWassersteinDistanceGrad(cobjs, scales, grad, len(others), abyss_depth_exp, abyss_depth_the)
        for cobj in cobjs:
            isoFFI.clib.deleteFixedEnvelope(cobj, True)
        self.mass_sorted = True
        self.prob_sorted = False
        for other in others:
            other.mass_sorted = True
            other.prob_sorted = False
        return ret

    def wassersteinMatch(self, other, flow_dist, other_scale = 1.0):
        x = self._get_cobject()
        y = other._get_cobject()
        ret = isoFFI.clib.wassersteinMatch(x, y, flow_dist, other_scale)
        isoFFI.clib.deleteFixedEnvelope(x, True)
        isoFFI.clib.deleteFixedEnvelope(y, True)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        return (ret.res1, ret.res2, ret.flow)

    def binned(self, width = 1.0, middle = 0.0):
        co = self._get_cobject()
        cbo = isoFFI.clib.binnedEnvelope(co, width, middle)
        isoFFI.clib.deleteFixedEnvelope(co, True)
        ret = IsoDistribution(cobject = cbo)
        isoFFI.clib.deleteFixedEnvelope(cbo, False)
        self.mass_sorted = True
        self.prob_sorted = False
        ret.mass_sorted = True
        ret.prob_sorted = False
        return ret


    @staticmethod
    def LinearCombination(envelopes, intensities):
        envelope_objs = [x._get_cobject() for x in envelopes]
        if type(intensities) != list:
            intensities = list(intensities)
        cobject = isoFFI.clib.linearCombination(envelope_objs, intensities, len(envelope_objs))
        for x in envelope_objs:
            isoFFI.clib.deleteFixedEnvelope(x, True)
        ret = IsoDistribution(cobject = cobject)
        isoFFI.clib.deleteFixedEnvelope(cobject, False)
        return ret


    def plot(self, plot_type = "bars", show = True, **matplotlib_args):
        """Plot the isotopic distribution.

        Args:
            **matplotlib_args: arguments for matplotlib plot.
        """
        try:
            from matplotlib import pyplot as plt
        except ImportError as e:
            raise ImportError(str(e) + "\nPlotting spectra requires matplotlib to be installed.")
        if plot_type == "bars":
            if "linewidth" not in matplotlib_args:
                matplotlib_args['linewidth'] = 1.0
            plt.vlines(list(self.masses), [0], list(self.probs), **matplotlib_args)
        elif plot_type == "profile":
            self.sort_by_mass()
            plt.plot(list(self.masses), list(self.probs), **matplotlib_args)
        plt.xlabel("Mass (Da)")
        plt.ylabel("Intensity (relative)")
        if show:
            plt.show()




def IsoThreshold(threshold,
                 formula="",
                 absolute=False,
                 get_confs=False,
                 **kwargs):
    """Initialize the IsoDistribution isotopic distribution by threshold.

    Args:
        threshold (float): value of the absolute or relative threshold.
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
        absolute (boolean): should we report peaks with probabilities above an absolute probability threshold, or above a relative threshold amounting to a given proportion of the most probable peak?
        get_confs (boolean): should we report counts of isotopologues?
        **kwds: named arguments to IsoSpectrum.
    """
    iso = Iso(formula = formula, get_confs = get_confs, **kwargs)
    tabulator = isoFFI.clib.setupThresholdFixedEnvelope(iso.iso, threshold, absolute, get_confs)
    ido = IsoDistribution(cobject = tabulator, get_confs = get_confs, iso = iso)
    isoFFI.clib.deleteFixedEnvelope(tabulator, False)
    return ido


def IsoTotalProb(prob_to_cover,
                 formula="",
                 get_minimal_pset=True,
                 get_confs=False,
                 **kwargs):
        """Initialize the IsoDistribution isotopic distribution by total probability.

        Args:
            prob_to_cover (float): minimal total probability of the reported peaks.
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            get_minimal_pset (boolean): should we trim the last calculated layer of isotopologues so that the reported result is as small as possible?
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        iso = Iso(formula=formula, get_confs=get_confs, **kwargs)
        tabulator = isoFFI.clib.setupTotalProbFixedEnvelope(iso.iso, prob_to_cover, get_minimal_pset, get_confs)
        ido = IsoDistribution(cobject = tabulator, get_confs = get_confs, iso = iso)
        isoFFI.clib.deleteFixedEnvelope(tabulator, False)
        return ido


def IsoStochastic(no_molecules,
                 formula="",
                 precision=0.9999,
                 beta_bias=5.0,
                 get_confs=False,
                 **kwargs):
        """Initialize the IsoDistribution isotopic distribution by total probability.

        Args:
            no_molecules (uint): ionic current in instrument
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            precision (float): passed to IsoTotalProbGenerator. Between 0.0 and 1.0.
            beta_bias (float, nonnegative): fiddling with this parameter does not change the result, but might make computations slightly faster (or likely, much, much slower is you screw it up...)
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        iso = Iso(formula=formula, get_confs=get_confs, **kwargs)
        tabulator = isoFFI.clib.setupStochasticFixedEnvelope(iso.iso, no_molecules, precision, beta_bias, get_confs)
        ido = IsoDistribution(cobject = tabulator, get_confs = get_confs, iso = iso)
        isoFFI.clib.deleteFixedEnvelope(tabulator, False)
        return ido


def IsoBinned(bin_width,
                 formula="",
                 target_total_prob=0.9999,
                 bin_middle=0.0,
                 **kwargs):
        """Initialize the IsoDistribution isotopic distribution by total probability.

        Args:
            no_molecules (uint): ionic current in instrument
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            precision (float): passed to IsoTotalProbGenerator. Between 0.0 and 1.0.
            beta_bias (float, nonnegative): fiddling with this parameter does not change the result, but might make computations slightly faster (or likely, much, much slower is you screw it up...)
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        iso = Iso(formula=formula, get_confs=False, **kwargs)
        tabulator = isoFFI.clib.setupBinnedFixedEnvelope(iso.iso, target_total_prob, bin_width, bin_middle)
        ido = IsoDistribution(cobject = tabulator, get_confs = False, iso = iso)
        isoFFI.clib.deleteFixedEnvelope(tabulator, False)
        return ido


class IsoGenerator(Iso):
    """Virtual class alowing memory-efficient iteration over the isotopic distribution.

    This iterator will stop only after enumerating all isotopologues.
    """
    def __init__(self, formula="", get_confs=False, **kwargs):
        """Initialize the IsoGenerator.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        self.cgen = None
        super(IsoGenerator, self).__init__(formula=formula, get_confs=get_confs, **kwargs)
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
    """Class alowing memory-efficient iteration over the isotopic distribution up till some probability threshold.

    This iterator will stop only after reaching a probability threshold.
    """
    def __init__(self, threshold, formula="", absolute=False, get_confs=False, reorder_marginals = True, **kwargs):
        """Initialize IsoThresholdGenerator.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            absolute (boolean): should we report peaks with probabilities above an absolute probability threshold, or above a relative threshold amounting to a given proportion of the most probable peak?
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        super(IsoThresholdGenerator, self).__init__(formula=formula, get_confs=get_confs, **kwargs)
        self.threshold = threshold
        self.absolute = absolute
        self.cgen = self.ffi.setupIsoThresholdGenerator(self.iso,
                                                        threshold,
                                                        absolute,
                                                        1000,
                                                        1000,
                                                        reorder_marginals)
        self.advancer = self.ffi.advanceToNextConfigurationIsoThresholdGenerator
        self.xprob_getter = self.ffi.probIsoThresholdGenerator
        self.mass_getter = self.ffi.massIsoThresholdGenerator
        self.conf_getter = self.ffi.get_conf_signatureIsoThresholdGenerator

    def __del__(self):
        """Destructor."""
        try:
            if self.cgen is not None:
                self.ffi.deleteIsoThresholdGenerator(self.cgen)
                self.cgen = None
        except AttributeError:
            pass
        super(IsoThresholdGenerator, self).__del__()


class IsoLayeredGenerator(IsoGenerator):
    """Class alowing memory-efficient iteration over the isotopic distribution up till some joint probability of the reported peaks."""
    def __init__(self, formula="", get_confs=False, reorder_marginals = True, t_prob_hint = 0.99, **kwargs):
        """Initialize IsoThresholdGenerator.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            absolute (boolean): should we report peaks with probabilities above an absolute probability threshold, or above a relative threshold amounting to a given proportion of the most probable peak?
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        super(IsoLayeredGenerator, self).__init__(formula=formula, get_confs=get_confs, **kwargs)
        self.cgen = self.ffi.setupIsoLayeredGenerator(self.iso, 1000, 1000, reorder_marginals, t_prob_hint)
        self.advancer = self.ffi.advanceToNextConfigurationIsoLayeredGenerator
        self.xprob_getter = self.ffi.probIsoLayeredGenerator
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
    """Class representing an isotopic distribution with peaks ordered by descending probability.

    This generator return probilities ordered with descending probability.
    It it not optimal to do so, but it might be useful.

    WARNING! This algorithm work in O(N*log(N)) vs O(N) of the threshold and layered algorithms.
    Also, the order of descending probability will most likely not reflect the order of ascending masses.
    """
    def __init__(self, formula="", get_confs=False, **kwargs):
        """Initialize IsoOrderedGenerator.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        super(IsoOrderedGenerator, self).__init__(formula=formula, get_confs=get_confs, **kwargs)
        self.cgen = self.ffi.setupIsoOrderedGenerator(self.iso, 1000, 1000)
        self.advancer = self.ffi.advanceToNextConfigurationIsoOrderedGenerator
        self.xprob_getter = self.ffi.probIsoOrderedGenerator
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


class IsoStochasticGenerator(IsoGenerator):
    def __init__(self, no_molecules, formula="", precision=0.9999, beta_bias = 1.0, get_confs=False, **kwargs):
        super(IsoStochasticGenerator, self).__init__(formula=formula, get_confs=get_confs, **kwargs)
        self.threshold = precision
        self.no_molecules = no_molecules
        self.cgen = self.ffi.setupIsoStochasticGenerator(self.iso,
                                                        no_molecules,
                                                        precision,
                                                        beta_bias)
        self.advancer = self.ffi.advanceToNextConfigurationIsoStochasticGenerator
        self.xprob_getter = self.ffi.probIsoStochasticGenerator
        self.mass_getter = self.ffi.massIsoStochasticGenerator
        self.conf_getter = self.ffi.get_conf_signatureIsoStochasticGenerator

    def __del__(self):
        """Destructor."""
        try:
            if self.cgen is not None:
                self.ffi.deleteIsoStochasticGenerator(self.cgen)
                self.cgen = None
        except AttributeError:
            pass
        super(IsoStochasticGenerator, self).__del__()

