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

"""
IsoSpec Python bindings using nanobind - maintains the same API as the CFFI version.
"""

import re
import math
import numpy as np
from collections import namedtuple, OrderedDict
from . import PeriodicTbl
from .confs_passthrough import ConfsPassthrough

try:
    from . import _isospec_nb as nb
except ImportError:
    import _isospec_nb as nb

regex_pattern = re.compile("([A-Z][a-z]?)(-?[0-9]*)")
ParsedFormula = namedtuple("ParsedFormula", "atomCounts masses probs elems")


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
            raise ValueError(
                """Invalid formula: {} (repeating element: "{}")""".format(
                    formula, elem
                )
            )
        ret[elem] = int(cnt) if cnt != "" else 1
        if last != match.start():
            raise ValueError(
                """Invalid formula: {}  (garbage inside: "{}")""".format(
                    formula, formula[last : match.start()]
                )
            )
        if elem not in PeriodicTbl.symbol_to_masses:
            raise ValueError(
                """Invalid formula: {} (unknown element symbol: "{}")""".format(
                    formula, elem
                )
            )
        last = match.end()

    if len(formula) != last:
        raise ValueError(
            """Invalid formula: {}  (trailing garbage: "{}")""".format(
                formula, formula[last:]
            )
        )

    if len(ret) == 0:
        raise ValueError("Invalid formula (empty)")

    return ret


def ParseFASTA(fasta):
    if isinstance(fasta, str):
        fasta = fasta.encode("ascii")
    atomCounts = nb.parse_fasta(fasta.decode("ascii"))
    elements = list("CHNOS")
    if atomCounts[5] > 0:
        elements.append("Se")
    od = OrderedDict()
    for i in range(len(elements)):
        od[elements[i]] = atomCounts[i]
    return od


def IsoParamsFromDict(formula, use_nominal_masses=False):
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
        probs = [PeriodicTbl.symbol_to_probs[s] for s in symbols]
    except KeyError:
        raise ValueError("Invalid formula")

    return ParsedFormula(atomCounts, masses, probs, symbols)


def IsoParamsFromFormula(formula, use_nominal_masses=False):
    """Produces a set of IsoSpec parameters from a chemical formula.

    Args:
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O"
        use_nominal_masses (boolean): use masses of elements rounded to integer numbers (nominal masses)

    Returns:
        ParsedFormula: a tuple containing atomCounts, masses and marginal probabilities of elements in the parsed formula.
    """
    parsed = ParseFormula(formula)
    return IsoParamsFromDict(parsed, use_nominal_masses=use_nominal_masses)


class Iso(object):
    """Virtual class representing an isotopic distribution."""

    def __init__(
        self,
        formula="",
        get_confs=False,
        atomCounts=None,
        isotopeMasses=None,
        isotopeProbabilities=None,
        use_nominal_masses=False,
        fasta="",
        charge=1.0,
    ):
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
        self.get_confs = get_confs
        self.charge = charge

        if len(fasta) > 0:
            molecule = ParseFASTA(fasta)
        else:
            molecule = OrderedDict()

        if len(formula) > 0:
            if isinstance(formula, dict):
                df = formula
            else:
                df = ParseFormula(formula)
            for elem, count in df.items():
                molecule[elem] = molecule.get(elem, 0) + count

        if atomCounts is not None:
            assert isotopeMasses is not None and isotopeProbabilities is not None
            assert len(atomCounts) == len(isotopeMasses) == len(isotopeProbabilities)

            # Flatten the data structures
            dimNumber = len(atomCounts)
            isotopeNumbers = [len(masses) for masses in isotopeMasses]

            flat_masses = []
            flat_probs = []
            for masses, probs in zip(isotopeMasses, isotopeProbabilities):
                flat_masses.extend(masses)
                flat_probs.extend(probs)

            # Create numpy arrays
            isotope_nums_arr = np.array(isotopeNumbers, dtype=np.int32)
            atom_counts_arr = np.array(atomCounts, dtype=np.int32)
            masses_arr = np.array(flat_masses, dtype=np.float64)
            probs_arr = np.array(flat_probs, dtype=np.float64)

            self.iso = nb.Iso(
                dimNumber, isotope_nums_arr, atom_counts_arr, masses_arr, probs_arr
            )

            self.dimNumber = dimNumber
            self.isotopeNumbers = isotopeNumbers
            self.atomCounts = list(atomCounts)
            self.isotopeMasses = [list(m) for m in isotopeMasses]
            self.isotopeProbabilities = [list(p) for p in isotopeProbabilities]

        elif len(molecule) > 0:
            params = IsoParamsFromDict(molecule, use_nominal_masses=use_nominal_masses)

            dimNumber = len(params.atomCounts)
            isotopeNumbers = [len(masses) for masses in params.masses]

            flat_masses = []
            flat_probs = []
            for masses, probs in zip(params.masses, params.probs):
                flat_masses.extend(masses)
                flat_probs.extend(probs)

            isotope_nums_arr = np.array(isotopeNumbers, dtype=np.int32)
            atom_counts_arr = np.array(params.atomCounts, dtype=np.int32)
            masses_arr = np.array(flat_masses, dtype=np.float64)
            probs_arr = np.array(flat_probs, dtype=np.float64)

            self.iso = nb.Iso(
                dimNumber, isotope_nums_arr, atom_counts_arr, masses_arr, probs_arr
            )

            self.dimNumber = dimNumber
            self.isotopeNumbers = isotopeNumbers
            self.atomCounts = list(params.atomCounts)
            self.isotopeMasses = [list(m) for m in params.masses]
            self.isotopeProbabilities = [list(p) for p in params.probs]

        if self.iso is not None and charge != 1.0:
            # Charge handling would need to be implemented in masses retrieval
            self._charge_divisor = charge
        else:
            self._charge_divisor = 1.0

        # For configuration parsing
        if self.iso is not None:
            offsets = []
            current_offset = 0
            for isoNum in self.isotopeNumbers:
                offsets.append(tuple(range(current_offset, current_offset + isoNum)))
                current_offset += isoNum
            self.offsets = tuple(offsets)

    def __del__(self):
        # Nanobind handles deletion automatically
        self.iso = None

    def getLightestPeakMass(self):
        """Get the mass of the lightest peak in the isotopic distribution."""
        return self.iso.get_lightest_peak_mass() / self._charge_divisor

    def getLightestPeakLProb(self):
        """Get the log probability of the lightest peak in the isotopic distribution."""
        return self.iso.get_lightest_peak_lprob()

    def getLightestPeakConf(self):
        """Get the isotopic configuration of the lightest peak in the isotopic distribution."""
        signature = self.iso.get_lightest_peak_signature()
        return self.parse_conf(signature)

    def getHeaviestPeakMass(self):
        """Get the mass of the heaviest peak in the isotopic distribution."""
        return self.iso.get_heaviest_peak_mass() / self._charge_divisor

    def getHeaviestPeakLProb(self):
        """Get the log probability of the heaviest peak in the isotopic distribution."""
        return self.iso.get_heaviest_peak_lprob()

    def getHeaviestPeakConf(self):
        """Get the isotopic configuration of the heaviest peak in the isotopic distribution."""
        signature = self.iso.get_heaviest_peak_signature()
        return self.parse_conf(signature)

    def getMonoisotopicPeakMass(self):
        """Get the mass of the monoisotopic peak in the isotopic distribution."""
        return self.iso.get_monoisotopic_peak_mass() / self._charge_divisor

    def getMonoisotopicPeakLProb(self):
        """Get the log probability of the monoisotopic peak in the isotopic distribution."""
        return self.iso.get_monoisotopic_peak_lprob()

    def getMonoisotopicPeakConf(self):
        """Get the isotopic configuration of the monoisotopic peak in the isotopic distribution."""
        signature = self.iso.get_monoisotopic_peak_signature()
        return self.parse_conf(signature)

    def getModeLProb(self):
        """Get the log probability of the most probable peak(s) in the isotopic distribution."""
        return self.iso.get_mode_lprob()

    def getModeMass(self):
        """Get the mass of the most probable peak.

        If there are more, return only the mass of one of them."""
        return self.iso.get_mode_mass() / self._charge_divisor

    def getTheoreticalAverageMass(self):
        return self.iso.get_theoretical_average_mass() / self._charge_divisor

    def variance(self):
        return self.iso.variance() / (self._charge_divisor**2)

    def stddev(self):
        return self.iso.stddev() / self._charge_divisor

    def getMarginalLogSizeEstimates(self, prob):
        return self.iso.get_marginal_log_size_estimates(prob)

    def parse_conf(self, cptr, starting_with=0):
        return tuple(tuple(cptr[i + starting_with] for i in o) for o in self.offsets)

    def _get_parse_conf_fun(self):
        # Can't just use the above function as lambda, as the entire class instance will be in closure
        offsets = self.offsets

        def pc(cptr, starting_with=0):
            return tuple(tuple(cptr[i + starting_with] for i in o) for o in offsets)

        return pc


class IsoDistribution(object):
    """Isotopoic distribution with precomputed vector of masses and probabilities."""

    def np_masses(self):
        """Return computed masses as a numpy array."""
        return np.array(self.masses)

    def np_probs(self):
        """Return computed probabilities as a numpy array."""
        return np.array(self.probs)

    def __iter__(self):
        if hasattr(self, "confs") and self.confs is not None:
            for i in range(self.size):
                yield (self.masses[i], self.probs[i], self.confs[i])
        else:
            for i in range(self.size):
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
        return self.parse_conf(
            self.raw_confs, starting_with=sum(self.iso.isotopeNumbers) * idx
        )

    def __del__(self):
        # Nanobind handles cleanup automatically
        pass

    def __init__(
        self, envelope_obj=None, probs=None, masses=None, get_confs=False, iso=None
    ):
        self.mass_sorted = False
        self.prob_sorted = False
        self._total_prob = float("nan")
        self.envelope = envelope_obj

        if envelope_obj is not None:
            self.size = envelope_obj.confs_no()
            self.masses = envelope_obj.masses()
            self.probs = envelope_obj.probs()
            self._total_prob = envelope_obj.get_total_prob()

            if get_confs and iso is not None:
                confs_list = envelope_obj.confs()
                if confs_list:
                    self.raw_confs = [
                        item for sublist in confs_list for item in sublist
                    ]
                    self.sum_isotope_numbers = sum(iso.isotopeNumbers)
                    self.confs = ConfsPassthrough(
                        lambda idx: self._get_conf(idx), self.size
                    )
                    self.parse_conf = iso._get_parse_conf_fun()
                    self.iso = iso

        elif probs is not None or masses is not None:
            assert probs is not None and masses is not None
            assert len(probs) == len(masses)
            self.size = len(probs)
            self.probs = list(probs)
            self.masses = list(masses)

        elif envelope_obj == probs == masses == get_confs == iso == None:
            self.size = 0
            self.probs = []
            self.masses = []
            self._total_prob = 0.0
            self.mass_sorted = True
            self.prob_sorted = True
        else:
            raise RuntimeError("Invalid arguments for IsoDistribution constructor")

    def _get_envelope(self):
        if self.envelope is not None:
            return self.envelope
        masses_arr = np.array(self.masses, dtype=np.float64)
        probs_arr = np.array(self.probs, dtype=np.float64)
        return nb.FixedEnvelope(
            masses_arr, probs_arr, self.mass_sorted, self.prob_sorted, self._total_prob
        )

    def copy(self):
        env = self._get_envelope()
        new_env = nb.FixedEnvelope(env)
        ret = IsoDistribution(envelope_obj=new_env)
        ret._total_prob = self._total_prob
        ret.mass_sorted = self.mass_sorted
        ret.prob_sorted = self.prob_sorted
        return ret

    def __add__(self, other):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        result_env = env1 + env2
        ret = IsoDistribution(envelope_obj=result_env)
        return ret

    def __mul__(self, other):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        result_env = env1 * env2
        ret = IsoDistribution(envelope_obj=result_env)
        return ret

    def total_prob(self):
        if math.isnan(self._total_prob):
            env = self._get_envelope()
            self._total_prob = env.get_total_prob()
        return self._total_prob

    def normalize(self):
        env = self._get_envelope()
        # If total_prob is 0 or invalid, calculate it first
        if env.get_total_prob() <= 0.0 or math.isnan(env.get_total_prob()):
            # Calculate actual total_prob as sum of probabilities
            actual_total = sum(env.probs())
            if actual_total > 0.0:
                # Create a new envelope with correct total_prob
                masses_arr = np.array(env.masses(), dtype=np.float64)
                probs_arr = np.array(env.probs(), dtype=np.float64)
                env = nb.FixedEnvelope(
                    masses_arr,
                    probs_arr,
                    self.mass_sorted,
                    self.prob_sorted,
                    actual_total,
                )
        env.normalize()
        self.masses = env.masses()
        self.probs = env.probs()
        self.envelope = env
        self._total_prob = 1.0

    def normalized(self):
        ret = self.copy()
        ret.normalize()
        return ret

    def add_mass(self, d_mass):
        self.masses = [m + d_mass for m in self.masses]
        self.envelope = None

    def mul_mass(self, d_mass):
        self.masses = [m * d_mass for m in self.masses]
        self.envelope = None

    def add_mul_mass(self, add, mul):
        self.masses = [(m * mul) + (add * mul) for m in self.masses]
        self.envelope = None

    def mul_add_mass(self, mul, add):
        self.masses = [(m * mul) + add for m in self.masses]
        self.envelope = None

    def scale(self, factor):
        """Multiplies the probabilities of spectrum by factor. Works in place."""
        env = self._get_envelope()
        env.scale(factor)
        self.probs = env.probs()
        self.envelope = env
        self._total_prob *= factor

    def scaled(self, factor):
        """Returns a copy of the spectrum where each probability was multiplied by factor."""
        ret = self.copy()
        ret.scale(factor)
        return ret

    def resample(self, ionic_current, beta_bias=1.0):
        env = self._get_envelope()
        env.resample(ionic_current, beta_bias)
        self.probs = env.probs()
        self.envelope = env
        self._total_prob = float(ionic_current)
        self.prob_sorted = False

    def sort_by_prob(self):
        if not self.prob_sorted:
            env = self._get_envelope()
            env.sort_by_prob()
            self.masses = env.masses()
            self.probs = env.probs()
            self.envelope = env
            self.mass_sorted = False
            self.prob_sorted = True

    def sort_by_mass(self):
        if not self.mass_sorted:
            env = self._get_envelope()
            env.sort_by_mass()
            self.masses = env.masses()
            self.probs = env.probs()
            self.envelope = env
            self.mass_sorted = True
            self.prob_sorted = False

    def _recalculate_everything(self):
        self._total_prob = float("nan")
        self.mass_sorted = False
        self.prob_sorted = False
        self.envelope = None

    def empiric_average_mass(self):
        env = self._get_envelope()
        return env.empiric_average_mass()

    def empiric_variance(self):
        env = self._get_envelope()
        return env.empiric_variance()

    def empiric_stddev(self):
        env = self._get_envelope()
        return env.empiric_stddev()

    def wassersteinDistance(self, other):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        ret = env1.wasserstein_distance(env2)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        if math.isnan(ret):
            raise ValueError(
                "Both spectra must be normalized before Wasserstein distance can be computed."
            )
        return ret

    def orientedWassersteinDistance(self, other):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        ret = env1.oriented_wasserstein_distance(env2)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        if math.isnan(ret):
            raise ValueError(
                "Both spectra must be normalized before Wasserstein distance can be computed."
            )
        return ret

    def abyssalWassersteinDistance(self, other, abyss_depth, other_scale=1.0):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        ret = env1.abyssal_wasserstein_distance(env2, abyss_depth, other_scale)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        return ret

    def wassersteinMatch(self, other, flow_dist, other_scale=1.0):
        env1 = self._get_envelope()
        env2 = other._get_envelope()
        ret = env1.wasserstein_match(env2, flow_dist, other_scale)
        self.mass_sorted = True
        self.prob_sorted = False
        other.mass_sorted = True
        other.prob_sorted = False
        return ret

    def binned(self, width=1.0, middle=0.0):
        env = self._get_envelope()
        binned_env = env.bin(width, middle)
        ret = IsoDistribution(envelope_obj=binned_env)
        self.mass_sorted = True
        self.prob_sorted = False
        ret.mass_sorted = True
        ret.prob_sorted = False
        return ret

    @staticmethod
    def LinearCombination(envelopes, intensities):
        envelope_objs = [x._get_envelope() for x in envelopes]
        result_env = nb.FixedEnvelope.linear_combination(
            envelope_objs, list(intensities)
        )
        ret = IsoDistribution(envelope_obj=result_env)
        return ret

    def plot(self, plot_type="bars", show=True, **matplotlib_args):
        """Plot the isotopic distribution.

        Args:
            **matplotlib_args: arguments for matplotlib plot.
        """
        try:
            from matplotlib import pyplot as plt
        except ImportError as e:
            raise ImportError(
                str(e) + "\nPlotting spectra requires matplotlib to be installed."
            )
        if plot_type == "bars":
            if "linewidth" not in matplotlib_args:
                matplotlib_args["linewidth"] = 1.0
            plt.vlines(list(self.masses), [0], list(self.probs), **matplotlib_args)
        elif plot_type == "profile":
            self.sort_by_mass()
            plt.plot(list(self.masses), list(self.probs), **matplotlib_args)
        plt.xlabel("Mass (Da)")
        plt.ylabel("Intensity (relative)")
        if show:
            plt.show()


def IsoThreshold(threshold, formula="", absolute=False, get_confs=False, **kwargs):
    """Initialize the IsoDistribution isotopic distribution by threshold.

    Args:
        threshold (float): value of the absolute or relative threshold.
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
        absolute (boolean): should we report peaks with probabilities above an absolute probability threshold, or above a relative threshold amounting to a given proportion of the most probable peak?
        get_confs (boolean): should we report counts of isotopologues?
        **kwds: named arguments to IsoSpectrum.
    """
    iso = Iso(formula=formula, get_confs=get_confs, **kwargs)
    envelope = nb.FixedEnvelope.from_threshold(iso.iso, threshold, absolute, get_confs)
    ido = IsoDistribution(envelope_obj=envelope, get_confs=get_confs, iso=iso)
    return ido


def IsoTotalProb(
    prob_to_cover, formula="", get_minimal_pset=True, get_confs=False, **kwargs
):
    """Initialize the IsoDistribution isotopic distribution by total probability.

    Args:
        prob_to_cover (float): minimal total probability of the reported peaks.
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
        get_minimal_pset (boolean): should we trim the last calculated layer of isotopologues so that the reported result is as small as possible?
        get_confs (boolean): should we report the counts of isotopologues?
        **kwargs: named arguments to the superclass.
    """
    iso = Iso(formula=formula, get_confs=get_confs, **kwargs)
    envelope = nb.FixedEnvelope.from_total_prob(
        iso.iso, prob_to_cover, get_minimal_pset, get_confs
    )
    ido = IsoDistribution(envelope_obj=envelope, get_confs=get_confs, iso=iso)
    return ido


def IsoStochastic(
    no_molecules, formula="", precision=0.9999, beta_bias=5.0, get_confs=False, **kwargs
):
    """Initialize the IsoDistribution isotopic distribution by stochastic sampling.

    Args:
        no_molecules (uint): ionic current in instrument
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
        precision (float): passed to IsoTotalProbGenerator. Between 0.0 and 1.0.
        beta_bias (float, nonnegative): fiddling with this parameter does not change the result, but might make computations slightly faster (or likely, much, much slower is you screw it up...)
        get_confs (boolean): should we report the counts of isotopologues?
        **kwargs: named arguments to the superclass.
    """
    iso = Iso(formula=formula, get_confs=get_confs, **kwargs)
    envelope = nb.FixedEnvelope.from_stochastic(
        iso.iso, no_molecules, precision, beta_bias, get_confs
    )
    ido = IsoDistribution(envelope_obj=envelope, get_confs=get_confs, iso=iso)
    return ido


def IsoBinned(
    bin_width, formula="", target_total_prob=0.9999, bin_middle=0.0, **kwargs
):
    """Initialize the IsoDistribution isotopic distribution by binning.

    Args:
        bin_width (float): width of bins
        formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
        target_total_prob (float): target total probability coverage
        bin_middle (float): position of bin centers
        **kwargs: named arguments to the superclass.
    """
    iso = Iso(formula=formula, get_confs=False, **kwargs)
    envelope = nb.FixedEnvelope.binned(
        iso.iso, target_total_prob, bin_width, bin_middle
    )
    ido = IsoDistribution(envelope_obj=envelope, get_confs=False, iso=iso)
    return ido


class IsoGenerator(Iso):
    """Virtual class allowing memory-efficient iteration over the isotopic distribution.

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
        super(IsoGenerator, self).__init__(
            formula=formula, get_confs=get_confs, **kwargs
        )
        self.firstuse = True

    def __iter__(self):
        if not self.firstuse:
            raise NotImplementedError(
                "Multiple iterations through the same IsoGenerator object are not supported. Either create a new (identical) generator for a second loop-through, or use one of the non-generator classes, which do support being re-used."
            )
        self.firstuse = False
        cgen = self.cgen
        if self.get_confs:
            while cgen.advance_to_next_configuration():
                conf = cgen.get_conf_signature()
                yield (
                    cgen.mass() / self._charge_divisor,
                    cgen.prob(),
                    self.parse_conf(conf),
                )
        else:
            while cgen.advance_to_next_configuration():
                yield (cgen.mass() / self._charge_divisor, cgen.prob())

    def __del__(self):
        super(IsoGenerator, self).__del__()
        self.cgen = None


class IsoThresholdGenerator(IsoGenerator):
    """Class allowing memory-efficient iteration over the isotopic distribution up till some probability threshold.

    This iterator will stop only after reaching a probability threshold.
    """

    def __init__(
        self,
        threshold,
        formula="",
        absolute=False,
        get_confs=False,
        reorder_marginals=True,
        **kwargs,
    ):
        """Initialize IsoThresholdGenerator.

        Args:
            threshold (float): threshold value
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            absolute (boolean): should we report peaks with probabilities above an absolute probability threshold, or above a relative threshold amounting to a given proportion of the most probable peak?
            get_confs (boolean): should we report the counts of isotopologues?
            reorder_marginals (boolean): should marginals be reordered for efficiency?
            **kwargs: named arguments to the superclass.
        """
        super(IsoThresholdGenerator, self).__init__(
            formula=formula, get_confs=get_confs, **kwargs
        )
        self.threshold = threshold
        self.absolute = absolute
        self.cgen = nb.IsoThresholdGenerator(
            self.iso, threshold, absolute, 1000, 1000, reorder_marginals
        )

    def __del__(self):
        """Destructor."""
        self.cgen = None
        super(IsoThresholdGenerator, self).__del__()


class IsoLayeredGenerator(IsoGenerator):
    """Class allowing memory-efficient iteration over the isotopic distribution up till some joint probability of the reported peaks."""

    def __init__(
        self,
        formula="",
        get_confs=False,
        reorder_marginals=True,
        t_prob_hint=0.99,
        **kwargs,
    ):
        """Initialize IsoLayeredGenerator.

        Args:
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            get_confs (boolean): should we report the counts of isotopologues?
            reorder_marginals (boolean): should marginals be reordered for efficiency?
            t_prob_hint (float): hint for target probability
            **kwargs: named arguments to the superclass.
        """
        super(IsoLayeredGenerator, self).__init__(
            formula=formula, get_confs=get_confs, **kwargs
        )
        self.cgen = nb.IsoLayeredGenerator(
            self.iso, 1000, 1000, reorder_marginals, t_prob_hint
        )

    def __del__(self):
        self.cgen = None
        super(IsoLayeredGenerator, self).__del__()


class IsoOrderedGenerator(IsoGenerator):
    """Class representing an isotopic distribution with peaks ordered by descending probability.

    This generator return probabilities ordered with descending probability.
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
        super(IsoOrderedGenerator, self).__init__(
            formula=formula, get_confs=get_confs, **kwargs
        )
        self.cgen = nb.IsoOrderedGenerator(self.iso, 1000, 1000)

    def __del__(self):
        self.cgen = None
        super(IsoOrderedGenerator, self).__del__()


class IsoStochasticGenerator(IsoGenerator):
    """Class for stochastic generation of isotopic distributions."""

    def __init__(
        self,
        no_molecules,
        formula="",
        precision=0.9999,
        beta_bias=1.0,
        get_confs=False,
        **kwargs,
    ):
        """Initialize IsoStochasticGenerator.

        Args:
            no_molecules (int): number of molecules
            formula (str): a chemical formula, e.g. "C2H6O1" or "C2H6O".
            precision (float): precision parameter
            beta_bias (float): beta bias parameter
            get_confs (boolean): should we report the counts of isotopologues?
            **kwargs: named arguments to the superclass.
        """
        super(IsoStochasticGenerator, self).__init__(
            formula=formula, get_confs=get_confs, **kwargs
        )
        self.threshold = precision
        self.no_molecules = no_molecules
        self.cgen = nb.IsoStochasticGenerator(
            self.iso, no_molecules, precision, beta_bias
        )

    def __del__(self):
        """Destructor."""
        self.cgen = None
        super(IsoStochasticGenerator, self).__del__()
