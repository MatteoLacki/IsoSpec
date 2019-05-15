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


'''
Bunch of deprecated functions for 1.0.X compatibility.
Avoid using them: there is a considerable overhead associated
with using the old interface... The backward compatibility module
is also very rudimentary, somewhat incomplete and not very well
tested too...

The current API is implemented in __init__.py, use that instead
'''

try:
    xrange
except NameError:
    xrange = range

import re
from .IsoSpecPy import IsoTotalProb, IsoThreshold

class IsoSpec():
    def __init__(
                    self,
                    _atomCounts,
                    _isotopeMasses,
                    _isotopeProbabilities,
                    _stopCondition,
                    tabSize = 1000,     # ignored
                    hashSize = 1000,    # ignored
                    step = 0.3,         # ignored
                    trim = True,        # True not supported yet, treated as False anyway
                    method = 'layered'
                ):

        isoargs = {
            "formula"           : None,
            "get_confs"         : True,
            "atomCounts"        : _atomCounts,
            "isotopeMasses"     : _isotopeMasses,
            "isotopeProbabilities"  : _isotopeProbabilities,
        }


        self.dimNumber                 = len(_atomCounts)
        self._isotopeNumbers           = [len(x) for x in _isotopeMasses]
        self.allIsotopeNumber         = sum(self._isotopeNumbers)
        self._atomCounts               = _atomCounts
        self._isotopeMasses            = _isotopeMasses
        self._isotopeProbabilities     = _isotopeProbabilities
        self._stopCondition            = _stopCondition

        from .__init__ import IsoThreshold

        try:
            algo = { 'layered' : lambda total_prob: IsoTotalProb(total_prob, **isoargs),
              'ordered' : lambda total_prob: IsoTotalProb(total_prob, **isoargs),
              'threshold_absolute' : lambda threshold: IsoThreshold(threshold, True, **isoargs),
              'threshold_relative' : lambda threshold: IsoThreshold(threshold, False, **isoargs),
              'layered_estimating' : lambda total_prob: IsoTotalProb(total_prob, **isoargs)
            }[method]
        except KeyError:
            raise Exception("Invalid ISO method")

        # Reference to iso needs to be held in this object: it will deallocate masses/lprobs/etc arrays on C++ side if we 
        # allow GC to collect it prematurely
        self.iso = algo(_stopCondition)

        self.masses = self.iso.masses
        self.lprobs = self.iso.lprobs
        self.probs  = self.iso.probs
        self.confs  = self.iso.confs
        self.size   = self.iso.size

        if method == 'ordered' and len(self.masses) > 0:
            L = sorted(zip(self.masses, self.lprobs, self.probs, self.confs), key = lambda x: -x[1])
            self.masses, self.lprobs, self.probs, self.confs = zip(*L)

    @staticmethod
    def IsoFromFormula(formula, cutoff, tabSize = 1000, hashSize = 1000, classId = None, method = 'threshold_relative', step = 0.25, trim = True):
        # It's much easier to just parse it in python than to use the C parsing function
        # and retrieve back into Python the relevant object sizes
        symbols = re.findall("\D+", formula)
        atom_counts = [int(x) for x in re.findall("\d+", formula)]

        if not len(symbols) == len(atom_counts):
            raise ValueError("Invalid formula")

        from .PeriodicTbl import symbol_to_masses, symbol_to_probs
        try:
            masses = tuple(symbol_to_masses[symbol] for symbol in symbols)
            probs = tuple(symbol_to_probs[symbol] for symbol in symbols)
        except KeyError:
            raise ValueError("Invalid formula")

        return IsoSpec(atom_counts, masses, probs, cutoff, tabSize, hashSize, step, trim, method)



    def __len__(self):
        return self.size

    def getConfsRaw(self):
        return (self.masses, self.lprobs, self.confs)

#    def get_conf_by_no(self, clist, idx):
#        idx *= self.allIsotopeNumber
#        ret = []
#        for ison in self._isotopeNumbers:
#            ret.append(tuple(clist[idx:idx+ison]))
#            idx += ison
#        return tuple(ret)


    def getConfs(self):
        masses, logProbs, isoCounts = self.getConfsRaw()
        rows_no = len(masses)
        masses  = list(masses)
        logProbs= list(logProbs)
        confs = []
        for i in xrange(rows_no):
            confs.append([x for sublist in isoCounts[i] for x in sublist])
        return masses, logProbs, confs


    def getConfsNumpy(self):
        import numpy as np
        masses, logProbs, configurations = self.getConfsRaw()
        rows_no = len(masses)
        if rows_no == 0:
            return np.array([])
        cols_no = len(configurations)//len(masses)
        masses  = np.array(list(masses))
        logProbs= np.array(list(logProbs))
        configurations = np.array(list(configurations)).reshape((rows_no,cols_no))
        return masses, logProbs, configurations


    def splitConf(self, l, offset = 0):
        conf = []
        idx = self.allIsotopeNumber * offset
        for i in xrange(self.dimNumber):
            conf.append(tuple(l[idx:idx+self._isotopeNumbers[i]]))
            idx += self._isotopeNumbers[i]
        return tuple(conf)

    def confStr(self, conf):
        return '\t'.join([' '.join([str(x) for x in y]) for y in conf])

    def printConfs(self):
        masses, logProbs, isoCounts = self.getConfsRaw()
        confs = []
        step = sum(self._isotopeNumbers)
        for i in xrange(len(masses)):
            confs.append((masses[i], logProbs[i], self.splitConf(isoCounts, i)))

        for conf in confs:
            print(("Mass = {0}\t and log-prob = {1} and prob = {2}\t and configuration"\
                  "=\t{3}").format(conf[0], conf[1], math.exp(conf[1]), self.confStr(conf[2])))



class IsoPlot(dict):
    def __init__(self, iso, bin_w):
        self.iso = iso
        self.bin_w = bin_w
        masses, logProbs, _isoCounts = iso.getConfsRaw()
        dd = defaultdict(Summator)
        for i in xrange(len(masses)):
            dd[float(int(masses[i]/bin_w))*bin_w].add(math.exp(logProbs[i]))
        for key, val in dd.items():
            self[key] = val.get()


def IsoSpecify( formula,
                cutoff,
                method= 'layered',
                output_format = 'numpy_arrays',
                trim  = True,
                _step = 0.25,
                _trim = True,
                _tabSize  = 1000,
                _hashSize = 1000    ):
    """
    Call IsoSpec on a formula with a given cutoff.

    This function wraps around the IsoSpec class.

    Parameters
    ----------
    formula : char
        a string of a form '< Element Tag 1 >< Count 1 > ... ',
        e.g. 'C100H202'. Using IUPAC conventions to name elements.

    cutoff : float
        The cutoff value. See description of the method argument.

    method : char
        Can take one of the following values: 'layered',
        'layered_estimating', 'threshold_absolute',
        'threshold_relative', 'ordered'.

        The threshold versions of the algorithm rely on user
        providing a precise lower bound on the reported peak heights.
        This can be specified in absolute terms ('threshold_absolute'),
        i.e. in terms of the limiting probability of the isotopologue,
        or as a percentage of the heighest peak ('threshold_relative').

        The layered versions of the algorithm rely on calculating
        consecutive values of peak thresholds on flight.
        The ultimate goal is to reach a peak probability that assures
        that the sum of probabilities of the more probable isotopologues
        exceeds the provided cutoff value.
        The sequence of consecutive thresholds can be generated in
        two ways. The default way, 'layered_estimating', estimates
        the threshold to joint probability function by a progressive
        linear spline, check Anal Chem. 2017 Mar 21;89(6):3272-3277.
        doi: 10.1021/acs.analchem.6b01459. Epub 2017 Mar 8.
        The other way, 'layered', estimates a threshold as a 30%%
        quantile of the probabilities gathered in the fringe set, i.e.
        isotopologues that are direct neighbours of the previously
        accepted layer. Finally, choosing the 'ordered' version will
        provide a loglinear version of the algorithm that relies on
        the priority queue. This version automatically sorts
        the configurations by their probability.

    trim
        while using a layered method, should one discard superfluously
        obtained isotopologues, i.e. such that without them the set of
        reported isotopologues already is an optimal p-set.

    output_format
        Should the output be presented as lists ('lists'),
        or as numpy arrays ('numpy_arrays').

    Returns
    -------
    masses
        masses of isotopologues, either a list or a numpy array.

    logProbs
        logarithms of probabilities (theoretical heights) of isotopologues,
        either a list or a numpy array.

    confs
        counts of isotopologues (extended chemical formulas that
        include counts of isotopes of elements)
        """


    assert output_format in ('lists', 'numpy_arrays'), "Wrong value of output_format. Should be either 'lists' or 'numpy_arrays'."

    assert method in ('layered', 'ordered', 'threshold_absolute', 'threshold_relative', 'layered_estimating'), "Wrong value of method. Should be among 'layered', 'ordered', 'threshold_absolute', 'threshold_relative', or 'layered_estimating'."

    assert isinstance(cutoff, float), "Provided cut off ain't a float."

    assert isinstance(formula, str), "Provided formula off ain't a string."

    iso = IsoSpec.IsoFromFormula(   formula,
                                    cutoff,
                                    tabSize = 1000,
                                    hashSize = 1000,
                                    classId = None,
                                    method = method,
                                    step = 0.25,
                                    trim = trim )

    if output_format == 'lists':
        masses, logProbs, confs = iso.getConfs()
    else:
        masses, logProbs, confs = iso.getConfsNumpy()

    # print 'Rev Startek is a silly old chump and his mother dresses up silly.'
    return masses, logProbs, confs


