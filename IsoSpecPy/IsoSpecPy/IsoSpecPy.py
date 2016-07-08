# -*- coding: utf-8 -*-
#
#    Copyright (C) 2015 Mateusz £±cki and Micha³ Startek.
#
#    This file is part of IsoSpec.
#
#    IsoSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License
#    version 3, as published by the Free Software Foundation.
#
#    IsoSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with IsoSpec.  If not, see <http://www.gnu.org/licenses/>.
#

import cffi
import itertools
import math
import re
import os
import glob
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from kahan import Summator
from collections import defaultdict

try:
    xrange
except NameError:
    xrange = range

class IsoFFI:
    def __init__(self):
        self.ffi = cffi.FFI()
        self.ffi.cdef('''void* setupMarginal(
                                const double* masses,
                                const double* probs,
                                int isotopeNo,
                                int atomCnt,
                                const int tabSize,
                                const int hashSize
                                );
                        int probeConfigurationIdx(void* MT, int idx);
                        int getConfNo(void* marginals);
                        int getIsotopesNo(void* iso);

                        void getConfs(
                                        int howmany,
                                        void* marginals,
                                        double* masses,
                                        double* logprobs,
                                        int* configurations
                                      );
                        void destroyConf(void* marginals);

                        void* setupIso( int             _dimNumber,
                                        const int*      _isotopeNumbers,
                                        const int*      _atomCounts,
                                        const double*   _isotopeMasses,
                                        const double*   _isotopeProbabilities,
                                        const double    _StopCondition,
                                        int             algo,
                                        int             tabSize,
                                        int             hashSize,
                                        double          step
                        );


                        void* IsoFromFormula(const char* formula, double cutoff, int tabsize, int hashsize);

                        int processMTUntilCutoff(void* MT, double cutoff);

                        int getConfMT(void* MT, int idx, double* mass, double* logProb, int* configuration);

                        int getIsoConfNo(void* iso);

                        void getIsoConfs(   void* iso,
                                            double* res_mass,
                                            double* res_logProb,
                                            int* res_isoCounts
                                        );

                        void destroyIso(void* iso);


                        #define NUMBER_OF_ISOTOPIC_ENTRIES 288

                        extern const int elem_table_atomicNo[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const double elem_table_probability[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const double elem_table_mass[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const int elem_table_massNo[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const int elem_table_extraNeutrons[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const char* elem_table_element[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const char* elem_table_symbol[NUMBER_OF_ISOTOPIC_ENTRIES];
                        extern const bool elem_table_Radioactive[NUMBER_OF_ISOTOPIC_ENTRIES];


                        ''');

        mod_dir = os.path.dirname(os.path.abspath(__file__))

        if os.path.exists(os.path.join(mod_dir, '..', 'setup.py')):
            raise Exception('''Attempted to load IsoSpecPy module from its build directory. This usually
won't work and is generally a Bad Idea. Please cd somewhere else, or, if you're really
sure you want to do that, edit the source and disable this check.''')

        libpaths = [
            'IsoSpecCppPy*.so',
            os.path.join(mod_dir, 'IsoSpecCppPy*.so'),
            os.path.join(mod_dir, '..', 'IsoSpecCppPy*.so'),
            'libIsoSpec++*.so',
            os.path.join(mod_dir, 'libIsoSpec++*.so'),
            os.path.join(mod_dir, '..', 'libIsoSpec++*.so'),
            os.path.join('..', 'IsoSpec++', 'libIsoSpec++*.so')
        ]

        self.clib = None
        for libpath in libpaths:
            try:
                fn = glob.glob(libpath)[0]
                self.clib = self.ffi.dlopen(fn)
                break
            except (OSError, IndexError):
                pass

        if self.clib == None:
            raise Exception("Cannot find IsoSpecCppPy.so")



isoFFI = IsoFFI()



class MarginalDistribution:
    def __init__(
                    self,
                    masses,
                    probs,
                    atomCnt,
                    tabSize = 1000,
                    hashSize = 1000,
                    burnIn = 1000
                ):
        self.clib = isoFFI.clib #because you can't use global vars in destructor, which is UTTERLY STUPID
        self.masses = masses
        self.probs = probs
        self.isotopeNo = len(masses)
        self.atomCnt = atomCnt
        self.tabSize = tabSize
        self.hashSize = hashSize
        self.marginal = isoFFI.clib.setupMarginal(
                                masses,
                                probs,
                                self.isotopeNo,
                                atomCnt,
                                tabSize,
                                hashSize)
        self.cconf_mass = isoFFI.ffi.new("double[1]")
        self.cconf_lProb = isoFFI.ffi.new("double[1]")
        self.cconf = isoFFI.ffi.new("int[{0}]".format(self.isotopeNo))
        self.cprobs = []
        self.summator = Summator()


    def __del__(self):
        self.clib.destroyConf(self.marginal)

    def getConfsRaw(self):
        masses = isoFFI.ffi.new("double[{0}]".format(len(self)))
        logProbs = isoFFI.ffi.new("double[{0}]".format(len(self)))
        confs = isoFFI.ffi.new("int[{0}]".format(len(self)*self.isotopeNo))
        isoFFI.clib.getConfs(len(self), self.marginal, masses, logProbs, confs)
        return (masses, logProbs, confs)

    def getConfs(self):
        ret = []
        masses, logProbs, confs = self.getConfsRaw()
        for i in xrange(len(masses)):
            ret.append((masses[i], logProbs[i], list(confs[self.isotopeNo*i:(i+1)*self.isotopeNo])))
        return ret

    def getConf(self, idx):
        if isoFFI.clib.getConfMT(self.marginal, idx, self.cconf_mass, self.cconf_lProb, self.cconf) == 1:
            return (self.cconf_mass[0], self.cconf_lProb[0], tuple(self.cconf))
        else:
            raise IndexError

    def probeConf(self, idx):
        return self.clib.probeConfigurationIdx(self.marginal, idx) > 0

    def CumulativeProb(self, idx):
        if len(self.cprobs) > idx:
            return self.cprobs[idx]
        else:
            for i in xrange(len(self.cprobs), idx+1):
                self.summator.add(math.exp(self.getConf(i)[1]))
                self.cprobs.append(self.summator.get())
            return self.cprobs[idx]



class IsoSpec:
    def __init__(
                    self,
                    _atomCounts,
                    _isotopeMasses,
                    _isotopeProbabilities,
                    _stopCondition,
                    tabSize = 1000,
                    hashSize = 1000,
                    step = 0.3,
                    method = 'layered'
                ):
        self.clib = isoFFI.clib #can't use global vars in destructor, again...
        self.dimNumber                 = len(_atomCounts)
        self._isotopeNumbers           = [len(x) for x in _isotopeMasses]
        self.allIsotopeNumber         = sum(self._isotopeNumbers)
        self._atomCounts               = _atomCounts
        self._isotopeMasses            = _isotopeMasses
        self._isotopeProbabilities     = _isotopeProbabilities
        self._stopCondition            = _stopCondition
        self.tabSize                   = tabSize
        self.hashSize                  = hashSize

        try:
            self.algo = { 'layered' : 0,
              'ordered' : 1,
              'threshold_absolute' : 2,
              'threshold_relative' : 3,
              'layered_estimating' : 4,
            }[method]
        except KeyError:
            raise Exception("Invalid ISO method")

        self.iso = isoFFI.clib.setupIso(
                                self.dimNumber,
                                self._isotopeNumbers,
                                _atomCounts,
                                list(itertools.chain.from_iterable(_isotopeMasses)),
                                list(itertools.chain.from_iterable(_isotopeProbabilities)),
                                _stopCondition,
                                self.algo,
                                tabSize,
                                hashSize,
                                step
                            )
        try:
            if not self.iso.__nonzero__():
                raise MemoryError()
        except AttributeError: # Python3 doesn't have __nonzero__...
            pass



    @staticmethod
    def IsoFromFormula(formula, cutoff, tabSize = 1000, hashSize = 1000, classId = None, method = 'layered', step = 0.25):
        # It's much easier to just parse it in python than to use the C parsing function
        # and retrieve back into Python the relevant object sizes
        symbols = re.findall("\D+", formula)
        atom_counts = [int(x) for x in re.findall("\d+", formula)]

        if not len(symbols) == len(atom_counts):
            raise ValueError("Invalid formula")

        indexes = [[x for x in xrange(isoFFI.clib.NUMBER_OF_ISOTOPIC_ENTRIES)
                        if isoFFI.ffi.string(isoFFI.clib.elem_table_symbol[x]) == symbol.encode('latin1')]
                    for symbol in symbols]

        if any([len(x) == 0 for x in indexes]):
            raise ValueError("Invalid formula")

        masses  = [[isoFFI.clib.elem_table_mass[idx] for idx in idxs] for idxs in indexes]
        probs   = [[isoFFI.clib.elem_table_probability[idx] for idx in idxs] for idxs in indexes]

        if classId == None:
            return IsoSpec(atom_counts, masses, probs, cutoff, tabSize, hashSize, step, method)
        else:
            return classId(atom_counts, masses, probs, cutoff, tabSize, hashSize)



    def __del__(self):
        self.cleanup()

    def cleanup(self):
        if self.iso is not None:
            self.clib.destroyIso(self.iso)
            self.iso = None

    def __len__(self):
        return isoFFI.clib.getIsoConfNo(self.iso)

    def getConfsRaw(self):
        masses = isoFFI.ffi.new("double[{0}]".format(len(self)))
        logProbs = isoFFI.ffi.new("double[{0}]".format(len(self)))
        isoCounts = isoFFI.ffi.new("int[{0}]".format(len(self)*sum(self._isotopeNumbers)))
        isoFFI.clib.getIsoConfs(self.iso, masses, logProbs, isoCounts)
        return (masses, logProbs, isoCounts)

    def get_conf_by_no(self, clist, idx):
        idx *= self.allIsotopeNumber
        ret = []
        for ison in self._isotopeNumbers:
            ret.append(tuple(clist[idx:idx+ison]))
            idx += ison
        return tuple(ret)


    def getConfs(self):
        masses, logProbs, isoCounts = self.getConfsRaw()
        return [(
                masses[i],
                logProbs[i],
                self.get_conf_by_no(isoCounts, i))
                for i in xrange(len(masses))]


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




