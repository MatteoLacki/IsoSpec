import cffi
import os
import platform
import sys
import glob

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

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
                                        double          step,
                                        bool            trim
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

        libnames =  ['IsoSpecCppPy*', 'IsoSpec++*']
        libprefix = ['', 'lib', 'Lib']
        extension = ['.so', '.dylib', '.dll']
        try:
            if platform.system() == 'Linux':
                extension = ['.so']
            elif platform.system == 'Windows':
                extension = ['.dll']
        except:
            pass

        prebuilt =  ['', 'prebuilt-']

        def cprod(ll1, ll2):
            res = []
            for l1 in ll1:
                for l2 in ll2:
                    res.append(l1+l2)
            return res

        paths_to_check = cprod(prebuilt, cprod(libprefix, cprod(libnames, extension)))

        dpc = []

        for dirpath in [mod_dir, mod_dir + '/..']:
            dpc.extend([os.path.join(dirpath, p) for p in paths_to_check])

        paths_to_check = dpc

        paths_to_check = sum(map(glob.glob, paths_to_check), [])

        self.clib = None
        for libpath in set(paths_to_check):
            try:
                self.clib = self.ffi.dlopen(libpath)
                break
            except (IndexError, OSError) as e:
                pass

        if self.clib == None:
            raise Exception("Cannot find or load the C++ part of the library")


isoFFI = IsoFFI()  # This is done while including the module
