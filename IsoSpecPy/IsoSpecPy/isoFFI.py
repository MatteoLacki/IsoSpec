import cffi
import os
import platform
import sys
import glob


class IsoFFI:
    def __init__(self):
        self.ffi = cffi.FFI()
        self.ffi.cdef('''
        void * setupIso(int dimNumber,
                const int* isotopeNumbers,
                const int* atomCounts,
                const double* isotopeMasses,
                const double* isotopeProbabilities);

        double getLightestPeakMassIso(void* iso);
        double getHeaviestPeakMassIso(void* iso);
        double getMonoisotopicPeakMassIso(void* iso);
        double getModeLProbIso(void* iso);
        double getModeMassIso(void* iso);
        double getTheoreticalAverageMassIso(void* iso);

        void deleteIso(void* iso);

        void* setupIsoThresholdGenerator(void* iso,
                                         double threshold,
                                         bool _absolute,
                                         int _tabSize,
                                         int _hashSize);
        double massIsoThresholdGenerator(void* generator); double lprobIsoThresholdGenerator(void* generator); double probIsoThresholdGenerator(void* generator); void methodIsoThresholdGenerator(void* generator); bool advanceToNextConfigurationIsoThresholdGenerator(void* generator); void deleteIsoThresholdGenerator(void* generator); void get_conf_signatureIsoThresholdGenerator(void* generator, int* space);



        void* setupIsoLayeredGenerator(void* iso,
                                       double _target_coverage,
                                       double _percentage_to_expand,
                                       int _tabSize,
                                       int _hashSize,
                                       bool _do_trim);
        double massIsoLayeredGenerator(void* generator); double lprobIsoLayeredGenerator(void* generator); double probIsoLayeredGenerator(void* generator); void methodIsoLayeredGenerator(void* generator); bool advanceToNextConfigurationIsoLayeredGenerator(void* generator); void deleteIsoLayeredGenerator(void* generator); void get_conf_signatureIsoLayeredGenerator(void* generator, int* space);


        void* setupIsoOrderedGenerator(void* iso,
                                       int _tabSize,
                                       int _hashSize);
        double massIsoOrderedGenerator(void* generator); double lprobIsoOrderedGenerator(void* generator); double probIsoOrderedGenerator(void* generator); void methodIsoOrderedGenerator(void* generator); bool advanceToNextConfigurationIsoOrderedGenerator(void* generator); void deleteIsoOrderedGenerator(void* generator); void get_conf_signatureIsoOrderedGenerator(void* generator, int* space);

        void* setupThresholdTabulator(void* generator,
                                      bool get_masses,
                                      bool get_probs,
                                      bool get_lprobs,
                                      bool get_confs);

        void deleteThresholdTabulator(void* tabulator);

        const double* massesThresholdTabulator(void* tabulator);
        const double* lprobsThresholdTabulator(void* tabulator);
        const double* probsThresholdTabulator(void* tabulator);
        const int* confsThresholdTabulator(void* tabulator);
        int confs_noThresholdTabulator(void* tabulator);


        void* setupLayeredTabulator(void* generator,
                                      bool get_masses,
                                      bool get_probs,
                                      bool get_lprobs,
                                      bool get_confs);

        void deleteLayeredTabulator(void* tabulator);

        const double* massesLayeredTabulator(void* tabulator);
        const double* lprobsLayeredTabulator(void* tabulator);
        const double* probsLayeredTabulator(void* tabulator);
        const int* confsLayeredTabulator(void* tabulator);
        int confs_noLayeredTabulator(void* tabulator);

        void freeReleasedArray(void* array);

        #define NUMBER_OF_ISOTOPIC_ENTRIES 287
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

        libnames  = ['IsoSpecCppPy*', 'IsoSpec++*']
        libprefix = ['', 'lib', 'Lib']
        extension = ['.so', '.dylib', '.dll']
        try:
            if platform.system() == 'Linux':
                extension = ['.so']
            elif platform.system() == 'Windows':
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
