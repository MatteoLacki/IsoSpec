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
                                       int _tabSize,
                                       int _hashSize
                                       );
        double massIsoLayeredGenerator(void* generator); double lprobIsoLayeredGenerator(void* generator); double probIsoLayeredGenerator(void* generator); void methodIsoLayeredGenerator(void* generator); bool advanceToNextConfigurationIsoLayeredGenerator(void* generator); void deleteIsoLayeredGenerator(void* generator); void get_conf_signatureIsoLayeredGenerator(void* generator, int* space);


        void* setupIsoOrderedGenerator(void* iso,
                                       int _tabSize,
                                       int _hashSize);
        double massIsoOrderedGenerator(void* generator); double lprobIsoOrderedGenerator(void* generator); double probIsoOrderedGenerator(void* generator); void methodIsoOrderedGenerator(void* generator); bool advanceToNextConfigurationIsoOrderedGenerator(void* generator); void deleteIsoOrderedGenerator(void* generator); void get_conf_signatureIsoOrderedGenerator(void* generator, int* space);

        void* setupThresholdFixedEnvelope(void* iso,
                                    double threshold,
                                    bool absolute,
                                    bool get_confs,
                                    bool get_lprobs,
                                    bool get_masses,
                                    bool get_probs);

        void deleteThresholdFixedEnvelope(void* tabulator);

        const double* massesThresholdFixedEnvelope(void* tabulator);
        const double* lprobsThresholdFixedEnvelope(void* tabulator);
        const double* probsThresholdFixedEnvelope(void* tabulator);
        const int* confsThresholdFixedEnvelope(void* tabulator);
        int confs_noThresholdFixedEnvelope(void* tabulator);


        void* setupTotalProbFixedEnvelope(void* iso,
                                      double taget_coverage,
                                      bool optimize,
                                      bool get_confs,
                                      bool get_lprobs,
                                      bool get_masses,
                                      bool get_probs);

        void deleteTotalProbFixedEnvelope(void* tabulator);

        const double* massesTotalProbFixedEnvelope(void* tabulator);
        const double* lprobsTotalProbFixedEnvelope(void* tabulator);
        const double* probsTotalProbFixedEnvelope(void* tabulator);
        const int* confsTotalProbFixedEnvelope(void* tabulator);
        int confs_noTotalProbFixedEnvelope(void* tabulator);

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
