import cffi
import os
import platform
import sys
import glob
if False:
    import IsoSpecCppPy


class IsoFFI:
    def __init__(self):
        self.ffi = cffi.FFI()
        self.ffi.cdef('''
        void * setupIso(int dimNumber,
                const int* isotopeNumbers,
                const int* atomCounts,
                const double* isotopeMasses,
                const double* isotopeProbabilities);

        void * isoFromFasta(const char* fasta, bool use_nominal_masses, bool add_water);

        double getLightestPeakMassIso(void* iso);
        double getHeaviestPeakMassIso(void* iso);
        double getMonoisotopicPeakMassIso(void* iso);
        double getModeLProbIso(void* iso);
        double getModeMassIso(void* iso);
        double getTheoreticalAverageMassIso(void* iso);
        double getIsoVariance(void* iso);
        double getIsoStddev(void* iso);
        double* getMarginalLogSizeEstimates(void* iso, double target_total_prob);

        void deleteIso(void* iso);

        void* setupIsoThresholdGenerator(void* iso,
                                         double threshold,
                                         bool _absolute,
                                         int _tabSize,
                                         int _hashSize,
                                         bool reorder_marginals);
        double massIsoThresholdGenerator(void* generator); double lprobIsoThresholdGenerator(void* generator); double probIsoThresholdGenerator(void* generator); void methodIsoThresholdGenerator(void* generator); bool advanceToNextConfigurationIsoThresholdGenerator(void* generator); void deleteIsoThresholdGenerator(void* generator); void get_conf_signatureIsoThresholdGenerator(void* generator, int* space);



        void* setupIsoLayeredGenerator(void* iso,
                                       int _tabSize,
                                       int _hashSize,
                                       bool reorder_marginals,
                                       double t_prob_hint);
        double massIsoLayeredGenerator(void* generator); double lprobIsoLayeredGenerator(void* generator); double probIsoLayeredGenerator(void* generator); void methodIsoLayeredGenerator(void* generator); bool advanceToNextConfigurationIsoLayeredGenerator(void* generator); void deleteIsoLayeredGenerator(void* generator); void get_conf_signatureIsoLayeredGenerator(void* generator, int* space);


        void* setupIsoOrderedGenerator(void* iso,
                                       int _tabSize,
                                       int _hashSize);
        double massIsoOrderedGenerator(void* generator); double lprobIsoOrderedGenerator(void* generator); double probIsoOrderedGenerator(void* generator); void methodIsoOrderedGenerator(void* generator); bool advanceToNextConfigurationIsoOrderedGenerator(void* generator); void deleteIsoOrderedGenerator(void* generator); void get_conf_signatureIsoOrderedGenerator(void* generator, int* space);

        void* setupIsoStochasticGenerator(void* iso,
                                   size_t no_molecules,
                                   double precision,
                                   double beta_bias);
        double massIsoStochasticGenerator(void* generator); double lprobIsoStochasticGenerator(void* generator); double probIsoStochasticGenerator(void* generator); void methodIsoStochasticGenerator(void* generator); bool advanceToNextConfigurationIsoStochasticGenerator(void* generator); void deleteIsoStochasticGenerator(void* generator); void get_conf_signatureIsoStochasticGenerator(void* generator, int* space);

        void* setupThresholdFixedEnvelope(void* iso,
                                    double threshold,
                                    bool absolute,
                                    bool get_confs);

        void* setupTotalProbFixedEnvelope(void* iso,
                                      double taget_coverage,
                                      bool optimize,
                                      bool get_confs);

        void* setupStochasticFixedEnvelope(void* iso,
                                     size_t no_molecules,
                                     double precision,
                                     double beta_bias,
                                     bool get_confs);

        void* setupBinnedFixedEnvelope(void* iso,
                            double target_total_prob,
                            double bin_width,
                            double bin_middle);

        void* setupFixedEnvelope(double* masses, double* probs, size_t size, bool mass_sorted, bool prob_sorted, double total_prob);
        void* copyFixedEnvelope(void* other);
        void deleteFixedEnvelope(void* tabulator, bool releaseEverything);

        const double* massesFixedEnvelope(void* tabulator);
        const double* probsFixedEnvelope(void* tabulator);
        const int*    confsFixedEnvelope(void* tabulator);
        size_t confs_noFixedEnvelope(void* tabulator);

        double empiricAverageMass(void* tabulator);
        double empiricVariance(void* tabulator);
        double empiricStddev(void* tabulator);

        double wassersteinDistance(void* tabulator1, void* tabulator2);
        double orientedWassersteinDistance(void* tabulator1, void* tabulator2);
        double abyssalWassersteinDistance(void* tabulator1, void* tabulator2, double abyss_depth, double other_scale);
//        double abyssalWassersteinDistanceGrad(void* const* envelopes, const double* scales, double* ret_gradient, size_t N, double abyss_depth_exp, double abyss_depth_the);

        struct ws_match_res{
        double res1;
        double res2;
        double flow;
        };

        struct ws_match_res wassersteinMatch(void* tabulator1, void* tabulator2, double flow_dist, double other_scale);

        void* addEnvelopes(void* tabulator1, void* tabulator2);
        void* convolveEnvelopes(void* tabulator1, void* tabulator2);

        double getTotalProbOfEnvelope(void* envelope);
        void scaleEnvelope(void* envelope, double factor);
        void normalizeEnvelope(void* envelope);
        void shiftMassEnvelope(void* envelope, double d_mass);
        void resampleEnvelope(void* envelope, size_t ionic_current, double beta_bias);
        void* binnedEnvelope(void* envelope, double width, double middle);
        void* linearCombination(void* const * const envelopes, const double* intensities, size_t count);

        void sortEnvelopeByMass(void* envelope);
        void sortEnvelopeByProb(void* envelope);

        void freeReleasedArray(void* array);
        void array_add(double* array, size_t N, double what);
        void array_mul(double* array, size_t N, double what);
        void array_fma(double* array, size_t N, double mul, double add);

        void parse_fasta_c(const char* fasta, int atomCounts[6]);

        #define NUMBER_OF_ISOTOPIC_ENTRIES 292
        extern const size_t isospec_number_of_isotopic_entries;
        extern const int elem_table_atomicNo[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const double elem_table_probability[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const double elem_table_mass[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const double elem_table_massNo[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const int elem_table_extraNeutrons[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const char* elem_table_element[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const char* elem_table_symbol[NUMBER_OF_ISOTOPIC_ENTRIES];
        extern const bool elem_table_Radioactive[NUMBER_OF_ISOTOPIC_ENTRIES];
                        ''');

        mod_dir = os.path.dirname(os.path.abspath(__file__))

        if os.path.exists(os.path.join(mod_dir, '..', 'setup.py')):
            raise ImportError('''Attempted to load IsoSpecPy module from its build directory. This usually
won't work and is generally a Bad Idea. Please cd somewhere else, or, if you're really
sure you want to do that, edit the source and disable this check.''')

        libnames  = ['IsoSpecCppPy*', 'IsoSpec++*']
        libprefix = ['', 'lib', 'Lib']
        extension = ['.so', '.dylib', '.dll']
        try:
            if platform.system() == 'Linux':
                extension = ['.so', 'pyd']
            elif platform.system() == 'Windows':
                extension = ['.dll', '.pyd']
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

        try:
            import importlib
            paths_to_check.insert(0, importlib.util.find_spec("IsoSpecCppPy").origin)
        except (ImportError, AttributeError):
            pass

        errors = []

        self.clib = None
        for libpath in set(paths_to_check):
            try:
                self.clib = self.ffi.dlopen(libpath)
                break
            except (IndexError, OSError) as e:
                errmsg = "Load libIsoSpec++.so, tried: " + libpath + '\n' + "Got error: " + str(type(e)) + ": " + str(e)
                errors.append(errmsg)

        if self.clib == None:
            raise ImportError("Cannot find or load the C++ part of the library\n" + '\n'.join(errors))


isoFFI = IsoFFI()  # This is done while including the module
