/*
 *   Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/ndarray.h>

#include "../IsoSpec++/isoSpec++.h"
#include "../IsoSpec++/fixedEnvelopes.h"
#include "../IsoSpec++/fasta.h"
#include "../IsoSpec++/element_tables.h"

namespace nb = nanobind;
using namespace nb::literals;
using namespace IsoSpec;

// Helper function to convert Python buffer to int array
std::vector<int> buffer_to_int_vec(nb::ndarray<int, nb::ndim<1>, nb::c_contig> arr) {
    size_t size = arr.shape(0);
    std::vector<int> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = arr(i);
    }
    return result;
}

// Helper function to convert Python buffer to double array
std::vector<double> buffer_to_double_vec(nb::ndarray<double, nb::ndim<1>, nb::c_contig> arr) {
    size_t size = arr.shape(0);
    std::vector<double> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = arr(i);
    }
    return result;
}

NB_MODULE(_isospec_nb, m) {
    m.doc() = "IsoSpec C++ Python bindings using nanobind";

    // Iso class
    nb::class_<Iso>(m, "Iso")
        .def("__init__", [](Iso* self, int dimNumber,
                        nb::ndarray<int, nb::ndim<1>, nb::c_contig> isotopeNumbers,
                        nb::ndarray<int, nb::ndim<1>, nb::c_contig> atomCounts,
                        nb::ndarray<double, nb::ndim<1>, nb::c_contig> isotopeMasses,
                        nb::ndarray<double, nb::ndim<1>, nb::c_contig> isotopeProbabilities) {
            auto isoNums = buffer_to_int_vec(isotopeNumbers);
            auto atomCts = buffer_to_int_vec(atomCounts);
            auto isoMasses = buffer_to_double_vec(isotopeMasses);
            auto isoProbs = buffer_to_double_vec(isotopeProbabilities);
            new (self) Iso(dimNumber, isoNums.data(), atomCts.data(),
                          isoMasses.data(), isoProbs.data());
        }, "dimNumber"_a, "isotopeNumbers"_a, "atomCounts"_a,
            "isotopeMasses"_a, "isotopeProbabilities"_a)
        .def_static("from_fasta",
            [](const std::string& fasta, bool use_nominal_masses, bool add_water) {
                return new Iso(Iso::FromFASTA(fasta, use_nominal_masses, add_water));
            }, "fasta"_a, "use_nominal_masses"_a = false, "add_water"_a = true)
        .def("get_lightest_peak_mass", &Iso::getLightestPeakMass)
        .def("get_lightest_peak_lprob", &Iso::getLightestPeakLProb)
        .def("get_lightest_peak_signature", [](Iso& self) {
            std::vector<int> result(self.getAllDim());
            self.getLightestPeakSignature(result.data());
            return result;
        })
        .def("get_heaviest_peak_mass", &Iso::getHeaviestPeakMass)
        .def("get_heaviest_peak_lprob", &Iso::getHeaviestPeakLProb)
        .def("get_heaviest_peak_signature", [](Iso& self) {
            std::vector<int> result(self.getAllDim());
            self.getHeaviestPeakSignature(result.data());
            return result;
        })
        .def("get_monoisotopic_peak_mass", &Iso::getMonoisotopicPeakMass)
        .def("get_monoisotopic_peak_lprob", &Iso::getMonoisotopicPeakLProb)
        .def("get_monoisotopic_peak_signature", [](Iso& self) {
            std::vector<int> result(self.getAllDim());
            self.getMonoisotopicPeakSignature(result.data());
            return result;
        })
        .def("get_mode_lprob", &Iso::getModeLProb)
        .def("get_mode_mass", &Iso::getModeMass)
        .def("get_theoretical_average_mass", &Iso::getTheoreticalAverageMass)
        .def("variance", &Iso::variance)
        .def("stddev", &Iso::stddev)
        .def("get_dim_number", &Iso::getDimNumber)
        .def("get_all_dim", &Iso::getAllDim)
        .def("get_marginal_log_size_estimates", [](Iso& self, double target_total_prob) {
            std::vector<double> result(self.getDimNumber());
            self.saveMarginalLogSizeEstimates(result.data(), target_total_prob);
            return result;
        }, "target_total_prob"_a);

    // IsoGenerator base class
    nb::class_<IsoGenerator, Iso>(m, "IsoGenerator")
        .def("advance_to_next_configuration", &IsoGenerator::advanceToNextConfiguration)
        .def("lprob", &IsoGenerator::lprob)
        .def("mass", &IsoGenerator::mass)
        .def("prob", &IsoGenerator::prob)
        .def("get_conf_signature", [](IsoGenerator& self) {
            std::vector<int> result(self.getAllDim());
            self.get_conf_signature(result.data());
            return result;
        })
        .def_prop_ro("mode_lprob", [](const IsoGenerator& self) { return self.mode_lprob; });

    // IsoThresholdGenerator
    nb::class_<IsoThresholdGenerator, IsoGenerator>(m, "IsoThresholdGenerator")
        .def("__init__", [](IsoThresholdGenerator* self, Iso* iso, double threshold, bool absolute,
                        int tabSize, int hashSize, bool reorder_marginals) {
            // Create a copy of iso for the generator to consume via move
            Iso iso_copy(*iso, true);
            new (self) IsoThresholdGenerator(std::move(iso_copy), threshold, absolute,
                                            tabSize, hashSize, reorder_marginals);
        }, "iso"_a, "threshold"_a, "absolute"_a = true,
            "tabSize"_a = 1000, "hashSize"_a = 1000, "reorder_marginals"_a = true);

    // IsoLayeredGenerator
    nb::class_<IsoLayeredGenerator, IsoGenerator>(m, "IsoLayeredGenerator")
        .def("__init__", [](IsoLayeredGenerator* self, Iso* iso, int tabSize, int hashSize,
                        bool reorder_marginals, double t_prob_hint) {
            Iso iso_copy(*iso, true);
            new (self) IsoLayeredGenerator(std::move(iso_copy), tabSize, hashSize,
                                          reorder_marginals, t_prob_hint);
        }, "iso"_a, "tabSize"_a = 1000, "hashSize"_a = 1000,
            "reorder_marginals"_a = true, "t_prob_hint"_a = 0.01);

    // IsoOrderedGenerator
    nb::class_<IsoOrderedGenerator, IsoGenerator>(m, "IsoOrderedGenerator")
        .def("__init__", [](IsoOrderedGenerator* self, Iso* iso, int tabSize, int hashSize) {
            Iso iso_copy(*iso, true);
            new (self) IsoOrderedGenerator(std::move(iso_copy), tabSize, hashSize);
        }, "iso"_a, "tabSize"_a = 1000, "hashSize"_a = 1000);

    // IsoStochasticGenerator
    nb::class_<IsoStochasticGenerator, IsoGenerator>(m, "IsoStochasticGenerator")
        .def("__init__", [](IsoStochasticGenerator* self, Iso* iso, size_t no_molecules,
                        double precision, double beta_bias) {
            Iso iso_copy(*iso, true);
            new (self) IsoStochasticGenerator(std::move(iso_copy), no_molecules,
                                             precision, beta_bias);
        }, "iso"_a, "no_molecules"_a, "precision"_a = 0.9999, "beta_bias"_a = 5.0);

    // FixedEnvelope class
    nb::class_<FixedEnvelope>(m, "FixedEnvelope")
        .def("__init__", [](FixedEnvelope* self, nb::ndarray<double, nb::ndim<1>, nb::c_contig> masses,
                        nb::ndarray<double, nb::ndim<1>, nb::c_contig> probs,
                        bool mass_sorted, bool prob_sorted, double total_prob) {
            size_t size = masses.shape(0);
            auto masses_vec = buffer_to_double_vec(masses);
            auto probs_vec = buffer_to_double_vec(probs);
            double* m_ptr = (double*)malloc(size * sizeof(double));
            double* p_ptr = (double*)malloc(size * sizeof(double));
            std::copy(masses_vec.begin(), masses_vec.end(), m_ptr);
            std::copy(probs_vec.begin(), probs_vec.end(), p_ptr);
            new (self) FixedEnvelope(m_ptr, p_ptr, size, mass_sorted, prob_sorted, total_prob);
        }, "masses"_a, "probs"_a, "mass_sorted"_a = false,
            "prob_sorted"_a = false, "total_prob"_a = NAN)
        .def("__init__", [](FixedEnvelope* self, nb::ndarray<double, nb::ndim<1>, nb::c_contig> masses,
                        nb::ndarray<double, nb::ndim<1>, nb::c_contig> probs,
                        nb::ndarray<int, nb::ndim<1>, nb::c_contig> confs,
                        int allDim, bool mass_sorted, bool prob_sorted, double total_prob) {
            size_t size = masses.shape(0);
            auto masses_vec = buffer_to_double_vec(masses);
            auto probs_vec = buffer_to_double_vec(probs);
            double* m_ptr = (double*)malloc(size * sizeof(double));
            double* p_ptr = (double*)malloc(size * sizeof(double));
            int* c_ptr = (int*)malloc(size * allDim * sizeof(int));
            std::copy(masses_vec.begin(), masses_vec.end(), m_ptr);
            std::copy(probs_vec.begin(), probs_vec.end(), p_ptr);
            // Copy 2D confs array (flattened)
            for (size_t i = 0; i < size; i++) {
                for (int j = 0; j < allDim; j++) {
                    c_ptr[i * allDim + j] = confs(i * allDim + j);
                }
            }
            new (self) FixedEnvelope(m_ptr, p_ptr, c_ptr, size, allDim,
                                    mass_sorted, prob_sorted, total_prob);
        }, "masses"_a, "probs"_a, "confs"_a, "allDim"_a,
            "mass_sorted"_a = false, "prob_sorted"_a = false, "total_prob"_a = NAN)
        .def(nb::init<const FixedEnvelope&>())
        .def_static("from_threshold",
            [](Iso* iso, double threshold, bool absolute, bool get_confs) {
                return new FixedEnvelope(
                    FixedEnvelope::FromThreshold(
                        Iso(*iso, true), threshold, absolute, get_confs));
            }, "iso"_a, "threshold"_a, "absolute"_a = true, "get_confs"_a = false)
        .def_static("from_total_prob",
            [](Iso* iso, double target_coverage, bool optimize, bool get_confs) {
                return new FixedEnvelope(
                    FixedEnvelope::FromTotalProb(
                        Iso(*iso, true), target_coverage, optimize, get_confs));
            }, "iso"_a, "target_coverage"_a, "optimize"_a = true, "get_confs"_a = false)
        .def_static("from_stochastic",
            [](Iso* iso, size_t no_molecules, double precision,
               double beta_bias, bool get_confs) {
                return new FixedEnvelope(
                    FixedEnvelope::FromStochastic(
                        Iso(*iso, true), no_molecules, precision, beta_bias, get_confs));
            }, "iso"_a, "no_molecules"_a, "precision"_a = 0.9999,
               "beta_bias"_a = 5.0, "get_confs"_a = false)
        .def_static("binned",
            [](Iso* iso, double target_total_prob, double bin_width, double bin_middle) {
                return new FixedEnvelope(
                    FixedEnvelope::Binned(Iso(*iso, true), target_total_prob,
                                         bin_width, bin_middle));
            }, "iso"_a, "target_total_prob"_a, "bin_width"_a, "bin_middle"_a)
        .def("confs_no", &FixedEnvelope::confs_no)
        .def("get_all_dim", &FixedEnvelope::getAllDim)
        .def("masses", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            const double* masses = self.masses();
            std::vector<double> result(masses, masses + size);
            return result;
        })
        .def("probs", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            const double* probs = self.probs();
            std::vector<double> result(probs, probs + size);
            return result;
        })
        .def("confs", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            int allDim = self.getAllDim();
            const int* confs = self.confs();
            if (confs == nullptr) {
                return std::vector<std::vector<int>>();
            }
            std::vector<std::vector<int>> result;
            for (size_t i = 0; i < size; i++) {
                std::vector<int> conf(allDim);
                for (int j = 0; j < allDim; j++) {
                    conf[j] = confs[i * allDim + j];
                }
                result.push_back(conf);
            }
            return result;
        })
        .def("release_masses", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            double* masses = self.release_masses();
            if (masses == nullptr) {
                return std::vector<double>();
            }
            std::vector<double> result(masses, masses + size);
            free(masses);
            return result;
        })
        .def("release_probs", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            double* probs = self.release_probs();
            if (probs == nullptr) {
                return std::vector<double>();
            }
            std::vector<double> result(probs, probs + size);
            free(probs);
            return result;
        })
        .def("release_confs", [](FixedEnvelope& self) {
            size_t size = self.confs_no();
            int allDim = self.getAllDim();
            int* confs = self.release_confs();
            if (confs == nullptr) {
                return std::vector<std::vector<int>>();
            }
            std::vector<std::vector<int>> result;
            for (size_t i = 0; i < size; i++) {
                std::vector<int> conf(allDim);
                for (int j = 0; j < allDim; j++) {
                    conf[j] = confs[i * allDim + j];
                }
                result.push_back(conf);
            }
            free(confs);
            return result;
        })
        .def("sort_by_mass", &FixedEnvelope::sort_by_mass)
        .def("sort_by_prob", &FixedEnvelope::sort_by_prob)
        .def("get_total_prob", &FixedEnvelope::get_total_prob)
        .def("scale", &FixedEnvelope::scale, "factor"_a)
        .def("normalize", &FixedEnvelope::normalize)
        .def("shift_mass", &FixedEnvelope::shift_mass, "shift"_a)
        .def("resample", &FixedEnvelope::resample,
             "ionic_current"_a, "beta_bias"_a = 1.0)
        .def("empiric_average_mass", &FixedEnvelope::empiric_average_mass)
        .def("empiric_variance", &FixedEnvelope::empiric_variance)
        .def("empiric_stddev", &FixedEnvelope::empiric_stddev)
        .def("wasserstein_distance", &FixedEnvelope::WassersteinDistance, "other"_a)
        .def("oriented_wasserstein_distance",
             &FixedEnvelope::OrientedWassersteinDistance, "other"_a)
        .def("abyssal_wasserstein_distance",
             &FixedEnvelope::AbyssalWassersteinDistance,
             "other"_a, "abyss_depth"_a, "other_scale"_a = 1.0)
        .def("wasserstein_match", &FixedEnvelope::WassersteinMatch,
             "other"_a, "flow_distance"_a, "other_scale"_a = 1.0)
        .def("bin", &FixedEnvelope::bin, "bin_width"_a = 1.0, "middle"_a = 0.0)
        .def("__add__", &FixedEnvelope::operator+)
        .def("__mul__", &FixedEnvelope::operator*)
        .def_static("linear_combination",
            [](const std::vector<FixedEnvelope*>& envelopes,
               const std::vector<double>& intensities) {
                std::vector<const FixedEnvelope*> const_envs;
                const_envs.reserve(envelopes.size());
                for (auto* env : envelopes) {
                    const_envs.push_back(env);
                }
                return new FixedEnvelope(
                    FixedEnvelope::LinearCombination(const_envs, intensities));
            }, "envelopes"_a, "intensities"_a);

    // FASTA parsing helper
    m.def("parse_fasta", [](const std::string& fasta) {
        int atomCounts[6] = {0};
        IsoSpec::parse_fasta(fasta.c_str(), atomCounts);
        return std::vector<int>(atomCounts, atomCounts + 6);
    }, "fasta"_a);

    // Periodic table data accessors
    m.attr("NUMBER_OF_ISOTOPIC_ENTRIES") = (int)ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES;

    m.def("elem_table_symbol", [](int idx) -> std::string {
        if (idx < 0 || idx >= ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES)
            throw std::out_of_range("Index out of range");
        return std::string(elem_table_symbol[idx]);
    }, "idx"_a);

    m.def("elem_table_mass", [](int idx) -> double {
        if (idx < 0 || idx >= ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES)
            throw std::out_of_range("Index out of range");
        return elem_table_mass[idx];
    }, "idx"_a);

    m.def("elem_table_massNo", [](int idx) -> double {
        if (idx < 0 || idx >= ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES)
            throw std::out_of_range("Index out of range");
        return elem_table_massNo[idx];
    }, "idx"_a);

    m.def("elem_table_probability", [](int idx) -> double {
        if (idx < 0 || idx >= ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES)
            throw std::out_of_range("Index out of range");
        return elem_table_probability[idx];
    }, "idx"_a);

    m.def("elem_table_atomicNo", [](int idx) -> int {
        if (idx < 0 || idx >= ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES)
            throw std::out_of_range("Index out of range");
        return elem_table_atomicNo[idx];
    }, "idx"_a);
}
