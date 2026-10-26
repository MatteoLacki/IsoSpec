/*
 *   Copyright (C) 2015-2018 Mateusz Łącki and Michał Startek.
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


#include <Rcpp.h>
#include "cwrapper.h"
#include "misc.h"
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "fasta_mods.h"
#include "unimod.h"
#include "element_tables.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <map>
#include <cmath>
#include <cstdint>
#include <random>

using namespace Rcpp;
using namespace IsoSpec;

FixedEnvelope mkIsoG(Iso& iso, int algo, double stopCondition, bool trim, bool get_confs)
{
    switch(algo)
    {
        case ISOSPEC_ALGO_LAYERED_ESTIMATE: // Not used anymore, just fall through to the next case
        case ISOSPEC_ALGO_LAYERED: return FixedEnvelope::FromTotalProb(std::move(iso), stopCondition, trim, get_confs);
        case ISOSPEC_ALGO_ORDERED: return FixedEnvelope::FromTotalProb(std::move(iso), stopCondition, true, get_confs); // Using the ordered algo in R is *completely* pointless today
                                                                                                                  // The only point of ordered algo is when one is using the generator
                                                                                                                  // interface, which we are not exposing in R
                                                                                                                  // We'll just do layered, trim and sort it afterwards - it's equivalent
                                                                                                                  // and much faster
        case ISOSPEC_ALGO_THRESHOLD_ABSOLUTE:
        case ISOSPEC_ALGO_THRESHOLD_RELATIVE: throw std::logic_error("");
    }
    throw std::logic_error("Invalid algorithm selected");
}

// The molecule as IsoSpecify() and the other molecule-level functions take
// it: a named vector of atom counts, each element's isotopes looked up by
// symbol in the `isotopes` data frame (masses divided by `charge`). With
// `isotopeTags`, the labels of the isotopes are appended to it, in the order
// of a configuration's counts.
static Iso iso_from_R(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  charge,
    CharacterVector*        isotopeTags = nullptr
){
    unsigned int dimNumber = molecule.size();
    std::vector<int>     stdIsotopeNumbers;
    std::vector<double> stdIsotopeMasses;
    std::vector<double> stdIsotopeProbabilities;

    const CharacterVector& element = isotopes["element"];
    const CharacterVector& isotope = isotopes["isotope"];
    const NumericVector& mass      = isotopes["mass"];
    const NumericVector& abundance = isotopes["abundance"];
    const CharacterVector& molecule_names = molecule.names();

    for (unsigned int i=0; i<dimNumber; i++)
    {
        // The sibling of the counter == 0 check below, and reachable the same
        // way: RParsePeptideSequence returns a *signed* composition (a
        // modification's delta can remove more atoms of an element than the
        // bare sequence provides -- "G[UNIMOD:11]" nets H-1 S-1), and the
        // documented two-call pattern feeds its result straight back in here.
        // A negative atom count reaching Iso is undefined behaviour in
        // exactly the same way a zero-isotope dimension is; observed on this
        // build as a silently wrong answer (the negative elements' mass just
        // missing from every peak), and as a heap-buffer-overflow under ASan
        // on the equivalent C++ path.
        if( molecule[i] < 0 )
            throw std::invalid_argument(
                "Negative atom count for element '" + Rcpp::as<std::string>(molecule_names[i]) +
                "' -- a molecule can't contain a negative number of atoms. If this vector came "
                "from RParsePeptideSequence, the sequence's modifications remove more of this "
                "element than the bare sequence provides.");
        unsigned int counter = 0;
        for (int j=0; j<element.size(); j++)
            if( element[j] == molecule_names[i] )
            {
                counter++;
                stdIsotopeMasses.push_back( mass[j] / charge );
                stdIsotopeProbabilities.push_back( abundance[j] );
                if( isotopeTags != nullptr )
                    isotopeTags->push_back( isotope[j] );
            }
        // A molecule element with zero matching rows in `isotopes` previously
        // went through silently as a zero-isotope dimension -- Iso's
        // constructor then had no known mass/probability data for it at all,
        // so that atom's contribution was just dropped from the output
        // (wrong mass, no error) rather than erroring, and in at least one
        // reported case (github.com/MatteoLacki/IsoSpec/issues/49) corrupted
        // memory badly enough to crash the R session outright on repeated
        // calls. Fail loudly and name the symbol instead: either it's a typo,
        // or (e.g. an element with no stable natural isotopes, like Tc or Ac)
        // it genuinely isn't in this table and the caller needs to supply its
        // own isotope data via the `isotopes=` parameter -- see the
        // radiolabelling example for that pattern. Rcpp turns an escaped C++
        // exception into an R error automatically, no manual try/catch here.
        if( counter == 0 )
            throw std::invalid_argument(
                "No isotope data for element '" + Rcpp::as<std::string>(molecule_names[i]) +
                "' in the supplied `isotopes` table -- check for a typo, or supply this "
                "element's isotope data yourself via the `isotopes` parameter (see the "
                "radiolabelling example) if it genuinely has none in the default table "
                "(e.g. an element with no stable natural isotopes).");
        stdIsotopeNumbers.push_back(counter);
    }

    std::vector<double*> IM;
    std::vector<double*> IP;
    size_t tot = 0;
    for(size_t ii = 0; ii<stdIsotopeNumbers.size(); ii++)
    {
        IM.push_back(&(stdIsotopeMasses.data()[tot]));
        IP.push_back(&(stdIsotopeProbabilities.data()[tot]));
        tot += stdIsotopeNumbers[ii];
    }
    return Iso(dimNumber, stdIsotopeNumbers.data(), Rcpp::as<std::vector<int> >( molecule).data(), IM.data(), IP.data());

}

// [[Rcpp::export]]
NumericMatrix Rinterface(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  stopCondition,
    int                     algo = 0,
    int                     tabSize = 64,
    int                     hashSize = 1000,
    double                  step = .3,
    bool                    showCounts = false,
    bool                    trim = true,
    double                  charge = 1.0
){

    CharacterVector stdIsotopeTags = CharacterVector::create("mass", "prob");
    Iso iso = iso_from_R(molecule, isotopes, charge, showCounts ? &stdIsotopeTags : nullptr);

    unsigned int columnsNo = stdIsotopeTags.size(); // standard

    unsigned int isotopesNo = iso.getAllDim();

    if(algo == ISOSPEC_ALGO_THRESHOLD_ABSOLUTE || algo == ISOSPEC_ALGO_THRESHOLD_RELATIVE)
    {
        IsoThresholdGenerator ITG(std::move(iso), stopCondition, (algo == ISOSPEC_ALGO_THRESHOLD_ABSOLUTE));

        size_t no_confs = ITG.count_confs();
        size_t ii = 0;

        NumericMatrix res(no_confs, columnsNo);

        if(showCounts)
        {
            int* conf_sig = new int[isotopesNo];
            while(ITG.advanceToNextConfiguration())
            {
                res(ii,0) = ITG.mass();
                res(ii,1) = ITG.prob();
                ITG.get_conf_signature(conf_sig);
                for(size_t jj = 0; jj < isotopesNo; jj++)
                    res(ii, 2+jj) = conf_sig[jj];
                ii++;
            }
            delete[] conf_sig;
        }
        else
            while(ITG.advanceToNextConfiguration())
            {
                res(ii,0) = ITG.mass();
                res(ii,1) = ITG.prob();
                ii++;
            }

	colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.
        return res;
    }

    // The remaining (layered) algos
    FixedEnvelope TAB = mkIsoG(iso, algo, stopCondition, trim, showCounts);

    const double* probs = TAB.probs();
    const double* masses = TAB.masses();
    const int* confs = TAB.confs();
    std::vector<size_t> ordering;
    size_t confs_no = TAB.confs_no();


    const bool needs_sorting = (ISOSPEC_ALGO_ORDERED == algo);

    const unsigned int isotopesNoplus2 = isotopesNo + 2;

    NumericMatrix res(confs_no, columnsNo);

    if(needs_sorting)
    {
        // We need to sort the confs for backward compatibility
        ordering.reserve(confs_no);
        for(size_t i = 0; i < confs_no; i++)
            ordering.push_back(i);
        std::sort(ordering.begin(), ordering.end(), [&probs](size_t idx1, size_t idx2) -> bool { return probs[idx1] > probs[idx2]; });

        unsigned int idx;
        for(size_t i = 0; i < confs_no; i++)
        {
            idx = ordering[i];
            res(i,0) = masses[idx];
            res(i,1) = probs[idx];
        }

        if(showCounts)
        {
            size_t confs_idx;
            for(size_t i = 0; i < confs_no; i++)
            {
                // Configurations are flat rows of isotopesNo counts.
                confs_idx = ordering[i] * static_cast<size_t>(isotopesNo);
                for(size_t j = 2; j < isotopesNoplus2; j++)
                {
                    res(i,j) = confs[confs_idx];
                    confs_idx++;
                }
            }
        }
    }
    else
    {
        for(size_t i = 0; i < confs_no; i++)
        {
            res(i,0) = masses[i];
            res(i,1) = probs[i];
        }

        if(showCounts)
        {
            size_t confs_idx = 0;
            for(size_t i = 0; i < confs_no; i++)
            {
                for(size_t j = 2; j < isotopesNoplus2; j++)
                {
                    res(i,j) = confs[confs_idx];
                    confs_idx++;
                }
            }
        }
    }

    colnames(res) = stdIsotopeTags; //This is RCPP sugar. It sucks.

    return(res);
}

// Unimod-modification-aware peptide sequence parsing -- see
// docs/ai/unimod.md. Unlike Rinterface above (which takes a raw atom-count
// vector), this binds straight to the C++ classes the same way Rinterface
// itself does, not through cwrapper.h -- R had no FASTA/sequence entry point
// at all before this, so there is nothing here to stay compatible with.
//
// [[Rcpp::export]]
IntegerVector RParsePeptideSequence(
    const std::string&      sequence,
    std::string             unimod_db_path = ""
){
    const UnimodTable& mods = unimod_db_path.empty() ? embedded_unimod_table()
                                                       : unimod_table_for_path(unimod_db_path);
    // Rcpp translates an escaped C++ exception (std::invalid_argument on a
    // malformed bracket or an unknown/excluded id) into an R error
    // automatically -- no manual try/catch needed here.
    ElementComposition composition = parse_fasta_with_mods(sequence.c_str(), mods);

    IntegerVector res(composition.count.size());
    CharacterVector symbols(composition.count.size());
    for (size_t i = 0; i < composition.count.size(); i++)
    {
        res[i] = composition.count[i];
        symbols[i] = elem_table_symbol[composition.element_first_index[i]];
    }
    res.names() = symbols;

    return(res);
}

// Reads the core element/isotope table directly out of element_tables.cpp --
// the single authoritative source every other binding (Python's PeriodicTbl,
// built the same way via cffi) already resolves against -- instead of R
// carrying its own hand-maintained copy. That duplication is exactly what
// caused the #51 (deuterium) fix to require four separate edits, two of
// which (data/isotopicData.rda and R/sysdata.rda) had already silently
// drifted from each other before this function existed. See IsoSpecR.R's
// IsoSpecify(), which now builds its default `isotopes` table from this.
//
// Excludes the synthetic charge-state pseudo-elements E/Me/Pn (elem_table_ID
// 1000-1002) by symbol, not by "ID >= 1000", because D (deuterium, ID 1003)
// is exactly as synthetic by that measure and must NOT be excluded -- it's
// the whole point. isotopicData$IsoSpec has always had 0 rows for E/Me/Pn
// (verified against data/isotopicData.rda at commit ada5ae5); this preserves
// that.
//
// Isotope labels follow this table's own existing convention (element
// symbol + mass number, e.g. "H1", "H2", "C12") -- but for a synthetic
// pseudo-element like D, the label uses its *parent* element's symbol (the
// real, non-synthetic entry sharing its atomic number: D's atomicNo is 1,
// hydrogen's), not its own ("H2", not "D2"). This matches the existing "D"
// row already present in both data/isotopicData.rda and the external
// enviPat package's own data (enviPat's D row predates this session and was
// left untouched) -- confirmed byte-for-byte against both this session.
//
// [[Rcpp::export]]
DataFrame RIsotopicTable()
{
    static const std::string kExcludedSymbols[] = {"E", "Me", "Pn"};
    auto is_excluded = [](const std::string& sym) {
        for (const std::string& excl : kExcludedSymbols)
            if (sym == excl)
                return true;
        return false;
    };

    // atomicNo -> the real (non-synthetic) element symbol sharing it, used
    // to label a synthetic pseudo-element's isotope by its parent's symbol.
    std::map<int, std::string> parent_symbol_by_atomic_no;
    for (int i = 0; i < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; i++)
    {
        std::string sym = elem_table_symbol[i];
        if (is_excluded(sym) || elem_table_ID[i] >= 1000)
            continue;
        parent_symbol_by_atomic_no.emplace(elem_table_atomicNo[i], sym);
    }

    std::vector<std::string> element, isotope;
    std::vector<double> mass, abundance;

    for (int i = 0; i < ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES; i++)
    {
        std::string sym = elem_table_symbol[i];
        if (is_excluded(sym))
            continue;

        std::string label_symbol = sym;
        if (elem_table_ID[i] >= 1000)
        {
            auto it = parent_symbol_by_atomic_no.find(elem_table_atomicNo[i]);
            if (it != parent_symbol_by_atomic_no.end())
                label_symbol = it->second;
        }

        element.push_back(sym);
        isotope.push_back(label_symbol + std::to_string(static_cast<int>(elem_table_massNo[i])));
        mass.push_back(elem_table_mass[i]);
        abundance.push_back(elem_table_probability[i]);
    }

    return DataFrame::create(
        Named("element") = element,
        Named("isotope") = isotope,
        Named("mass") = mass,
        Named("abundance") = abundance,
        Named("stringsAsFactors") = false
    );
}

// ---------------------------------------------------------------------------
// Binned and stochastic distributions, molecule properties, and the envelope
// operations. An envelope crosses the boundary as two plain vectors (masses,
// probabilities) and comes back as a two-column matrix, "mass" and "prob",
// the shape IsoSpecify() returns: R code holds nothing of the C++ side between
// calls. Operations that are one line of R on such a matrix (shifting masses,
// scaling or normalising probabilities, sorting, weighted means) are not
// duplicated here.

static NumericMatrix envelope_to_R(const FixedEnvelope& env)
{
    const size_t n = env.confs_no();
    NumericMatrix res(n, 2);
    for(size_t i = 0; i < n; i++)
    {
        res(i, 0) = env.mass(i);
        res(i, 1) = env.prob(i);
    }
    colnames(res) = CharacterVector::create("mass", "prob");
    return res;
}

static FixedEnvelope envelope_from_R(const NumericVector& masses, const NumericVector& probs)
{
    if(masses.size() != probs.size())
        throw std::invalid_argument("An envelope needs as many masses as probabilities.");
    // Copies, and checks what every operation relies on (finite masses,
    // finite nonnegative probabilities).
    return FixedEnvelope::FromArrays(masses.begin(), probs.begin(), masses.size());
}

// NULL: the calling thread's randomly seeded generator. A number: a seed,
// giving the same stream as the C++, Python and Rust interfaces for that seed.
static std::mt19937 rng_from_R(const Nullable<NumericVector>& seed, bool& seeded)
{
    seeded = seed.isNotNull();
    if(!seeded)
        return std::mt19937();
    NumericVector s(seed);
    // 2^53: beyond it a double no longer holds every integer exactly.
    if(s.size() != 1 || !std::isfinite(s[0]) || s[0] < 0 || s[0] != std::floor(s[0]) || s[0] > 9007199254740992.0)
        throw std::invalid_argument("seed must be NULL or a single whole number between 0 and 2^53.");
    return seeded_rng(static_cast<std::uint64_t>(s[0]));
}

static size_t count_from_R(double count, const char* name)
{
    // Up to 2^53 so that the double holds it exactly.
    if(!std::isfinite(count) || count < 0 || count != std::floor(count) || count > 9007199254740992.0)
        throw std::invalid_argument(std::string(name) + " must be a whole number between 0 and 2^53.");
    return static_cast<size_t>(count);
}

// [[Rcpp::export]]
NumericMatrix RBinned(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  stopCondition,
    double                  binWidth,
    double                  binMiddle,
    double                  charge
){
    Iso iso = iso_from_R(molecule, isotopes, charge);
    return envelope_to_R(FixedEnvelope::Binned(std::move(iso), stopCondition, binWidth, binMiddle));
}

// [[Rcpp::export]]
NumericMatrix RStochastic(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  ionicCurrent,
    double                  precision,
    double                  betaBias,
    Nullable<NumericVector> seed,
    bool                    showCounts,
    double                  charge
){
    CharacterVector tags = CharacterVector::create("mass", "prob");
    Iso iso = iso_from_R(molecule, isotopes, charge, showCounts ? &tags : nullptr);
    const size_t molecules = count_from_R(ionicCurrent, "ionicCurrent");
    bool seeded;
    std::mt19937 rng = rng_from_R(seed, seeded);
    FixedEnvelope env = FixedEnvelope::FromStochastic(std::move(iso), molecules, precision, betaBias,
                                                      showCounts, seeded ? rng : random_gen);
    if(!showCounts)
        return envelope_to_R(env);
    const size_t n = env.confs_no();
    const int dim = env.getAllDim();
    NumericMatrix res(n, 2 + dim);
    for(size_t i = 0; i < n; i++)
    {
        res(i, 0) = env.mass(i);
        res(i, 1) = env.prob(i);
        const int* conf = env.conf(i);
        for(int j = 0; j < dim; j++)
            res(i, 2 + j) = conf[j];
    }
    colnames(res) = tags;
    return res;
}

// [[Rcpp::export]]
List RIsoProperties(
    const IntegerVector&    molecule,
    const DataFrame&        isotopes,
    double                  charge
){
    CharacterVector tags;
    Iso iso = iso_from_R(molecule, isotopes, charge, &tags);
    const int dim = iso.getAllDim();
    auto conf = [&](void (Iso::*getter)(int*) const) {
        IntegerVector counts(dim);
        (iso.*getter)(counts.begin());
        counts.names() = tags;
        return counts;
    };
    return List::create(
        Named("monoisotopicMass")   = iso.getMonoisotopicPeakMass(),
        Named("monoisotopicLogProb") = iso.getMonoisotopicPeakLProb(),
        Named("monoisotopicCounts") = conf(&Iso::getMonoisotopicPeakSignature),
        Named("lightestMass")       = iso.getLightestPeakMass(),
        Named("lightestLogProb")    = iso.getLightestPeakLProb(),
        Named("lightestCounts")     = conf(&Iso::getLightestPeakSignature),
        Named("heaviestMass")       = iso.getHeaviestPeakMass(),
        Named("heaviestLogProb")    = iso.getHeaviestPeakLProb(),
        Named("heaviestCounts")     = conf(&Iso::getHeaviestPeakSignature),
        Named("modeMass")           = iso.getModeMass(),
        Named("modeLogProb")        = iso.getModeLProb(),
        Named("averageMass")        = iso.getTheoreticalAverageMass(),
        Named("variance")           = iso.variance());
}

// [[Rcpp::export]]
NumericMatrix REnvelopeConvolve(const NumericVector& masses1, const NumericVector& probs1,
                                const NumericVector& masses2, const NumericVector& probs2)
{
    return envelope_to_R(envelope_from_R(masses1, probs1) * envelope_from_R(masses2, probs2));
}

// [[Rcpp::export]]
NumericMatrix REnvelopeLinearCombination(const List& masses, const List& probs, const NumericVector& intensities)
{
    if(masses.size() != probs.size() || masses.size() != intensities.size())
        throw std::invalid_argument("As many intensities as envelopes are needed.");
    std::vector<FixedEnvelope> envelopes;
    envelopes.reserve(masses.size());
    for(R_xlen_t i = 0; i < masses.size(); i++)
        envelopes.push_back(envelope_from_R(masses[i], probs[i]));
    std::vector<const FixedEnvelope*> pointers;
    for(const FixedEnvelope& env : envelopes)
        pointers.push_back(&env);
    return envelope_to_R(FixedEnvelope::LinearCombination(
        pointers, std::vector<double>(intensities.begin(), intensities.end())));
}

// [[Rcpp::export]]
NumericMatrix REnvelopeBin(const NumericVector& masses, const NumericVector& probs, double width, double middle)
{
    return envelope_to_R(envelope_from_R(masses, probs).bin(width, middle));
}

// [[Rcpp::export]]
NumericMatrix REnvelopeResample(const NumericVector& masses, const NumericVector& probs,
                                double ionicCurrent, double betaBias, Nullable<NumericVector> seed)
{
    FixedEnvelope env = envelope_from_R(masses, probs);
    const size_t molecules = count_from_R(ionicCurrent, "ionicCurrent");
    bool seeded;
    std::mt19937 rng = rng_from_R(seed, seeded);
    env.resample(molecules, betaBias, seeded ? rng : random_gen);
    return envelope_to_R(env);
}

// kind: 0 plain, 1 oriented, 2 abyssal.
// [[Rcpp::export]]
double RWasserstein(const NumericVector& masses1, const NumericVector& probs1,
                    const NumericVector& masses2, const NumericVector& probs2,
                    int kind, double abyssDepth, double otherScale)
{
    FixedEnvelope a = envelope_from_R(masses1, probs1);
    FixedEnvelope b = envelope_from_R(masses2, probs2);
    switch(kind)
    {
        case 0: return a.WassersteinDistance(b);
        case 1: return a.OrientedWassersteinDistance(b);
        case 2: return a.AbyssalWassersteinDistance(b, abyssDepth, otherScale);
    }
    throw std::logic_error("Invalid Wasserstein distance kind");
}

// [[Rcpp::export]]
NumericVector RWassersteinMatch(const NumericVector& masses1, const NumericVector& probs1,
                                const NumericVector& masses2, const NumericVector& probs2,
                                double flowDistance, double otherScale)
{
    FixedEnvelope a = envelope_from_R(masses1, probs1);
    FixedEnvelope b = envelope_from_R(masses2, probs2);
    const auto [unmatched_this, unmatched_other, flow] = a.WassersteinMatch(b, flowDistance, otherScale);
    return NumericVector::create(Named("unmatchedThis") = unmatched_this,
                                 Named("unmatchedOther") = unmatched_other,
                                 Named("flow") = flow);
}
