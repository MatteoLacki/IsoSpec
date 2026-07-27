// FASTA (aminoacid sequence) input.
//
// The oracle here is the standard table of aminoacid *residue* compositions
// from biochemistry, written out below — not another IsoSpec code path.

#include <string>
#include <vector>

#include "doctest.h"
#include "fasta.h"
#include "test_helpers.h"

using namespace IsoSpec;
using namespace test_helpers;

namespace {

// Residue compositions in the library's C,H,N,O,S,Se order.  A residue is the
// aminoacid minus the water lost in forming the peptide bond.
struct Residue {
    char code;
    int C, H, N, O, S, Se;
};

const Residue kResidues[] = {
    {'G', 2,  3, 1, 1, 0, 0},
    {'A', 3,  5, 1, 1, 0, 0},
    {'S', 3,  5, 1, 2, 0, 0},
    {'P', 5,  7, 1, 1, 0, 0},
    {'V', 5,  9, 1, 1, 0, 0},
    {'T', 4,  7, 1, 2, 0, 0},
    {'C', 3,  5, 1, 1, 1, 0},
    {'L', 6, 11, 1, 1, 0, 0},
    {'I', 6, 11, 1, 1, 0, 0},
    {'J', 6, 11, 1, 1, 0, 0},  // xleucine: L or I, same formula
    {'N', 4,  6, 2, 2, 0, 0},
    {'D', 4,  5, 1, 3, 0, 0},
    {'Q', 5,  8, 2, 2, 0, 0},
    {'K', 6, 12, 2, 1, 0, 0},
    {'E', 5,  7, 1, 3, 0, 0},
    {'M', 5,  9, 1, 1, 1, 0},
    {'H', 6,  7, 3, 1, 0, 0},
    {'F', 9,  9, 1, 1, 0, 0},
    {'R', 6, 12, 4, 1, 0, 0},
    {'Y', 9,  9, 1, 2, 0, 0},
    {'W', 11, 10, 2, 1, 0, 0},
    {'U', 3,  5, 1, 1, 0, 1},  // selenocysteine
};

std::string formula_of(const int counts[6]) {
    std::string f;
    const char* const syms[6] = {"C", "H", "N", "O", "S", "Se"};
    for (int i = 0; i < 6; ++i)
        if (counts[i] > 0) f += std::string(syms[i]) + std::to_string(counts[i]);
    return f;
}

}  // namespace

TEST_CASE("every aminoacid code has the standard residue composition") {
    for (const Residue& r : kResidues) {
        INFO("aminoacid=" << r.code);
        int counts[6];
        const char seq[2] = {r.code, '\0'};
        parse_fasta(seq, counts);
        CHECK(counts[0] == r.C);
        CHECK(counts[1] == r.H);
        CHECK(counts[2] == r.N);
        CHECK(counts[3] == r.O);
        CHECK(counts[4] == r.S);
        CHECK(counts[5] == r.Se);

        // Lower case must resolve to the same residue.
        int lower[6];
        const char lseq[2] = {static_cast<char>(r.code + ('a' - 'A')), '\0'};
        parse_fasta(lseq, lower);
        for (int i = 0; i < 6; ++i) CHECK(lower[i] == counts[i]);
    }
}

TEST_CASE("residue counts are additive over a sequence") {
    const std::string seq = "MKWVTFISLLLLFSSAYSRGV";
    int counts[6];
    parse_fasta(seq.c_str(), counts);

    int expected[6] = {0, 0, 0, 0, 0, 0};
    for (char c : seq)
        for (const Residue& r : kResidues)
            if (r.code == c) {
                expected[0] += r.C; expected[1] += r.H; expected[2] += r.N;
                expected[3] += r.O; expected[4] += r.S; expected[5] += r.Se;
            }
    for (int i = 0; i < 6; ++i) CHECK(counts[i] == expected[i]);
}

TEST_CASE("parse_fasta_full adds the terminating water") {
    int residues[6], full[6];
    parse_fasta("PEPTIDE", residues);
    parse_fasta_full("PEPTIDE", full);
    CHECK(full[1] == residues[1] + 2);  // H
    CHECK(full[3] == residues[3] + 1);  // O
    for (int i : {0, 2, 4, 5}) CHECK(full[i] == residues[i]);
}

TEST_CASE("unrecognized characters are silently ignored") {
    int base[6];
    parse_fasta("AEDA", base);

    // Every variant holds the same multiset of residues (order is irrelevant to
    // the composition) plus characters of indeterminate formula, which drop out.
    for (const char* variant : {"AE-DA", "EAXXDA*", "AE DA", "ae\tda", "A E D A\n", "AEDA?!"}) {
        INFO("variant='" << variant << "'");
        int counts[6];
        parse_fasta(variant, counts);
        for (int i = 0; i < 6; ++i) CHECK(counts[i] == base[i]);
    }
}

TEST_CASE("empty sequence gives an empty composition") {
    int counts[6];
    parse_fasta("", counts);
    for (int i = 0; i < 6; ++i) CHECK(counts[i] == 0);

    parse_fasta_full("", counts);
    CHECK(counts[0] == 0);
    CHECK(counts[1] == 2);
    CHECK(counts[2] == 0);
    CHECK(counts[3] == 1);
    CHECK(counts[4] == 0);
    CHECK(counts[5] == 0);
}

TEST_CASE("Iso::FromFASTA matches the equivalent chemical formula") {
    for (const char* seq : {"A", "PEPTIDE", "MKWVTFISLLLLFSSAYSRGV", "GG", "ACDEFGHIKLMNPQRSTVWY"}) {
        INFO("sequence=" << seq);
        int with_water[6], without[6];
        parse_fasta_full(seq, with_water);
        parse_fasta(seq, without);

        Iso from_fasta = Iso::FromFASTA(seq);
        Iso from_formula(formula_of(with_water));
        CHECK(from_fasta.getMonoisotopicPeakMass() ==
              doctest::Approx(from_formula.getMonoisotopicPeakMass()));
        CHECK(from_fasta.getTheoreticalAverageMass() ==
              doctest::Approx(from_formula.getTheoreticalAverageMass()));
        CHECK(from_fasta.variance() == doctest::Approx(from_formula.variance()));

        Iso no_water = Iso::FromFASTA(seq, false, false);
        Iso no_water_formula(formula_of(without));
        CHECK(no_water.getMonoisotopicPeakMass() ==
              doctest::Approx(no_water_formula.getMonoisotopicPeakMass()));

        // The water is worth exactly one H2O.
        CHECK(from_fasta.getMonoisotopicPeakMass() - no_water.getMonoisotopicPeakMass() ==
              doctest::Approx(Iso("H2O1").getMonoisotopicPeakMass()));
    }
}

TEST_CASE("FromFASTA: std::string overload and nominal masses") {
    const std::string seq = "PEPTIDE";
    CHECK(Iso::FromFASTA(seq).getMonoisotopicPeakMass() ==
          Iso::FromFASTA(seq.c_str()).getMonoisotopicPeakMass());

    Iso nominal = Iso::FromFASTA(seq, true, true);
    const double nm = nominal.getMonoisotopicPeakMass();
    CHECK(nm == doctest::Approx(std::round(nm)));
    CHECK(nm == doctest::Approx(799.0));  // PEPTIDE, nominal monoisotopic mass
}

TEST_CASE("selenocysteine adds a sixth dimension") {
    Iso without = Iso::FromFASTA("PEPTIDE");
    Iso with = Iso::FromFASTA("PEPTIDU");
    CHECK(without.getDimNumber() == 5);
    CHECK(with.getDimNumber() == 6);
    CHECK(with.getAllDim() > without.getDimNumber());

    // Selenium is polyisotopic, so the U-containing peptide has a much richer
    // fine structure.
    FixedEnvelope e_without = FixedEnvelope::FromThreshold(std::move(without), 1e-6, true, false);
    FixedEnvelope e_with = FixedEnvelope::FromThreshold(std::move(with), 1e-6, true, false);
    CHECK(e_with.confs_no() > e_without.confs_no());
}

TEST_CASE("empty FASTA sequence yields plain water") {
    Iso water_from_fasta = Iso::FromFASTA("");
    Iso water("H2O1");
    CHECK(water_from_fasta.getMonoisotopicPeakMass() ==
          doctest::Approx(water.getMonoisotopicPeakMass()));
    CHECK(water_from_fasta.getTheoreticalAverageMass() ==
          doctest::Approx(water.getTheoreticalAverageMass()));

    // With no water added, an empty sequence is an empty (zero-mass) molecule.
    Iso nothing = Iso::FromFASTA("", false, false);
    CHECK(nothing.getMonoisotopicPeakMass() == doctest::Approx(0.0));
    FixedEnvelope env = FixedEnvelope::FromThreshold(std::move(nothing), 0.0, true, false);
    CHECK(env.confs_no() == 1);
    CHECK(env.probs()[0] == doctest::Approx(1.0));
}

TEST_CASE("known peptide monoisotopic masses") {
    // Reference values from standard peptide mass calculators (monoisotopic,
    // neutral, with terminal water).
    struct Ref { const char* seq; double mass; };
    const Ref refs[] = {
        {"G",        75.03203},
        {"GG",      132.05349},
        {"A",        89.04768},
        {"PEPTIDE", 799.35977},
        {"YGGFL",   555.26935},  // leu-enkephalin
    };
    for (const Ref& r : refs) {
        INFO("sequence=" << r.seq);
        CHECK(Iso::FromFASTA(r.seq).getMonoisotopicPeakMass() ==
              doctest::Approx(r.mass).epsilon(1e-6));
    }
}
