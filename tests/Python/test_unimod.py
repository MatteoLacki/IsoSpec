"""[UNIMOD:<id>] modification-aware peptide sequence parsing.

ParsePeptideSequence / ParseFASTA (the compatibility alias) / Iso's
peptide_sequence=/fasta= kwargs -- see docs/ai/unimod.md for the notation,
the shipped table's exclusions, and the fasta/peptide_sequence naming note.
"""

import math
import os

import pytest

import IsoSpecPy
from IsoSpecPy.IsoSpecPy import ParseFASTA, ParsePeptideSequence


def test_parse_peptide_sequence_matches_parse_fasta_for_plain_sequence():
    assert ParsePeptideSequence("PEPTIDE") == ParseFASTA("PEPTIDE")


def test_peptide_sequence_kwarg_is_alias_for_fasta_kwarg():
    from_fasta = IsoSpecPy.Iso(fasta="PEPTIDE")
    from_peptide_sequence = IsoSpecPy.Iso(peptide_sequence="PEPTIDE")
    assert math.isclose(from_fasta.getMonoisotopicPeakMass(),
                        from_peptide_sequence.getMonoisotopicPeakMass())


def test_fasta_and_peptide_sequence_kwargs_conflict_raises():
    with pytest.raises(ValueError):
        IsoSpecPy.Iso(fasta="PEPTIDE", peptide_sequence="OTHERSEQ")


def test_unimod_bracket_regression_reproduces_this_sessions_manual_check():
    # Before this feature, every letter in "UNIMOD" (U, N, I, M, O, D) was a
    # valid 1-letter amino-acid code, so a [UNIMOD:<id>] bracket was silently
    # misread as six phantom residues instead of a modification. Confirmed
    # empirically against the unfixed code this session:
    #   ParseFASTA('MPEPTC[UNIMOD:4]DEK') -> {'C': 76, 'H': 123, 'N': 19,
    #                                          'O': 27, 'S': 3, 'Se': 1}
    # (wrong) instead of the base composition plus Carbamidomethyl's delta.
    bracketed = ParseFASTA("MPEPTC[UNIMOD:4]DEK")
    base = ParseFASTA("MPEPTCDEK")

    assert "Se" not in bracketed  # the old bug's tell: phantom Se from misread "U"

    carbamidomethyl_delta = {"H": 3, "C": 2, "N": 1, "O": 1}
    expected = dict(base)
    for symbol, delta in carbamidomethyl_delta.items():
        expected[symbol] = expected.get(symbol, 0) + delta
    assert dict(bracketed) == expected


def test_formula_and_peptide_sequence_with_mods_combine():
    both = IsoSpecPy.Iso(formula="H2O1", peptide_sequence="PEPTC[UNIMOD:4]DEK")
    seq_only = IsoSpecPy.Iso(peptide_sequence="PEPTC[UNIMOD:4]DEK")
    water = IsoSpecPy.Iso("H2O1")
    assert math.isclose(
        both.getMonoisotopicPeakMass(),
        seq_only.getMonoisotopicPeakMass() + water.getMonoisotopicPeakMass())


def test_unknown_unimod_id_raises():
    with pytest.raises(ValueError):
        ParsePeptideSequence("PEPTC[UNIMOD:999999999]DEK")


def test_excluded_unimod_id_raises_same_as_unknown():
    # UNIMOD:9, ICAT-G:2H(8) -- isotope-labeled, deliberately excluded from
    # the shipped table (see scripts/build_unimod_table.py).
    with pytest.raises(ValueError):
        ParsePeptideSequence("PEPTC[UNIMOD:9]DEK")


def test_malformed_bracket_raises():
    with pytest.raises(ValueError):
        ParsePeptideSequence("PEPTC[UNIMOD:4DEK")


def test_unimod_db_path_override(tmp_path):
    override = tmp_path / "custom_unimod.csv"
    override.write_text("id,name,mono_mass,composition\n4,TestOverride,1.0,H1\n")

    overridden = ParsePeptideSequence("PEPTC[UNIMOD:4]DEK", unimod_db_path=str(override))
    base = ParseFASTA("PEPTCDEK")
    expected = dict(base)
    expected["H"] = expected.get("H", 0) + 1
    assert dict(overridden) == expected

    # The embedded default is untouched by having loaded an override.
    still_default = ParsePeptideSequence("PEPTC[UNIMOD:4]DEK")
    assert still_default["C"] - base["C"] == 2  # Carbamidomethyl's real C delta, not the override's 0


def test_id_past_32_bits_is_unknown_not_truncated():
    # The id is parsed as 64-bit but looked up as unsigned int; without a
    # range check 2**32 + 4 wraps to 4 and silently applies Carbamidomethyl.
    for bad_id in (2**32 + 4, 2**33 + 4):
        with pytest.raises(ValueError):
            ParsePeptideSequence("PEPTC[UNIMOD:{}]DEK".format(bad_id))
    # ...and the id it would have wrapped onto is unaffected.
    assert ParsePeptideSequence("PEPTC[UNIMOD:4]DEK")["C"] - ParseFASTA("PEPTCDEK")["C"] == 2


def test_zero_count_elements_are_omitted():
    # Documented behaviour change from the pre-Unimod ParseFASTA, which
    # always emitted C/H/N/O/S keys: only non-zero counts appear now.
    composition = ParseFASTA("A")
    assert "S" not in composition
    assert composition["C"] == 3


def test_negative_net_count_is_reported_by_the_parser_and_refused_by_Iso():
    # Met->Hsl (UNIMOD:11) removes H4 C1 S1; a lone glycine cannot pay for it.
    composition = ParsePeptideSequence("G[UNIMOD:11]")
    assert composition["H"] < 0
    with pytest.raises(Exception):
        IsoSpecPy.Iso(peptide_sequence="G[UNIMOD:11]")


def test_charge_keeps_its_positional_slot():
    # peptide_sequence/unimod_db_path were added after charge precisely so
    # this keeps working: formula, get_confs, atomCounts, isotopeMasses,
    # isotopeProbabilities, use_nominal_masses, fasta, charge.
    positional = IsoSpecPy.Iso("H2O1", False, None, None, None, False, "", 2.0)
    keyword = IsoSpecPy.Iso(formula="H2O1", charge=2.0)
    assert math.isclose(positional.getMonoisotopicPeakMass(),
                        keyword.getMonoisotopicPeakMass())
    assert math.isclose(keyword.getMonoisotopicPeakMass(),
                        IsoSpecPy.Iso(formula="H2O1").getMonoisotopicPeakMass() / 2.0)
