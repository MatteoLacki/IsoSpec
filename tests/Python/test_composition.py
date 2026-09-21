"""IsoParamsFromDict now resolves element symbols to isotopes in C++
(cwrapper.h's expandCompositionC -> fasta_mods.h's expand_composition_into)
rather than walking PeriodicTbl itself.

These tests pin two separate things: that the values are unchanged from the
PeriodicTbl oracle they used to come from, and that the Python-side copy is
genuinely gone rather than merely unused on the happy path.
"""

import pytest

import IsoSpecPy
from IsoSpecPy import PeriodicTbl
from IsoSpecPy.IsoSpecPy import IsoParamsFromDict, IsoParamsFromFormula


def test_matches_the_periodic_table_oracle():
    # Every element the library knows, one at a time: the C++ expansion must
    # reproduce exactly what PeriodicTbl (built from the same tables, but
    # assembled independently in Python) says.
    for symbol in sorted(PeriodicTbl.symbol_to_masses):
        parsed = IsoParamsFromDict({symbol: 1})
        assert parsed.elems == [symbol]
        assert parsed.atomCounts == [1]
        assert parsed.masses == [PeriodicTbl.symbol_to_masses[symbol]], symbol
        assert parsed.probs == [PeriodicTbl.symbol_to_probs[symbol]], symbol


def test_nominal_masses_match_the_oracle():
    for symbol in sorted(PeriodicTbl.symbol_to_massNo):
        parsed = IsoParamsFromDict({symbol: 1}, use_nominal_masses=True)
        assert parsed.masses == [PeriodicTbl.symbol_to_massNo[symbol]], symbol
        # Nominal masses are nucleon counts; probabilities are untouched by
        # the mass table choice.
        assert all(float(m).is_integer() for m in parsed.masses[0]), symbol
        assert parsed.probs == [PeriodicTbl.symbol_to_probs[symbol]], symbol


def test_preserves_dict_order_and_shape():
    parsed = IsoParamsFromDict({"O": 1, "H": 2, "Se": 3})

    assert parsed.elems == ["O", "H", "Se"]
    assert parsed.atomCounts == [1, 2, 3]
    # Per-element tuples, not one flat sequence: Advanced.py and
    # approximations.py both index masses/probs by element.
    assert [len(m) for m in parsed.masses] == [3, 2, 6]
    assert [len(p) for p in parsed.probs] == [3, 2, 6]
    assert all(isinstance(m, tuple) for m in parsed.masses)
    assert all(isinstance(p, tuple) for p in parsed.probs)


def test_resolution_no_longer_reads_periodic_tbl(monkeypatch):
    # The real regression guard. If a Python-side copy of the symbol ->
    # isotopes lookup ever comes back, emptying PeriodicTbl's dicts will
    # break this and nothing else will.
    monkeypatch.setattr(PeriodicTbl, "symbol_to_masses", {})
    monkeypatch.setattr(PeriodicTbl, "symbol_to_massNo", {})
    monkeypatch.setattr(PeriodicTbl, "symbol_to_probs", {})

    parsed = IsoParamsFromDict({"H": 2, "O": 1})
    assert parsed.atomCounts == [2, 1]
    assert [len(m) for m in parsed.masses] == [2, 3]
    assert parsed.masses[0][0] == pytest.approx(1.00782503227)

    # ...and an Iso built from a dict still works with PeriodicTbl gutted.
    iso = IsoSpecPy.Iso(formula={"H": 2, "O": 1})
    assert iso.getMonoisotopicPeakMass() == pytest.approx(18.0105646, abs=1e-6)


def test_unknown_symbol_raises_value_error():
    for bad in ("Xx", "Unobtainium", "", "c", "H2"):
        with pytest.raises(ValueError):
            IsoParamsFromDict({bad: 1})


def test_non_ascii_symbol_raises_value_error():
    with pytest.raises(ValueError):
        IsoParamsFromDict({"Ćw": 1})


def test_empty_composition():
    parsed = IsoParamsFromDict({})
    assert parsed.atomCounts == []
    assert parsed.masses == []
    assert parsed.probs == []
    assert parsed.elems == []


def test_negative_counts_pass_through_but_iso_rejects_them():
    # A composition is not a molecule: the expansion carries a negative count
    # through unchanged, and Iso() is what refuses it.
    parsed = IsoParamsFromDict({"H": -1, "S": -1})
    assert parsed.atomCounts == [-1, -1]
    assert [len(m) for m in parsed.masses] == [2, 4]

    with pytest.raises(Exception):
        IsoSpecPy.Iso(formula={"H": -1, "S": -1})


def test_formula_string_path_unchanged():
    # IsoParamsFromFormula goes through the same expansion; approximations.py
    # consumes its .probs/.elems directly.
    parsed = IsoParamsFromFormula("C6H12N2O1S1")
    assert parsed.elems == ["C", "H", "N", "O", "S"]
    assert parsed.atomCounts == [6, 12, 2, 1, 1]
    assert [len(p) for p in parsed.probs] == [2, 2, 2, 3, 4]


def test_end_to_end_masses_are_unchanged():
    # Belt and braces: the whole point is that nothing observable moved.
    for formula, mono in (
        ("H2O1", 18.0105646),
        ("C6H12O6", 180.0633881),
        ("C8H10N4O2", 194.0803756),
    ):
        iso = IsoSpecPy.Iso(formula=formula)
        assert iso.getMonoisotopicPeakMass() == pytest.approx(mono, abs=1e-5), formula


def test_peptide_sequence_plus_formula_still_composes():
    # The case that kept Python from simply forwarding to isoFromFastaWithMods:
    # a sequence composed with an extra formula, which only works because the
    # composition is merged in Python and expanded afterwards.
    seq = IsoSpecPy.Iso(peptide_sequence="PEPTIDE", formula="H2O1")
    residues = IsoSpecPy.Iso(peptide_sequence="PEPTIDE")
    assert seq.getMonoisotopicPeakMass() == pytest.approx(
        residues.getMonoisotopicPeakMass() + 18.0105646, abs=1e-6
    )


def test_repeated_calls_do_not_accumulate():
    # The C handle is freed on every call (finally: deleteExpandedCompositionC).
    # Nothing here can see a leak directly; what it can see is a wrong answer
    # from a handle reused or a buffer appended to across calls.
    first = IsoParamsFromDict({"C": 1})
    for _ in range(1000):
        again = IsoParamsFromDict({"C": 1})
        assert again == first
