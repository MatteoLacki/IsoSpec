"""The IsoSpecPy.Iso layer: construction, validation and the molecule accessors.

These check the Python wrapper itself — formula parsing, the several ways of
describing a molecule, the charge/nominal-mass options and the input validation
— against values computed from IsoSpecPy.PeriodicTbl rather than from another
IsoSpec call.
"""

import math

import pytest

import IsoSpecPy
from IsoSpecPy import PeriodicTbl
from IsoSpecPy.IsoSpecPy import ParseFASTA, ParseFormula, IsoParamsFromFormula


def element_isotopes(symbol):
    return list(zip(PeriodicTbl.symbol_to_masses[symbol],
                    PeriodicTbl.symbol_to_probs[symbol]))


def expected_monoisotopic(formula_counts):
    total = 0.0
    for symbol, count in formula_counts.items():
        isotopes = element_isotopes(symbol)
        total += count * max(isotopes, key=lambda mp: mp[1])[0]
    return total


def expected_average(formula_counts):
    total = 0.0
    for symbol, count in formula_counts.items():
        total += count * sum(m * p for m, p in element_isotopes(symbol))
    return total


# --------------------------------------------------------------------------
# Formula parsing
# --------------------------------------------------------------------------

def test_parse_formula_counts():
    assert dict(ParseFormula("C6H12O6")) == {"C": 6, "H": 12, "O": 6}
    # A missing count means one atom.
    assert dict(ParseFormula("CH4")) == {"C": 1, "H": 4}
    assert dict(ParseFormula("H2O")) == {"H": 2, "O": 1}
    # Two-letter symbols and large counts.
    assert dict(ParseFormula("Se2Sn10Cl100")) == {"Se": 2, "Sn": 10, "Cl": 100}
    # Order is preserved.
    assert list(ParseFormula("O1C1H1").keys()) == ["O", "C", "H"]


@pytest.mark.parametrize("formula", [
    "",             # empty
    "C6H12O6Xx1",   # unknown element
    "6C",           # leading garbage
    "C6 H12",       # space
    "C6H12O6!",     # trailing garbage
    "c6h12o6",      # wrong case
    "C6C6",         # repeated element
])
def test_parse_formula_rejects_bad_input(formula):
    with pytest.raises(ValueError):
        ParseFormula(formula)


def test_iso_params_from_formula():
    counts, masses, probs, elems = IsoParamsFromFormula("C2O1")
    assert counts == [2, 1]
    assert elems == ["C", "O"]
    assert len(masses) == len(probs) == 2
    assert len(masses[0]) == 2   # carbon
    assert len(masses[1]) == 3   # oxygen
    assert math.isclose(sum(probs[0]), 1.0, rel_tol=1e-9)

    # Nominal masses are integers.
    _, nominal, _, _ = IsoParamsFromFormula("C2O1", use_nominal_masses=True)
    assert all(float(m).is_integer() for element in nominal for m in element)


# --------------------------------------------------------------------------
# Iso construction
# --------------------------------------------------------------------------

def test_iso_accessors_match_the_periodic_table():
    counts = {"C": 6, "H": 12, "O": 6}
    iso = IsoSpecPy.Iso("C6H12O6")

    assert math.isclose(iso.getMonoisotopicPeakMass(), expected_monoisotopic(counts),
                        rel_tol=1e-12)
    assert math.isclose(iso.getTheoreticalAverageMass(), expected_average(counts),
                        rel_tol=1e-12)
    assert iso.getLightestPeakMass() <= iso.getMonoisotopicPeakMass()
    assert iso.getMonoisotopicPeakMass() <= iso.getHeaviestPeakMass()
    assert iso.variance() > 0.0
    assert math.isclose(iso.stddev(), math.sqrt(iso.variance()), rel_tol=1e-12)

    # Log-probabilities are negative and bounded by the mode.
    for lprob in (iso.getLightestPeakLProb(), iso.getHeaviestPeakLProb(),
                  iso.getMonoisotopicPeakLProb()):
        assert lprob <= iso.getModeLProb() + 1e-12
    assert iso.getModeLProb() <= 0.0


def test_iso_peak_configurations():
    iso = IsoSpecPy.Iso("C6H12O6")
    for conf in (iso.getLightestPeakConf(), iso.getHeaviestPeakConf(),
                 iso.getMonoisotopicPeakConf()):
        # One tuple per element, summing to that element's atom count.
        assert len(conf) == 3
        assert [sum(element) for element in conf] == [6, 12, 6]
        assert all(count >= 0 for element in conf for count in element)

    # The monoisotopic configuration puts every atom on the commonest isotope.
    mono = iso.getMonoisotopicPeakConf()
    assert mono == ((6, 0), (12, 0), (6, 0, 0))


def test_iso_from_dict_and_from_string_agree():
    from_string = IsoSpecPy.Iso("C6H12O6")
    from_dict = IsoSpecPy.Iso({"C": 6, "H": 12, "O": 6})
    assert math.isclose(from_string.getMonoisotopicPeakMass(),
                        from_dict.getMonoisotopicPeakMass())
    assert math.isclose(from_string.variance(), from_dict.variance())


def test_iso_from_fasta_matches_formula():
    counts = ParseFASTA("PEPTIDE")
    # ParseFASTA returns the residue skeleton; Iso(fasta=...) adds nothing more,
    # so the two must describe the same molecule.
    from_fasta = IsoSpecPy.Iso(fasta="PEPTIDE")
    from_counts = IsoSpecPy.Iso(counts)
    assert math.isclose(from_fasta.getMonoisotopicPeakMass(),
                        from_counts.getMonoisotopicPeakMass())

    # Selenocysteine brings in a sixth element.
    assert "Se" not in ParseFASTA("PEPTIDE")
    assert ParseFASTA("PEPTIDU")["Se"] == 1


def test_iso_fasta_and_formula_combine():
    both = IsoSpecPy.Iso(formula="H2O1", fasta="PEPTIDE")
    fasta_only = IsoSpecPy.Iso(fasta="PEPTIDE")
    water = IsoSpecPy.Iso("H2O1")
    assert math.isclose(
        both.getMonoisotopicPeakMass(),
        fasta_only.getMonoisotopicPeakMass() + water.getMonoisotopicPeakMass())


def test_iso_from_raw_arrays():
    # Two fictitious elements, given explicitly.
    iso = IsoSpecPy.Iso(atomCounts=[2, 1],
                        isotopeMasses=[[10.0, 11.0], [20.0, 21.0, 22.0]],
                        isotopeProbabilities=[[0.7, 0.3], [0.5, 0.3, 0.2]])
    assert math.isclose(iso.getMonoisotopicPeakMass(), 2 * 10.0 + 20.0)
    assert math.isclose(iso.getLightestPeakMass(), 2 * 10.0 + 20.0)
    assert math.isclose(iso.getHeaviestPeakMass(), 2 * 11.0 + 22.0)
    assert math.isclose(iso.getTheoreticalAverageMass(),
                        2 * (10.0 * 0.7 + 11.0 * 0.3) +
                        (20.0 * 0.5 + 21.0 * 0.3 + 22.0 * 0.2))

    # The distribution is the corresponding product multinomial: 3 * 6 configs.
    dist = IsoSpecPy.IsoThreshold(0.0, atomCounts=[2, 1],
                                  isotopeMasses=[[10.0, 11.0], [20.0, 21.0, 22.0]],
                                  isotopeProbabilities=[[0.7, 0.3], [0.5, 0.3, 0.2]])
    assert len(dist) == 3 * 3
    assert math.isclose(dist.total_prob(), 1.0, rel_tol=1e-12)


def test_iso_nominal_masses():
    real = IsoSpecPy.Iso("C6H12O6")
    nominal = IsoSpecPy.Iso("C6H12O6", use_nominal_masses=True)
    assert math.isclose(nominal.getMonoisotopicPeakMass(), 180.0)
    assert not math.isclose(real.getMonoisotopicPeakMass(), 180.0, abs_tol=1e-6)
    # Probabilities are unaffected.
    assert math.isclose(real.getModeLProb(), nominal.getModeLProb())


def test_iso_charge_divides_masses():
    neutral = IsoSpecPy.Iso("C6H12O6")
    for charge in (2.0, 3.0):
        charged = IsoSpecPy.Iso("C6H12O6", charge=charge)
        assert math.isclose(charged.getMonoisotopicPeakMass(),
                            neutral.getMonoisotopicPeakMass() / charge)
        assert math.isclose(charged.getTheoreticalAverageMass(),
                            neutral.getTheoreticalAverageMass() / charge)
        # Probabilities are unchanged by the m/z conversion.
        assert math.isclose(charged.getModeLProb(), neutral.getModeLProb())


def test_marginal_log_size_estimates():
    iso = IsoSpecPy.Iso("C100H2O1")
    estimates = iso.getMarginalLogSizeEstimates(0.99)
    assert len(estimates) == 3
    assert all(not math.isnan(e) for e in estimates)
    # C100 dominates H2 and O1.
    assert estimates[0] > estimates[1]
    assert estimates[0] > estimates[2]


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def test_iso_requires_a_molecule():
    with pytest.raises(Exception):
        IsoSpecPy.Iso()


def test_iso_rejects_negative_counts():
    with pytest.raises(Exception):
        IsoSpecPy.Iso("C-1H2")


@pytest.mark.parametrize("probs", [
    [[0.0, 1.0]],       # zero probability
    [[-0.5, 1.5]],      # negative probability
    [[1.5, 0.5]],       # probability above one
])
def test_iso_rejects_invalid_probabilities(probs):
    with pytest.raises(ValueError):
        IsoSpecPy.Iso(atomCounts=[2], isotopeMasses=[[10.0, 11.0]],
                      isotopeProbabilities=probs)


def test_iso_rejects_mismatched_arrays():
    with pytest.raises(ValueError):
        # Two masses, three probabilities.
        IsoSpecPy.Iso(atomCounts=[2], isotopeMasses=[[10.0, 11.0]],
                      isotopeProbabilities=[[0.5, 0.3, 0.2]])

    with pytest.raises(ValueError):
        # Two elements' worth of masses, one atom count.
        IsoSpecPy.Iso(atomCounts=[2], isotopeMasses=[[10.0, 11.0], [20.0, 21.0]],
                      isotopeProbabilities=[[0.5, 0.5], [0.5, 0.5]])


def test_zero_atom_counts_are_allowed():
    iso = IsoSpecPy.Iso("C0H2")
    water_free = IsoSpecPy.Iso("H2")
    assert math.isclose(iso.getMonoisotopicPeakMass(),
                        water_free.getMonoisotopicPeakMass())
