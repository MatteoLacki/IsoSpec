"""The IsoSpecPy generator classes.

The generators stream isotopologues instead of materializing them, so what
matters is the iteration contract (single use, correct stopping rule, matching
configurations) and that they agree with the materialized classes.
"""

import math

import pytest

import IsoSpecPy


def collect(generator):
    return [(mass, prob) for mass, prob in generator]


def merged(pairs, tol=1e-6):
    """Peaks sorted by mass, with near-degenerate masses merged.

    The generators sum the same terms in different orders, so masses can differ
    in the last bits; two isotopologues whose masses coincide to ~1e-13 can even
    come out in a different order.  Merging on a tolerance makes the comparison
    independent of both.
    """
    out = []
    for mass, prob in sorted(pairs):
        if out and abs(out[-1][0] - mass) <= tol:
            out[-1][1] += prob
        else:
            out.append([mass, prob])
    return out


def assert_peaks_close(a, b):
    a, b = merged(a), merged(b)
    assert len(a) == len(b)
    for (m1, p1), (m2, p2) in zip(a, b):
        assert math.isclose(m1, m2, rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(p1, p2, rel_tol=1e-9, abs_tol=1e-15)


# --------------------------------------------------------------------------
# Iteration contract
# --------------------------------------------------------------------------

def test_generators_are_single_use():
    for generator in (IsoSpecPy.IsoThresholdGenerator(0.01, "C6H12O6"),
                      IsoSpecPy.IsoLayeredGenerator("C6H12O6"),
                      IsoSpecPy.IsoOrderedGenerator("C6H12O6"),
                      IsoSpecPy.IsoStochasticGenerator(1000, "C6H12O6")):
        assert len(collect(generator)) > 0
        with pytest.raises(NotImplementedError):
            list(generator)


def test_threshold_generator_respects_its_threshold():
    absolute = IsoSpecPy.IsoThresholdGenerator(1e-4, "C50H100N20O20S5", absolute=True)
    peaks = collect(absolute)
    assert len(peaks) > 0
    assert all(prob >= 1e-4 for _, prob in peaks)

    # It enumerates exactly what the materialized version contains.
    materialized = IsoSpecPy.IsoThreshold(1e-4, "C50H100N20O20S5", absolute=True)
    assert len(peaks) == len(materialized)
    assert_peaks_close(peaks, zip(materialized.masses, materialized.probs))


def test_threshold_generator_relative_threshold():
    relative = collect(IsoSpecPy.IsoThresholdGenerator(1e-4, "C50H100N20O20S5",
                                                       absolute=False))
    absolute = collect(IsoSpecPy.IsoThresholdGenerator(1e-4, "C50H100N20O20S5",
                                                       absolute=True))
    assert len(relative) >= len(absolute)

    mode = math.exp(IsoSpecPy.Iso("C50H100N20O20S5").getModeLProb())
    assert all(prob >= 1e-4 * mode * (1 - 1e-9) for _, prob in relative)


def test_layered_generator_enumerates_everything():
    peaks = collect(IsoSpecPy.IsoLayeredGenerator("C10H20O2"))
    assert math.isclose(sum(p for _, p in peaks), 1.0, rel_tol=1e-9)

    # Same support as the threshold generator run with no threshold.
    threshold = collect(IsoSpecPy.IsoThresholdGenerator(0.0, "C10H20O2", absolute=True))
    assert len(peaks) == len(threshold)
    assert_peaks_close(peaks, threshold)


def test_ordered_generator_is_descending():
    probs = [prob for _, prob in IsoSpecPy.IsoOrderedGenerator("C6H12O6")]
    assert probs == sorted(probs, reverse=True)
    assert math.isclose(sum(probs), 1.0, rel_tol=1e-9)

    # The first peak is the mode.
    first = next(iter(IsoSpecPy.IsoOrderedGenerator("C6H12O6")))
    assert math.isclose(math.log(first[1]),
                        IsoSpecPy.Iso("C6H12O6").getModeLProb(), rel_tol=1e-9)


def test_stochastic_generator_counts():
    for n in (1, 1000, 100000):
        peaks = collect(IsoSpecPy.IsoStochasticGenerator(n, "C100H200N40O40"))
        assert len(peaks) > 0
        assert all(float(count).is_integer() and count >= 1 for _, count in peaks)
        assert math.isclose(sum(count for _, count in peaks), float(n))


def test_stochastic_generator_stays_in_the_mass_range():
    iso = IsoSpecPy.Iso("C100H200N40O40")
    lo, hi = iso.getLightestPeakMass(), iso.getHeaviestPeakMass()
    peaks = collect(IsoSpecPy.IsoStochasticGenerator(50000, "C100H200N40O40"))
    assert all(lo - 1e-9 <= mass <= hi + 1e-9 for mass, _ in peaks)

    total = sum(count for _, count in peaks)
    mean = sum(mass * count for mass, count in peaks) / total
    assert math.isclose(mean, iso.getTheoreticalAverageMass(), rel_tol=1e-4)


# --------------------------------------------------------------------------
# Configurations
# --------------------------------------------------------------------------

@pytest.mark.parametrize("factory", [
    lambda: IsoSpecPy.IsoThresholdGenerator(1e-6, "C6H12O6", get_confs=True),
    lambda: IsoSpecPy.IsoLayeredGenerator("C6H12O6", get_confs=True),
    lambda: IsoSpecPy.IsoOrderedGenerator("C6H12O6", get_confs=True),
])
def test_generator_configurations(factory):
    seen = set()
    count = 0
    for mass, prob, conf in factory():
        assert [sum(element) for element in conf] == [6, 12, 6]
        assert all(c >= 0 for element in conf for c in element)
        # Every configuration is reported exactly once.
        assert conf not in seen
        seen.add(conf)
        count += 1
    assert count > 10


def test_generator_configuration_masses_are_consistent():
    from IsoSpecPy import PeriodicTbl
    masses_per_element = [PeriodicTbl.symbol_to_masses[s] for s in ("C", "H", "O")]

    for mass, prob, conf in IsoSpecPy.IsoThresholdGenerator(1e-5, "C6H12O6",
                                                            get_confs=True):
        expected = sum(count * m
                       for element_masses, element_conf in zip(masses_per_element, conf)
                       for m, count in zip(element_masses, element_conf))
        assert math.isclose(mass, expected, rel_tol=1e-12)


# --------------------------------------------------------------------------
# Cross-checks
# --------------------------------------------------------------------------

@pytest.mark.parametrize("formula", ["H2O1", "C6H12O6", "S2Cl2", "Sn1", "F10"])
def test_all_generators_agree_on_the_full_support(formula):
    threshold = collect(IsoSpecPy.IsoThresholdGenerator(0.0, formula, absolute=True))
    layered = collect(IsoSpecPy.IsoLayeredGenerator(formula))
    ordered = collect(IsoSpecPy.IsoOrderedGenerator(formula))

    assert len(threshold) > 0
    assert len(threshold) == len(layered) == len(ordered)
    assert_peaks_close(threshold, layered)
    assert_peaks_close(threshold, ordered)


def test_generator_accepts_fasta_and_raw_arrays():
    from_fasta = collect(IsoSpecPy.IsoThresholdGenerator(1e-6, fasta="PEPTIDE"))
    assert len(from_fasta) > 0

    raw = collect(IsoSpecPy.IsoThresholdGenerator(
        0.0, atomCounts=[2], isotopeMasses=[[10.0, 11.0]],
        isotopeProbabilities=[[0.7, 0.3]], absolute=True))
    assert len(raw) == 3
    assert math.isclose(sum(p for _, p in raw), 1.0, rel_tol=1e-12)
    # Binomial(2, 0.3): 0.49, 0.42, 0.09.
    assert sorted(round(p, 12) for _, p in raw) == [0.09, 0.42, 0.49]


def test_charge_is_applied_to_generated_masses():
    neutral = sorted(mass for mass, _ in
                     IsoSpecPy.IsoThresholdGenerator(1e-6, "C6H12O6"))
    doubly = sorted(mass for mass, _ in
                    IsoSpecPy.IsoThresholdGenerator(1e-6, "C6H12O6", charge=2.0))
    assert len(neutral) == len(doubly)
    assert all(math.isclose(n / 2.0, d, rel_tol=1e-12) for n, d in zip(neutral, doubly))
