"""The IsoSpecPy.IsoDistribution layer: the materialized-spectrum object.

Construction (all four factory functions plus raw masses/probs), the container
protocol, arithmetic, in-place transforms and the statistics.  Small hand-built
spectra are used wherever the answer can be written down exactly.
"""

import math

import pytest

import IsoSpecPy
from IsoSpecPy import IsoDistribution


def peaks(dist):
    return sorted(zip(dist.masses, dist.probs))


# --------------------------------------------------------------------------
# Construction
# --------------------------------------------------------------------------

def test_from_masses_and_probs():
    d = IsoDistribution(masses=[1.0, 2.0, 3.0], probs=[0.2, 0.3, 0.5])
    assert len(d) == 3
    assert math.isclose(d.total_prob(), 1.0)
    assert list(d.masses) == [1.0, 2.0, 3.0]
    assert list(d.probs) == [0.2, 0.3, 0.5]


def test_construction_errors():
    with pytest.raises(ValueError):
        IsoDistribution()
    with pytest.raises(ValueError):
        IsoDistribution(masses=[1.0, 2.0])
    with pytest.raises(ValueError):
        IsoDistribution(probs=[1.0, 2.0])
    with pytest.raises(ValueError):
        IsoDistribution(masses=[1.0, 2.0], probs=[1.0])


def test_container_protocol():
    d = IsoDistribution(masses=[1.0, 2.0], probs=[0.25, 0.75])
    assert len(d) == 2
    assert d[0] == (1.0, 0.25)
    assert d[1] == (2.0, 0.75)
    assert list(d) == [(1.0, 0.25), (2.0, 0.75)]


def test_iso_threshold_absolute_and_relative():
    absolute = IsoSpecPy.IsoThreshold(1e-4, "C50H100N20O20S5", absolute=True)
    relative = IsoSpecPy.IsoThreshold(1e-4, "C50H100N20O20S5", absolute=False)

    assert all(p >= 1e-4 for p in absolute.probs)
    # The relative threshold is a fraction of the highest peak, a lower bar.
    assert len(relative) >= len(absolute)

    # A threshold of zero enumerates everything.
    full = IsoSpecPy.IsoThreshold(0.0, "C6H12O6")
    assert math.isclose(full.total_prob(), 1.0, rel_tol=1e-9)


def test_iso_total_prob_covers_the_target():
    for target in (0.5, 0.9, 0.99, 0.999):
        d = IsoSpecPy.IsoTotalProb(target, "C100H200N40O40")
        assert d.total_prob() >= target - 1e-9

    # Trimming to the minimal peak set can only shrink the result.
    trimmed = IsoSpecPy.IsoTotalProb(0.99, "C100H200N40O40", get_minimal_pset=True)
    untrimmed = IsoSpecPy.IsoTotalProb(0.99, "C100H200N40O40", get_minimal_pset=False)
    assert len(trimmed) <= len(untrimmed)
    assert trimmed.total_prob() >= 0.99 - 1e-9


def test_iso_stochastic_counts():
    for n in (1, 100, 100000):
        d = IsoSpecPy.IsoStochastic(n, "C100H200N40O40")
        assert len(d) > 0
        assert all(float(p).is_integer() and p > 0 for p in d.probs)
        assert math.isclose(sum(d.probs), float(n))


def test_iso_binned():
    d = IsoSpecPy.IsoBinned(1.0, "C6H12O6", target_total_prob=0.99)
    assert len(d) > 0
    assert d.total_prob() >= 0.99
    assert all(math.isclose(m, round(m)) for m in d.masses)
    assert list(d.masses) == sorted(d.masses)

    # A non-zero bin middle shifts the grid.
    shifted = IsoSpecPy.IsoBinned(1.0, "C6H12O6", target_total_prob=0.99, bin_middle=0.5)
    assert all(math.isclose(m - 0.5, round(m - 0.5)) for m in shifted.masses)


def test_get_confs():
    d = IsoSpecPy.IsoTotalProb(0.99, "C6H12O6", get_confs=True)
    assert len(d) > 0
    for mass, prob, conf in d:
        assert [sum(element) for element in conf] == [6, 12, 6]
        assert prob > 0.0
        assert mass > 0.0
    # Indexing gives the configuration too.
    assert len(d[0]) == 3


# --------------------------------------------------------------------------
# Transforms
# --------------------------------------------------------------------------

def test_copy_is_independent():
    original = IsoSpecPy.IsoTotalProb(0.99, "C6H12O6")
    duplicate = original.copy()
    assert len(duplicate) == len(original)
    assert peaks(duplicate) == peaks(original)

    duplicate.add_mass(100.0)
    assert not math.isclose(list(duplicate.masses)[0], list(original.masses)[0])


def test_scale_and_normalize():
    d = IsoDistribution(masses=[1.0, 2.0], probs=[0.25, 0.75])
    d.scale(4.0)
    assert math.isclose(d.total_prob(), 4.0)
    assert math.isclose(list(d.probs)[0], 1.0)

    d.normalize()
    assert math.isclose(d.total_prob(), 1.0)
    assert math.isclose(list(d.probs)[0], 0.25)

    # The non-mutating variants leave the original alone.
    scaled = d.scaled(2.0)
    assert math.isclose(scaled.total_prob(), 2.0)
    assert math.isclose(d.total_prob(), 1.0)

    truncated = IsoSpecPy.IsoTotalProb(0.9, "C100H200N40O40")
    normalized = truncated.normalized()
    assert truncated.total_prob() < 1.0
    assert math.isclose(normalized.total_prob(), 1.0)


def test_mass_transforms():
    d = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])

    d.add_mass(10.0)
    assert list(d.masses) == [11.0, 12.0]

    d.mul_mass(2.0)
    assert list(d.masses) == [22.0, 24.0]

    # add_mul_mass: (m + add) * mul
    d2 = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    d2.add_mul_mass(1.0, 3.0)
    assert list(d2.masses) == [6.0, 9.0]

    # mul_add_mass: m * mul + add
    d3 = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    d3.mul_add_mass(3.0, 1.0)
    assert list(d3.masses) == [4.0, 7.0]


def test_resample_redistributes_a_fixed_current():
    d = IsoSpecPy.IsoTotalProb(0.999, "C100H200N40O40")
    d.normalize()
    d.resample(100000)
    assert all(float(p).is_integer() and p >= 0 for p in d.probs)
    assert math.isclose(sum(d.probs), 100000.0)
    d._recalculate_everything()
    assert math.isclose(d.total_prob(), 100000.0)


def test_sorting():
    d = IsoSpecPy.IsoTotalProb(0.999, "C6H12O6")
    before = set(peaks(d))

    d.sort_by_mass()
    assert list(d.masses) == sorted(d.masses)
    assert set(peaks(d)) == before

    d.sort_by_prob()
    assert list(d.probs) == sorted(d.probs)
    assert set(peaks(d)) == before

    # Sorting back by mass must not be short-circuited by the cached flag.
    d.sort_by_mass()
    assert list(d.masses) == sorted(d.masses)


def test_addition_concatenates():
    a = IsoDistribution(masses=[1.0, 2.0], probs=[0.3, 0.7])
    b = IsoDistribution(masses=[5.0], probs=[1.0])
    s = a + b
    assert len(s) == 3
    assert math.isclose(s.total_prob(), 2.0)
    assert peaks(s) == [(1.0, 0.3), (2.0, 0.7), (5.0, 1.0)]


def test_multiplication_convolves():
    a = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    b = IsoDistribution(masses=[10.0, 20.0], probs=[0.25, 0.75])
    p = a * b
    assert len(p) == 4
    assert math.isclose(p.total_prob(), 1.0)
    assert peaks(p) == sorted([(11.0, 0.125), (21.0, 0.375),
                               (12.0, 0.125), (22.0, 0.375)])


def test_convolution_of_elements_equals_the_molecule():
    carbon = IsoSpecPy.IsoThreshold(0.0, "C2")
    hydrogen = IsoSpecPy.IsoThreshold(0.0, "H6")
    direct = IsoSpecPy.IsoThreshold(0.0, "C2H6")

    product = carbon * hydrogen
    assert len(product) == len(direct)
    assert math.isclose(product.total_prob(), direct.total_prob(), rel_tol=1e-9)
    assert math.isclose(product.empiric_average_mass(),
                        direct.empiric_average_mass(), rel_tol=1e-12)


def test_linear_combination():
    a = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    b = IsoDistribution(masses=[3.0], probs=[1.0])
    mix = IsoDistribution.LinearCombination([a, b], [2.0, 10.0])
    assert len(mix) == 3
    assert math.isclose(mix.total_prob(), 12.0)
    assert peaks(mix) == [(1.0, 1.0), (2.0, 1.0), (3.0, 10.0)]


def test_binned_conserves_probability():
    d = IsoSpecPy.IsoTotalProb(0.999, "C50H100N20O20S5")
    total = d.total_prob()
    for width in (0.01, 0.1, 1.0, 10.0):
        binned = IsoSpecPy.IsoTotalProb(0.999, "C50H100N20O20S5").binned(width=width)
        binned._recalculate_everything()
        assert math.isclose(binned.total_prob(), total, rel_tol=1e-9)
        assert len(binned) <= len(d)
        assert list(binned.masses) == sorted(binned.masses)


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def test_empiric_statistics_by_hand():
    d = IsoDistribution(masses=[10.0, 20.0, 30.0], probs=[1.0, 2.0, 1.0])
    mean = (10.0 + 2 * 20.0 + 30.0) / 4.0
    variance = (100.0 + 0.0 + 100.0) / 4.0
    assert math.isclose(d.empiric_average_mass(), mean)
    assert math.isclose(d.empiric_variance(), variance)
    assert math.isclose(d.empiric_stddev(), math.sqrt(variance))


def test_empiric_statistics_match_the_theoretical_ones():
    for formula in ("C6H12O6", "C10H20O2", "S4"):
        iso = IsoSpecPy.Iso(formula)
        full = IsoSpecPy.IsoThreshold(0.0, formula)
        assert math.isclose(full.total_prob(), 1.0, rel_tol=1e-9)
        assert math.isclose(full.empiric_average_mass(),
                            iso.getTheoreticalAverageMass(), rel_tol=1e-9)
        assert math.isclose(full.empiric_variance(), iso.variance(), rel_tol=1e-6)


def test_numpy_views_share_the_data():
    np = pytest.importorskip("numpy")
    d = IsoSpecPy.IsoTotalProb(0.99, "C6H12O6")
    masses = d.np_masses()
    probs = d.np_probs()
    assert len(masses) == len(d)
    assert len(probs) == len(d)
    assert np.isclose(probs.sum(), d.total_prob())
    assert np.allclose(masses, list(d.masses))
