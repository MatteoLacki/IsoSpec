"""Spectrum distances through the Python layer.

The Wasserstein family is exposed as IsoDistribution methods; besides the
numerics (checked against closed forms) these tests pin the Python-side
contract: an unnormalized pair raises ValueError rather than returning the NaN
the C ABI uses to signal the error.
"""

import math

import pytest

import IsoSpecPy
from IsoSpecPy import IsoDistribution


def two_point(mass_a, mass_b, prob_a=0.5, prob_b=0.5):
    return IsoDistribution(masses=[mass_a, mass_b], probs=[prob_a, prob_b])


def test_distance_to_self_is_zero():
    for formula in ("H2O1", "C6H12O6", "C50H100N20O20S5"):
        a = IsoSpecPy.IsoThreshold(1e-8, formula, absolute=True)
        b = IsoSpecPy.IsoThreshold(1e-8, formula, absolute=True)
        a.normalize()
        b.normalize()
        assert math.isclose(a.wassersteinDistance(b), 0.0, abs_tol=1e-12)
        assert math.isclose(a.orientedWassersteinDistance(b), 0.0, abs_tol=1e-12)
        assert math.isclose(a.abyssalWassersteinDistance(b, 1.0), 0.0, abs_tol=1e-12)


def test_two_point_distance_is_the_transport_cost():
    for gap in (0.5, 1.0, 17.25):
        a = IsoDistribution(masses=[0.0], probs=[1.0])
        b = IsoDistribution(masses=[gap], probs=[1.0])
        assert math.isclose(a.wassersteinDistance(b), gap, rel_tol=1e-12)
        assert math.isclose(b.wassersteinDistance(a), gap, rel_tol=1e-12)

    # Half the mass moves two units: cost one.
    assert math.isclose(two_point(0.0, 2.0).wassersteinDistance(two_point(0.0, 4.0)),
                        1.0, rel_tol=1e-12)


def test_shifting_a_spectrum_puts_it_at_that_distance():
    a = IsoSpecPy.IsoThreshold(1e-10, "C6H12O6", absolute=True)
    a.normalize()
    for shift in (0.001, 0.5, 3.0, 1000.0):
        b = IsoSpecPy.IsoThreshold(1e-10, "C6H12O6", absolute=True)
        b.normalize()
        b.add_mass(shift)
        assert math.isclose(a.wassersteinDistance(b), shift, rel_tol=1e-9)
        assert math.isclose(abs(a.orientedWassersteinDistance(b)), shift, rel_tol=1e-9)


def test_oriented_distance_is_antisymmetric():
    a = two_point(0.0, 1.0)
    b = two_point(2.0, 3.0)
    assert math.isclose(a.orientedWassersteinDistance(b),
                        -b.orientedWassersteinDistance(a), rel_tol=1e-12)
    assert math.isclose(abs(a.orientedWassersteinDistance(b)),
                        a.wassersteinDistance(b), rel_tol=1e-12)


def test_distance_is_a_metric_on_real_spectra():
    a = IsoSpecPy.IsoThreshold(1e-8, "C6H12O6", absolute=True)
    b = IsoSpecPy.IsoThreshold(1e-8, "C6H12O5N1", absolute=True)
    c = IsoSpecPy.IsoThreshold(1e-8, "C5H12O6", absolute=True)
    for d in (a, b, c):
        d.normalize()

    ab = a.wassersteinDistance(b)
    bc = b.wassersteinDistance(c)
    ac = a.wassersteinDistance(c)
    assert ab > 0.0 and bc > 0.0 and ac > 0.0
    assert ac <= ab + bc + 1e-9
    assert math.isclose(ab, b.wassersteinDistance(a), rel_tol=1e-12)


def test_unnormalized_spectra_raise():
    a = two_point(1.0, 2.0)
    b = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 5.0])
    with pytest.raises(ValueError):
        a.wassersteinDistance(b)
    with pytest.raises(ValueError):
        a.orientedWassersteinDistance(b)


def test_abyssal_distance_transports_or_discards():
    a = IsoDistribution(masses=[0.0], probs=[1.0])
    b = IsoDistribution(masses=[1.0], probs=[1.0])

    # A deep abyss makes transport the cheaper option: cost is the gap.
    assert math.isclose(a.abyssalWassersteinDistance(b, 10.0), 1.0, rel_tol=1e-12)

    # A shallow one makes discarding cheaper: cost is (mass discarded)*depth/2.
    depth = 0.25
    assert math.isclose(a.abyssalWassersteinDistance(b, depth), 2 * depth * 0.5,
                        rel_tol=1e-12)

    # Never above the plain Wasserstein distance, and monotone in the depth.
    x = IsoSpecPy.IsoThreshold(1e-8, "C6H12O6", absolute=True)
    y = IsoSpecPy.IsoThreshold(1e-8, "C6H12O5N1", absolute=True)
    x.normalize()
    y.normalize()
    plain = x.wassersteinDistance(y)
    previous = -1.0
    for depth in (0.1, 0.5, 1.0, 2.0, 5.0):
        value = x.abyssalWassersteinDistance(y, depth)
        assert value >= previous - 1e-9
        assert value <= plain + 1e-9
        previous = value


def test_abyssal_distance_other_scale():
    a = IsoDistribution(masses=[0.0], probs=[1.0])
    b = IsoDistribution(masses=[0.0], probs=[0.5])
    # Scaling b by two makes the spectra identical.
    assert math.isclose(a.abyssalWassersteinDistance(b, 10.0, 2.0), 0.0, abs_tol=1e-12)
    # Unscaled, half a unit of probability has nowhere to go.
    assert math.isclose(a.abyssalWassersteinDistance(b, 10.0, 1.0), 0.5 * 10.0 * 0.5,
                        rel_tol=1e-12)


def test_wasserstein_match():
    a = IsoSpecPy.IsoThreshold(1e-6, "C6H12O6", absolute=True)
    b = IsoSpecPy.IsoThreshold(1e-6, "C6H12O6", absolute=True)
    a.normalize()
    b.normalize()

    unmatched_a, unmatched_b, flow = a.wassersteinMatch(b, 0.5)
    assert math.isclose(flow, 1.0, rel_tol=1e-9)
    assert math.isclose(unmatched_a, 0.0, abs_tol=1e-9)
    assert math.isclose(unmatched_b, 0.0, abs_tol=1e-9)

    # Pushed far apart, nothing matches — and nothing reads out of bounds.
    far = IsoSpecPy.IsoThreshold(1e-6, "C6H12O6", absolute=True)
    far.normalize()
    far.add_mass(1000.0)
    unmatched_a, unmatched_far, flow = a.wassersteinMatch(far, 0.5)
    assert math.isclose(flow, 0.0, abs_tol=1e-9)
    assert math.isclose(unmatched_a, 1.0, rel_tol=1e-9)
    assert math.isclose(unmatched_far, 1.0, rel_tol=1e-9)

    # Matched plus unmatched accounts for all the probability.
    other = IsoSpecPy.IsoThreshold(1e-6, "C6H12O5N1", absolute=True)
    other.normalize()
    for flow_dist in (0.0, 0.01, 0.1, 1.0):
        unmatched_a, unmatched_other, flow = a.wassersteinMatch(other, flow_dist)
        assert flow >= 0.0
        assert math.isclose(flow + unmatched_a, a.total_prob(), rel_tol=1e-9)
        assert math.isclose(flow + unmatched_other, other.total_prob(), rel_tol=1e-9)


def test_wasserstein_match_flow_distance_matters():
    # 1.0 pairs with 0.9 only once the flow distance exceeds the 0.1 gap.
    a = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    b = IsoDistribution(masses=[0.9, 2.01], probs=[0.4, 0.6])
    unmatched_a, unmatched_b, flow = a.wassersteinMatch(b, 0.05)
    assert math.isclose(flow, 0.5, rel_tol=1e-12)
    assert math.isclose(unmatched_a, 0.5, rel_tol=1e-12)
    assert math.isclose(unmatched_b, 0.5, rel_tol=1e-12)

    a2 = IsoDistribution(masses=[1.0, 2.0], probs=[0.5, 0.5])
    b2 = IsoDistribution(masses=[0.9, 2.01], probs=[0.4, 0.6])
    _, _, wider_flow = a2.wassersteinMatch(b2, 0.2)
    assert wider_flow > flow
