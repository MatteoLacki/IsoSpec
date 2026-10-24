"""Scalar arguments that would hang the library or reach undefined behaviour
are rejected with ValueError, including where the C ABI call returns void."""
import math

import pytest

import IsoSpecPy

NAN = float("nan")
INF = float("inf")
F = "C10H20O5S1"


def env():
    return IsoSpecPy.IsoTotalProb(0.999, F)


def test_total_prob_nan_target_raises():
    with pytest.raises(ValueError):
        IsoSpecPy.IsoTotalProb(NAN, F)


def test_binned_nan_target_raises():
    with pytest.raises(ValueError):
        IsoSpecPy.IsoBinned(1.0, formula=F, target_total_prob=NAN)


def test_threshold_nan_raises():
    with pytest.raises(ValueError):
        IsoSpecPy.IsoThreshold(NAN, formula=F)


@pytest.mark.parametrize("precision", [NAN, -1.0, 0.0, 2.0, INF])
def test_stochastic_precision_raises(precision):
    with pytest.raises(ValueError):
        IsoSpecPy.IsoStochastic(1000, F, precision=precision)


@pytest.mark.parametrize("factor", [NAN, INF, -1.0])
def test_scale_raises_and_leaves_envelope_unchanged(factor):
    e = env()
    probs = list(e.probs)
    with pytest.raises(ValueError):
        e.scale(factor)
    assert list(e.probs) == probs


@pytest.mark.parametrize("beta_bias", [NAN, -1.0])
def test_resample_beta_bias_raises(beta_bias):
    with pytest.raises(ValueError):
        env().resample(1000, beta_bias)


@pytest.mark.parametrize("method", ["add_mass", "mul_mass"])
@pytest.mark.parametrize("value", [NAN, INF])
def test_mass_shift_raises(method, value):
    with pytest.raises(ValueError):
        getattr(env(), method)(value)


@pytest.mark.parametrize("value", [NAN, -1.0])
def test_wasserstein_parameters_raise(value):
    a, b = env(), IsoSpecPy.IsoTotalProb(0.999, "C10H22O5S1")
    with pytest.raises(ValueError):
        a.abyssalWassersteinDistance(b, value)
    with pytest.raises(ValueError):
        a.wassersteinMatch(b, value)
    with pytest.raises(ValueError):
        a.abyssalWassersteinDistance(b, 1.0, value)


def test_linear_combination_bad_intensity_raises():
    with pytest.raises(ValueError):
        IsoSpecPy.IsoDistribution.LinearCombination([env()], [NAN])


def test_binned_bad_width_raises():
    with pytest.raises(ValueError):
        env().binned(-1.0)


def test_infinite_depth_and_flow_still_accepted():
    a, b = env(), IsoSpecPy.IsoTotalProb(0.999, "C10H22O5S1")
    assert not math.isnan(a.abyssalWassersteinDistance(b, INF))
    assert a.wassersteinMatch(b, INF)[1] >= 0


def test_normalize_zero_total_raises_and_empty_is_noop():
    z = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs=[0.0, 0.0])
    with pytest.raises(ValueError):
        z.normalize()
    e = IsoSpecPy.IsoTotalProb(0.0, F)
    e.normalize()
    assert len(e.masses) == 0
