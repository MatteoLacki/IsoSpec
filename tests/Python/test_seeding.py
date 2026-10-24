"""Reproducible sampling: the seed= arguments of IsoStochastic,
IsoStochasticGenerator and IsoDistribution.resample, and the deterministic
output order of IsoTotalProb."""

import concurrent.futures

import pytest

import IsoSpecPy

FORMULA = "C520H817N139O147S8"


def _stochastic(seed):
    d = IsoSpecPy.IsoStochastic(100000, formula=FORMULA, seed=seed)
    return list(d.masses), list(d.probs)


def test_iso_stochastic_same_seed_same_sample():
    assert _stochastic(7) == _stochastic(7)
    assert _stochastic(7) != _stochastic(8)


def test_seeded_sample_is_the_same_on_another_thread():
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        assert pool.submit(_stochastic, 7).result() == _stochastic(7)


def test_stochastic_generator_same_seed_same_sample():
    def draw(seed):
        return list(IsoSpecPy.IsoStochasticGenerator(50000, formula="C100H202", seed=seed))
    assert draw(11) == draw(11)
    assert draw(11) != draw(12)


def test_resample_same_seed_same_counts():
    def counts(seed):
        d = IsoSpecPy.IsoTotalProb(0.999, formula="C100H200N40O40")
        d.normalize()
        d.resample(100000, seed=seed)
        return list(d.probs)
    assert counts(3) == counts(3)
    assert counts(3) != counts(4)


def test_largest_seed_is_accepted():
    assert _stochastic(2**64 - 1) == _stochastic(2**64 - 1)


@pytest.mark.parametrize("seed, error", [(-1, ValueError), (2**64, ValueError), (1.5, TypeError),
                                         ("7", TypeError), (True, TypeError)])
def test_bad_seeds_are_rejected(seed, error):
    with pytest.raises(error):
        IsoSpecPy.IsoStochastic(1000, formula="C100H202", seed=seed)
    with pytest.raises(error):
        IsoSpecPy.IsoStochasticGenerator(1000, formula="C100H202", seed=seed)


def test_iso_total_prob_order_is_the_same_every_time():
    # The trimming quickselect used to draw its pivots from a randomly seeded
    # generator, so the order of the result differed between calls.
    for optimize in (True, False):
        def run():
            d = IsoSpecPy.IsoTotalProb(0.9999, formula="C200H400O20", get_minimal_pset=optimize)
            return list(d.masses), list(d.probs)
        assert run() == run()
