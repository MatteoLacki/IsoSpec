"""Coverage for the binned envelope path (FixedEnvelope::Binned).

Binning redistributes peaks into fixed-width mass bins.  It must:
  * preserve total probability (mass moves, probability does not),
  * never produce more points than the unbinned envelope,
  * return mass-sorted output.

This exercises the path that the SIMD binned work will touch, across several
molecules, target probabilities and bin widths.
"""

import math

import pytest

import IsoSpecPy

MOLECULES = ["H2O1", "C6H12O6", "C50H100N20O20S5", "C100H200N40O40",
             "C520H817N139O147S8", "Fe2O3", "C1000"]
TARGETS = [0.5, 0.9, 0.99, 0.999]
WIDTHS = [0.01, 0.1, 1.0]


@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("molecule", MOLECULES)
def test_binned_preserves_total_prob(molecule, target, width):
    unbinned_tp = IsoSpecPy.IsoTotalProb(target, molecule).total_prob()

    binned = IsoSpecPy.IsoTotalProb(target, molecule).binned(width=width)
    binned._recalculate_everything()

    assert math.isclose(binned.total_prob(), unbinned_tp, rel_tol=1e-9)


@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("molecule", MOLECULES)
def test_binned_not_larger_and_sorted(molecule, width):
    env = IsoSpecPy.IsoTotalProb(0.999, molecule)
    n_unbinned = len(env)

    binned = IsoSpecPy.IsoTotalProb(0.999, molecule).binned(width=width)
    binned._recalculate_everything()

    assert len(binned) <= n_unbinned
    assert len(binned) > 0
    masses = list(binned.masses)
    assert masses == sorted(masses)


def test_binned_full_distribution():
    # target >= 1 routes through the threshold (full-enumeration) filler, a
    # different code path from the layered one exercised above.
    full = IsoSpecPy.IsoTotalProb(1.0, "C6H12O6")
    tp = full.total_prob()
    binned = IsoSpecPy.IsoTotalProb(1.0, "C6H12O6").binned(width=0.5)
    binned._recalculate_everything()
    assert math.isclose(binned.total_prob(), tp, rel_tol=1e-9)
