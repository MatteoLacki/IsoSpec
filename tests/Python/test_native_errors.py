"""Native error sentinels must become exceptions, never invalid dereferences."""

import subprocess
import sys

import pytest
import IsoSpecPy as I
from IsoSpecPy.isoFFI import isoFFI


class FailingLibrary:
    """Inject a C ABI allocation failure while retaining the real cleanup calls."""

    def __init__(self, library, failing_name):
        self.library = library
        self.failing_name = failing_name

    def __getattr__(self, name):
        if name == self.failing_name:
            return lambda *args: isoFFI.ffi.NULL
        return getattr(self.library, name)


def test_oversized_formula_raises_instead_of_crashing():
    # Keep a regression from killing pytest itself.
    proc = subprocess.run([sys.executable, "-c", """
import IsoSpecPy as I
try:
    I.IsoTotalProb(.99, formula='C2147483647')
except ValueError:
    pass
else:
    raise AssertionError('expected a construction error')
"""], capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr


@pytest.mark.parametrize("function,build", [
    ("setupIso", lambda: I.Iso("C2")),
    ("setupIsoThresholdGenerator", lambda: I.IsoThresholdGenerator(.01, "C2")),
    ("setupIsoLayeredGenerator", lambda: I.IsoLayeredGenerator("C2")),
    ("setupIsoOrderedGenerator", lambda: I.IsoOrderedGenerator("C2")),
    ("setupIsoStochasticGenerator", lambda: I.IsoStochasticGenerator(10, "C2")),
    ("setupThresholdFixedEnvelope", lambda: I.IsoThreshold(.01, "C2")),
    ("setupTotalProbFixedEnvelope", lambda: I.IsoTotalProb(.99, "C2")),
    ("setupStochasticFixedEnvelope", lambda: I.IsoStochastic(10, "C2")),
    ("setupBinnedFixedEnvelope", lambda: I.IsoBinned(1., "C2")),
    ("getMarginalLogSizeEstimates", lambda: I.Iso("C2").getMarginalLogSizeEstimates(.99)),
])
def test_null_native_constructors(monkeypatch, function, build):
    monkeypatch.setattr(isoFFI, "clib", FailingLibrary(isoFFI.clib, function))
    with pytest.raises(ValueError, match="failed"):
        build()


@pytest.mark.parametrize("function,operation", [
    ("setupFixedEnvelope", lambda d: d.total_prob()),
    ("setupFixedEnvelopeWithConfs", lambda d: d.total_prob()),
    ("copyFixedEnvelope", lambda d: d.copy()),
    ("addEnvelopes", lambda d: d + d),
    ("convolveEnvelopes", lambda d: d * d),
    ("binnedEnvelope", lambda d: d.binned()),
    ("linearCombination", lambda d: I.IsoDistribution.LinearCombination([d], [1.])),
])
def test_null_envelope_operations(monkeypatch, function, operation):
    d = I.IsoTotalProb(.99, "C2", get_confs=function == "setupFixedEnvelopeWithConfs")
    monkeypatch.setattr(isoFFI, "clib", FailingLibrary(isoFFI.clib, function))
    with pytest.raises(ValueError, match="failed"):
        operation(d)


def test_second_wrapper_failure_releases_first(monkeypatch):
    a = I.IsoDistribution(masses=[1.], probs=[1.])
    b = I.IsoDistribution(masses=[2.], probs=[1.])
    library = isoFFI.clib
    released = []

    class FailSecondWrapper:
        calls = 0

        def __getattr__(self, name):
            return getattr(library, name)

        def setupFixedEnvelope(self, *args):
            self.calls += 1
            return library.setupFixedEnvelope(*args) if self.calls == 1 else isoFFI.ffi.NULL

        def deleteFixedEnvelope(self, handle, release):
            released.append(release)
            library.deleteFixedEnvelope(handle, release)

    monkeypatch.setattr(isoFFI, "clib", FailSecondWrapper())
    with pytest.raises(ValueError):
        a + b
    assert released == [True]
    assert list(a) == [(1., 1.)]
    assert list(b) == [(2., 1.)]


def test_empty_envelopes_still_allow_null_data_buffers():
    d = I.IsoTotalProb(0., "C2", get_confs=True)
    assert list(d) == []
    assert list(d.copy()) == []
