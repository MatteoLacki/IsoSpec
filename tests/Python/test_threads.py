"""The library used from several Python threads at once.

cffi releases the GIL around every C call, so these threads really do run
the C++ concurrently -- on a free-threaded interpreter, everything else too.
Each test compares every thread's answers with a single-threaded reference.
The C++ suite's ThreadSanitizer build (tests/C++/test_threads.cpp) is what
proves the absence of data races; this proves the binding's path through the
C ABI gives the same answers under contention.
"""

import concurrent.futures

import IsoSpecPy
from IsoSpecPy.IsoSpecPy import ParsePeptideSequence

THREADS = 8
FORMULAS = ["C100H202", "C2H6O1", "C520H817N139O147S8", "Se2Fe3C40H60", "H2O1"]


def _run_everywhere(fn, jobs):
    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as pool:
        return list(pool.map(fn, jobs))


def _summary(formula):
    iso = IsoSpecPy.Iso(formula)
    threshold = IsoSpecPy.IsoThreshold(1e-4, formula=formula)
    total = IsoSpecPy.IsoTotalProb(0.999, formula=formula)
    binned = IsoSpecPy.IsoBinned(1.0, formula=formula, target_total_prob=0.999)
    return (
        iso.getModeLProb(),
        iso.getMonoisotopicPeakMass(),
        len(threshold),
        sum(threshold.masses),
        len(total),
        len(binned),
    )


def test_independent_formulas_from_many_threads():
    want = {f: _summary(f) for f in FORMULAS}
    jobs = FORMULAS * (4 * THREADS)
    for formula, got in zip(jobs, _run_everywhere(_summary, jobs)):
        assert got == want[formula]


def test_unimod_override_tables_from_many_threads(tmp_path):
    # Two override files giving the same id different compositions, so a
    # thread handed the other one's table gets a different answer.
    path_a = tmp_path / "a.csv"
    path_b = tmp_path / "b.csv"
    path_a.write_text("id,name,mono_mass,composition\n4,Carbamidomethyl,57.021464,H3C2N1O1\n")
    path_b.write_text("id,name,mono_mass,composition\n4,Heavier,58.0,H4C2N1O1\n")
    paths = [None, str(path_a), str(path_b)]
    sequence = "MPEPTC[UNIMOD:4]DEKC[UNIMOD:4]"

    def parse(path):
        return ParsePeptideSequence(sequence, unimod_db_path=path)

    want = {path: parse(path) for path in paths}
    assert want[str(path_a)] != want[str(path_b)]

    jobs = paths * (20 * THREADS)
    for path, got in zip(jobs, _run_everywhere(parse, jobs)):
        assert got == want[path]
