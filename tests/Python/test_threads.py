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


def test_generators_iterated_on_many_threads():
    def drain(kind):
        if kind == "threshold":
            gen = IsoSpecPy.IsoThresholdGenerator(1e-4, formula="C520H817N139O147S8")
        elif kind == "layered":
            gen = IsoSpecPy.IsoLayeredGenerator(formula="C200H300N50O60S2")
        else:
            gen = IsoSpecPy.IsoOrderedGenerator(formula="C100H202")
        out = []
        for mass, prob in gen:
            out.append((mass, prob))
            if kind != "threshold" and len(out) == 2000:
                break
        return kind, out

    kinds = ["threshold", "layered", "ordered"]
    want = dict(drain(k) for k in kinds)
    for kind, got in _run_everywhere(drain, kinds * (2 * THREADS)):
        assert got == want[kind]


def test_stochastic_generator_made_on_one_thread_used_on_another():
    # The generator must draw from whichever thread is sampling, not hold on
    # to the random-number state of the (pool) thread that built it.
    molecules = 100000

    def make(_):
        return IsoSpecPy.IsoStochasticGenerator(molecules, formula="C520H817N139O147S8")

    def sample_while_others_build(gen):
        return sum(count for _, count in gen)

    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as pool:
        generators = list(pool.map(make, range(2 * THREADS)))
        # Builders keep drawing for their own samples while the generators
        # they made are consumed elsewhere.
        busy = [pool.submit(lambda: len(IsoSpecPy.IsoStochastic(1000, formula="C100H202")))
                for _ in range(THREADS)]
        totals = list(pool.map(sample_while_others_build, generators))
        for b in busy:
            b.result()
    assert totals == [molecules] * len(generators)


def test_shared_objects_read_by_many_threads():
    iso = IsoSpecPy.Iso("C520H817N139O147S8")
    a = IsoSpecPy.IsoThreshold(1e-5, formula="C200H300N50O60S2")
    b = IsoSpecPy.IsoThreshold(1e-5, formula="H2O1")

    def read(_):
        # Read-only on the shared objects; anything that sorts or normalises
        # in place works on a private copy.
        mine = a.copy()
        mine.normalize()
        theirs = b.copy()
        theirs.normalize()
        return (
            iso.getModeLProb(),
            iso.getMonoisotopicPeakMass(),
            iso.getTheoreticalAverageMass(),
            len(a + b),
            len(a * b),
            mine.wassersteinDistance(theirs),
        )

    want = read(None)
    for got in _run_everywhere(read, range(4 * THREADS)):
        assert got == want


def test_errors_raised_on_many_threads():
    def fail(job):
        try:
            if job == "repeat":
                IsoSpecPy.Iso("C1C1")
            elif job == "symbol":
                IsoSpecPy.Iso("Xx5")
            elif job == "unimod":
                ParsePeptideSequence("PEP[UNIMOD:999999]TIDE")
            else:
                ParsePeptideSequence("PEPC[UNIMOD:4]", unimod_db_path="/nonexistent/unimod.csv")
        except ValueError:
            return True
        return False

    jobs = ["repeat", "symbol", "unimod", "missing"] * (5 * THREADS)
    assert all(_run_everywhere(fail, jobs))
    assert ParsePeptideSequence("PEPC[UNIMOD:4]")  # still works afterwards
