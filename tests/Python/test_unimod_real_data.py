"""Real-data process/benchmark tests for [UNIMOD:<id>]-aware peptide parsing.

Fixture: tests/Python/data/sample_peptide_sequences.txt -- 988 real peptide
sequences from this pipeline's own dump_peptides output, vendored here for
test coverage (see that file's header comment for provenance). A deterministic
subset is decorated with real Unimod ids (data/unimod.csv) to exercise the
mods-aware path against realistic sequences; the undecorated remainder is used
to compare the new general parser against the old fixed-CHNOSSe one it sits
alongside (parse_fasta_c / isoFromFasta, both unchanged by this feature).
"""

import csv
import math
import os
import random
import time

import IsoSpecPy
from IsoSpecPy.isoFFI import isoFFI
from IsoSpecPy.IsoSpecPy import _UNIMOD_CSV_PATH

HERE = os.path.dirname(__file__)
FIXTURE_PATH = os.path.join(HERE, "data", "sample_peptide_sequences.txt")
UNIMOD_CSV_PATH = os.path.join(HERE, "..", "..", "data", "unimod.csv")


def _load_sequences():
    with open(FIXTURE_PATH) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]


def _load_unimod_ids():
    with open(UNIMOD_CSV_PATH) as f:
        return [int(row["id"]) for row in csv.DictReader(f) if row["composition"]]


def _decorate(seq, mod_id, rng):
    """Insert one synthetic [UNIMOD:<mod_id>] at a random valid position."""
    kind = rng.choice(("internal", "nterm", "cterm"))
    if kind == "nterm":
        return "[UNIMOD:{}]-{}".format(mod_id, seq)
    if kind == "cterm":
        return "{}-[UNIMOD:{}]".format(seq, mod_id)
    idx = rng.randrange(len(seq))
    return seq[: idx + 1] + "[UNIMOD:{}]".format(mod_id) + seq[idx + 1 :]


def _build_mixed_fixture():
    """Real sequences, ~20% decorated with a real Unimod id at a random valid
    position. A randomly-chosen (sequence, mod) pairing can be chemically
    impossible (e.g. a mod whose delta removes a Sulfur applied to a
    Sulfur-free sequence) -- that's a correct rejection by Iso's existing
    non-negative-count validation, not a parser bug, but it's not what this
    fixture is for, so retry with the next candidate id (fixed, deterministic
    order) until one actually resolves."""
    sequences = _load_sequences()
    unimod_ids = _load_unimod_ids()
    rng = random.Random(0)
    plain, mixed = [], []
    for i, seq in enumerate(sequences):
        if i % 5 == 0:
            candidates = list(unimod_ids)
            rng.shuffle(candidates)
            for mod_id in candidates:
                decorated = _decorate(seq, mod_id, rng)
                try:
                    IsoSpecPy.Iso(peptide_sequence=decorated)
                except Exception:
                    continue
                mixed.append(decorated)
                break
            else:  # pragma: no cover -- should never happen with 979 candidate ids
                raise AssertionError("no compatible Unimod id found for " + seq)
        else:
            plain.append(seq)
            mixed.append(seq)
    return plain, mixed


# ---------------------------------------------------------------------------
# Old (fixed CHNOSSe, unchanged fasta.h/fasta.cpp) vs new (general, mods-aware)
# C ABI entry points, called directly -- bypasses any Python-level wrapper
# differences to compare the actual underlying implementations.
# ---------------------------------------------------------------------------


def _old_atom_counts(seq):
    buf = isoFFI.ffi.new("int[6]")
    isoFFI.clib.parse_fasta_c(seq.encode("ascii"), buf)
    return tuple(buf[i] for i in range(6))


def _new_atom_counts(seq):
    handle = isoFFI.clib.parseFastaWithModsC(seq.encode("ascii"), os.fsencode(_UNIMOD_CSV_PATH))
    assert handle != isoFFI.ffi.NULL, "parseFastaWithModsC failed on plain sequence: " + seq
    try:
        n = isoFFI.clib.compositionSizeC(handle)
        symbols = isoFFI.clib.compositionSymbolsC(handle)
        counts = isoFFI.clib.compositionCountsC(handle)
        return {
            isoFFI.ffi.string(symbols[i]).decode("ascii"): counts[i] for i in range(n)
        }
    finally:
        isoFFI.clib.deleteCompositionC(handle)


def _old_iso_masses(seq):
    # add_water=False: matches what the old Python Iso(fasta=...) actually did
    # (ParseFASTA -> parse_fasta_c, the bare no-water residue composition --
    # it never called this add_water=True C++ entry point at all). The new
    # Python Iso(peptide_sequence=...) is equally water-less (parseFastaWithModsC
    # is the same kind of bare-composition call), so this is the fair comparison.
    handle = isoFFI.clib.isoFromFasta(seq.encode("ascii"), False, False)
    assert handle != isoFFI.ffi.NULL, "isoFromFasta failed on plain sequence: " + seq
    try:
        return (
            isoFFI.clib.getMonoisotopicPeakMassIso(handle),
            isoFFI.clib.getTheoreticalAverageMassIso(handle),
        )
    finally:
        isoFFI.clib.deleteIso(handle)


def _new_iso_masses(seq):
    iso = IsoSpecPy.Iso(peptide_sequence=seq)
    return iso.getMonoisotopicPeakMass(), iso.getTheoreticalAverageMass()


CHNOSSE = ("C", "H", "N", "O", "S", "Se")


def test_process_all_real_sequences_including_unimod_decorated():
    """The mods-aware parser must handle every sequence in the real-data
    sample, plain and Unimod-decorated alike, with no exceptions."""
    _, mixed = _build_mixed_fixture()
    assert len(mixed) == 988

    failures = []
    for seq in mixed:
        try:
            mass = IsoSpecPy.Iso(peptide_sequence=seq).getMonoisotopicPeakMass()
        except Exception as e:  # noqa: BLE001 -- want to collect every failure, not stop at the first
            failures.append((seq, repr(e)))
            continue
        if not (math.isfinite(mass) and mass > 0):
            failures.append((seq, "non-finite or non-positive mass: {}".format(mass)))

    assert not failures, "{} / {} sequences failed:\n{}".format(
        len(failures), len(mixed), "\n".join("{}: {}".format(s, e) for s, e in failures[:20])
    )


def test_old_and_new_composition_match_exactly_on_plain_sequences():
    plain, _ = _build_mixed_fixture()
    mismatches = []
    for seq in plain:
        old = dict(zip(CHNOSSE, _old_atom_counts(seq)))
        old = {k: v for k, v in old.items() if v != 0}
        new = _new_atom_counts(seq)
        if old != new:
            mismatches.append((seq, old, new))
    assert not mismatches, "{} / {} plain sequences mismatched:\n{}".format(
        len(mismatches), len(plain), mismatches[:10]
    )


def test_old_and_new_iso_masses_match_exactly_on_plain_sequences():
    plain, _ = _build_mixed_fixture()
    mismatches = []
    for seq in plain:
        old_mono, old_avg = _old_iso_masses(seq)
        new_mono, new_avg = _new_iso_masses(seq)
        if old_mono != new_mono or old_avg != new_avg:
            mismatches.append((seq, (old_mono, old_avg), (new_mono, new_avg)))
    assert not mismatches, "{} / {} plain sequences' masses mismatched:\n{}".format(
        len(mismatches), len(plain), mismatches[:10]
    )


def test_old_vs_new_runtime_on_plain_sequences():
    """Not a strict perf gate (machine-dependent) -- reports real numbers and
    only fails on a gross regression (>10x slower per call on average)."""
    plain, _ = _build_mixed_fixture()
    reps = 3

    t0 = time.perf_counter()
    for _ in range(reps):
        for seq in plain:
            _old_iso_masses(seq)
    old_elapsed = time.perf_counter() - t0

    t0 = time.perf_counter()
    for _ in range(reps):
        for seq in plain:
            _new_iso_masses(seq)
    new_elapsed = time.perf_counter() - t0

    calls = reps * len(plain)
    old_mean_us = old_elapsed / calls * 1e6
    new_mean_us = new_elapsed / calls * 1e6
    ratio = new_elapsed / old_elapsed

    summary = (
        "old-vs-new runtime over {} plain sequences x {} reps: "
        "old={:.2f}us/call new={:.2f}us/call ratio(new/old)={:.2f}x".format(
            len(plain), reps, old_mean_us, new_mean_us, ratio
        )
    )
    print(summary)
    assert ratio < 10.0, "new path is >10x slower than old: " + summary
