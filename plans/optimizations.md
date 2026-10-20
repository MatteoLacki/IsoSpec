# Stochastic sampling: where the time goes, and what could be optimised

## Context

Question asked: do `log`/`exp`/`pow` cost a lot in random sampling, could dedicated
implementations speed them up, and how much SIMD does the project use?

Short answer: the transcendental maths is **not** the bottleneck. `pow` is ~8% of a large
`FromStochastic` run; `log`/`exp` (BTRD) are ~0%. About 90% of the time is spent walking
configurations in the layered generator.

## Measurement

- Tree: `fcc06ed`.
- Workload: `FixedEnvelope::FromStochastic`, formula `C23832H37816N6528O7031S170`,
  3e7 molecules, precision 0.999 (the reference case from the `rdvariate_unif01` comment
  in `isoMath.h`), via `tests/C++/cli/from_formula_stochastic.cpp`.
- Machine: AMD Ryzen 7 4800H (Zen 2), g++, `-O3 -march=native -g -static`.
- Uninstrumented wall time: ~10.0 s (3 runs: 10.01, 9.83, 10.35 s).
- `perf` unavailable (`kernel.perf_event_paranoid=4`), so profiled with gprof on a
  **statically linked** `-pg` build (static so libm's `pow`/`exp`/`log` land in the
  sampled text). `gprof -l` crashes on this binary ("somebody miscounted"), so the raw
  `gmon.out` PC histogram was decoded and attributed through inlining with
  `addr2line -i` (script in the appendix).
- 1213 samples (~±1 pp per share). gprof's own `mcount` overhead (~10%) removed and shares
  renormalised.

## Results

| Where | Share |
|---|---|
| Layered-generator walk (total) | **~90%** |
| — skip loop in `IsoStochasticGeneratorTemplate::advanceToNextConfiguration` (`isoSpec++.h:739-744`, inlined) | ~41% |
| — `IsoLayeredGeneratorTemplate::carry()` (272M calls) | ~26% |
| — `IsoLayeredGeneratorTemplate::prob()` | ~8% |
| — `IsoLayeredGeneratorTemplate::recalc()` | ~7% |
| — `IsoLayeredGeneratorTemplate::advanceToNextConfiguration` | ~4% |
| — marginal lookups (`get_lProb`, `get_mass`, `get_prob`) | ~3% |
| `pow` in `rdvariate_beta_1_b` (`__ieee754_pow_fma` + `__pow`) | **~8%** |
| `std::mt19937` + rest of the beta draw | ~1% |
| `store_conf` | ~1% |
| BTRD `log`/`exp` (`btrd.h`) | **~0%** (6057 `btrd` calls in the whole run) |

- Peptide-sized formulas (`C60H100N17O18S1`, `C250H400N70O75S2`, 3e7 molecules) finish in
  microseconds; nothing to optimise in sampling there.
- The RNG engine itself is ~0.5%, so swapping `mt19937` for a faster engine is a
  hygiene/quality change, not a speed one.

## Options for the maths (Amdahl ceiling ~1.09x on this workload)

Gains below are estimates, not measured.

- **Beta draw without `pow`.** Replace `1 - pow(u, 1/b)` with `-expm1(log(u)/b)`. `b` is
  `to_sample_left`, usually huge, so the argument is tiny and a 2–3 term polynomial can
  replace `expm1`; one `log` per draw remains. Est. ~4% end-to-end. Side benefit: removes
  the cancellation in `1 - u^(1/b)` for large `b` (relative error ~1e-9 at b=3e7;
  harmless today).
- **Ziggurat Exp(1) variates** (NumPy's approach): `beta(1,b) = -expm1(-E/b)` with `E`
  from a ziggurat, so no transcendental call in the common case. Est. ~6–7%. More code.
- **SIMD `log`/`exp`**: only pays with batching. `E = -log(u)` does not depend on `b`, so
  a buffer of Exp(1) variates could be pre-generated in vectors — heavy machinery for ≤8%.

## SIMD in the project today

- Two kernels only: `fill_run` and `bin_run` in `src/IsoSpec++/isa_kernel_impl.h`, using
  `std::experimental::simd`. They vectorise the innermost-marginal run of
  `IsoThresholdGenerator`: mass add + probability multiply (`fill_run`), plus the bin-index
  `floor` (`bin_run`; the scatter into bins stays scalar).
- Callers: `FixedEnvelope::threshold_init` (threshold mode) and `Binned`
  (`fill_bins_full`).
- Compiled at four ISA levels (baseline, AVX, AVX2+FMA, AVX-512) with runtime dispatch in
  `isa_dispatch.cpp`.
- No SIMD elsewhere: not in the layered, ordered or stochastic generators, nor in marginal
  setup. No vectorised `exp`/`log` anywhere — none is needed per configuration, since
  marginal probabilities are precomputed and each configuration costs only multiplies.
- Compiler auto-vectorisation elsewhere: not checked.

## Bigger lever (untested)

The skip loop at `isoSpec++.h:739-744` sums probabilities of zero-count configurations one
at a time until `confs_prob > chasing_prob`. Over the innermost marginal's run this is a
prefix-sum search over contiguous probabilities times a loop-invariant partial product —
the same shape `fill_run` already vectorises for the threshold generator. Whether a batched
kernel pays off here has not been tried.

## Reproduce

```bash
cd tests/C++
g++ -std=c++20 -I../../src/IsoSpec++ -O3 -march=native -g -pg -static \
    ../../src/IsoSpec++/unity-build.cpp cli/from_formula_stochastic.cpp -o stoch_pg
./stoch_pg C23832H37816N6528O7031S170 30000000 false   # writes gmon.out
gprof -b -p ./stoch_pg gmon.out | head -20             # flat profile
python3 hist.py                                        # inlined attribution (appendix)
```

### Appendix: `hist.py` (decode `gmon.out` histogram, attribute via `addr2line -i`)

```python
import struct, subprocess, collections, sys
d = open('gmon.out','rb').read()
assert d[:4] == b'gmon'
off = 20
bins = []
while off < len(d):
    tag = d[off]; off += 1
    if tag == 0:
        lo, hi, n, rate = struct.unpack_from('<QQiI', d, off); off += 8+8+4+4+15+1
        cnt = struct.unpack_from('<%dH' % n, d, off); off += 2*n
        step = (hi - lo) / n
        bins += [(int(lo + i*step), c) for i, c in enumerate(cnt) if c]
    elif tag == 1: off += 8+8+4
    elif tag == 2:
        k = struct.unpack_from('<I', d, off)[0]; off += 4 + k*16
    else: sys.exit('bad tag %d' % tag)
total = sum(c for _, c in bins)
out = subprocess.run(['addr2line','-f','-i','-C','-a','-e','stoch_pg'] + [hex(a) for a, _ in bins],
                     capture_output=True, text=True).stdout.splitlines()
blocks, cur = [], None
for l in out:
    if l.startswith('0x'):
        cur = []; blocks.append(cur)
    else: cur.append(l)
inner = collections.Counter(); line = collections.Counter()
for (a, c), b in zip(bins, blocks):
    fn = b[0] if b else '?'
    src = b[1].split('/')[-1].split(' ')[0] if len(b) > 1 else '?'
    inner[fn[:90]] += c
    line[src] += c
print('total samples', total)
print('--- innermost (inlined) function')
for k, v in inner.most_common(18): print('%6.2f%%  %s' % (100*v/total, k))
print('--- source line')
for k, v in line.most_common(25): print('%6.2f%%  %s' % (100*v/total, k))
```
