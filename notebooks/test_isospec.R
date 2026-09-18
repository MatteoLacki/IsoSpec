# R equivalent of notebooks/test_isospec.ipynb (Python):
#
#   import IsoSpecPy
#   env = IsoSpecPy.IsoTotalProb(prob_to_cover=0.999, peptide_sequence="AB[UNIMOD:42]C")
#   masses = env.np_masses()
#   probs = env.np_probs()
#
# R has no direct peptide_sequence= convenience on IsoSpecify (it only ever
# took a raw named-integer-vector molecule -- see git/isospec/plans/
# unimod_modification_parsing.md's R section) -- so this is a two-step
# equivalent: RParsePeptideSequence() to turn the [UNIMOD:<id>]-annotated
# sequence into a composition, then IsoSpecify() to compute the envelope
# from that composition. algo=0 (IsoSpecify's default) is the same
# total-probability/layered algorithm Python's IsoTotalProb wraps, so
# stopCondition here means the same thing prob_to_cover does there.

library(IsoSpecR)

composition <- RParsePeptideSequence("AB[UNIMOD:42]C")
print(composition)
#  H  C  N  O  S
# 22 14  2  3  3

res <- IsoSpecify(molecule = composition, stopCondition = 0.999)

masses <- res[, "mass"]
probs  <- res[, "prob"]

print(masses)
print(probs)

# Sanity check against the Python notebook's own printed output: the two
# should agree to float precision (verified this session -- first six rows
# matched Python's masses/probs to displayed precision, e.g. 362.0793 /
# 0.723716126 both sides).
