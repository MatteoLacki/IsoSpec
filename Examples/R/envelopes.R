library(IsoSpecR)

bovine_insulin = c(C=254, H=377, N=65, O=75, S=6)

# Properties of the distribution, computed without enumerating it: the
# monoisotopic peak is not the most probable one for a molecule this size.
props = IsoProperties(bovine_insulin)
print(props$monoisotopicMass)
print(props$modeMass)
print(props$averageMass)

# The distribution at unit resolution: isotopologues summed into bins 1 Da
# wide, centred on whole daltons, until 99.99% of the probability is covered.
# Cheaper than IsoSpecify() followed by binning, since the individual
# isotopologues are never stored.
BINNED = IsoBinned(bovine_insulin, binWidth = 1)
print(BINNED[order(BINNED[, "mass"]), ])

# A simulated measurement of 10000 ions. With a seed the counts are the same
# on every run and every platform (and the same as IsoSpecPy gives).
MEASURED = IsoStochastic(bovine_insulin, ionicCurrent = 10000, seed = 1)
print(sum(MEASURED[, "prob"]))

# Envelopes are plain matrices with columns "mass" and "prob": normalising,
# shifting or sorting them is ordinary R.
normalised = function(envelope) {
  envelope[, "prob"] = envelope[, "prob"] / sum(envelope[, "prob"])
  envelope
}
THEORY = normalised(IsoSpecify(bovine_insulin, .999))
MEASURED = normalised(MEASURED)

# How far the measurement is from the theory, in daltons of moved probability.
print(wassersteinDistance(THEORY, MEASURED))

# A dimer, by convolution, and a 2:1 mixture with the monomer.
DIMER = envelopeConvolve(THEORY, THEORY)
MIXTURE = envelopeLinearCombination(list(THEORY, DIMER), c(2, 1))
print(envelopeBin(MIXTURE, width = 1))
