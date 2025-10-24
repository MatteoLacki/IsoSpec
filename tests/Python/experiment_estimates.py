import IsoSpecPy
from math import exp
from IsoSpecPy.Formulas import *
from IsoSpecPy.approximations import approximate_subisotopologues

test_on = horse_myoglobin
test_prob = 0.9999

print("Formula:", test_on, "Probability:", test_prob)

i = IsoSpecPy.Iso(test_on)
print ("From C++ code:", list(map(exp,i.getMarginalLogSizeEstimates(test_prob))))

symbols, _ = IsoSpecPy.ParseFormula(test_on)
dct = approximate_subisotopologues(test_on, test_prob)
print ("From Python:  ", [dct[s] for s in symbols])

v = IsoSpecPy.IsoTotalProb(formula=test_on, prob_to_cover = test_prob, get_confs = True, get_minimal_pset = True)

acc = [set() for _ in range(v.dimNumber)]

for conf in v.confs:
    for i in range(v.dimNumber):
        acc[i].add(conf[i])

print ("Real:", list(map(len, acc)))
print (len(v), "total confs.")
