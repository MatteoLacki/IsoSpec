import IsoSpecPy
from math import exp
from IsoSpecPy.Formulas import *


test_on = bovine_insulin
test_prob = 0.9999
i = IsoSpecPy.Iso(test_on)
print map(exp,i.getMarginalLogSizeEstimates(test_prob))

v = IsoSpecPy.IsoTotalProb(formula=test_on, prob_to_cover = test_prob, get_confs = True, get_minimal_pset = True)
print len(v)

acc = [set() for _ in xrange(v.dimNumber)]

for conf in v.confs:
    for i in xrange(v.dimNumber):
        acc[i].add(conf[i])

print map(len, acc)

