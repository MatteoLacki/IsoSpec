from __future__ import print_function
import IsoSpecPy
from IsoSpecPy.Formulas import *
import math


try:
    math.isclose
except AttributeError:
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
    math.isclose = isclose

glu = IsoSpecPy.IsoThreshold(0.0, formula=glucose)
ca = IsoSpecPy.IsoThreshold(0.0, formula=caffeine)

print("Checking Wasserstein distance...", end=' ')
print(ca.wassersteinDistance(glu), end=' ')
assert(math.isclose(ca.wassersteinDistance(glu), 14.03495145836358))
print("OK!")

print("Checking normalization... ", end='')

ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
print(ubiq.total_prob(), end=' ')
assert(math.isclose(ubiq.total_prob(), 0.9999, rel_tol=0.01))
ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
ubiq.scale(0.5)
assert(math.isclose(ubiq.total_prob(), 0.9999*0.5, rel_tol=0.01))
ubiq._recalculate_everything()
assert(math.isclose(ubiq.total_prob(), 0.9999*0.5, rel_tol=0.01))
ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
ubiq.scale(0.5)
ubiq.normalize()
assert(math.isclose(ubiq.total_prob(), 1.0))
ubiq._recalculate_everything()
assert(math.isclose(ubiq.total_prob(), 1.0))
print("OK!")


print("Checking addition...", end=' ')
wa = IsoSpecPy.IsoThreshold(0.0, formula=water)
ox = IsoSpecPy.IsoThreshold(0.0, formula=oxygen)
s = wa+ox
assert(math.isclose(s.total_prob(), 2.0))
assert(len(list(s)) == len(list(wa)) + len(list(ox)) == 15)
print("OK!")


print("Checking sorting...", end=' ')
ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
ubiq.sort_by_mass()
assert(list(ubiq.masses) == sorted(ubiq.masses))
ubiq.sort_by_prob()
assert(list(ubiq.probs) == sorted(ubiq.probs))
print("OK!")

print("Checking binning...", end= ' ')
ubiq = IsoSpecPy.IsoTotalProb(0.999999, ubiquitin)
#ubiq.plot()
print(len(ubiq), end = ' -> ')
bu = ubiq.binned()
ubiq._recalculate_everything()
print(len(bu), end = ' ')
bu._recalculate_everything()
assert(math.isclose(ubiq.total_prob(), bu.total_prob()))
print("OK!")

print("Checking convolution...", end=' ')
o = IsoSpecPy.IsoThreshold(0.0, formula="H1")
sur = IsoSpecPy.IsoThreshold(0.0, formula=surcose)
#(glu*glu).plot()
#(sur*wa).plot()
WSD = (sur*wa).wassersteinDistance(glu*glu)
print(WSD, end=' ')
assert math.isclose(WSD, 0.0, abs_tol=1e-7)
print("OK!")


print("Checking negative formulas... ", end="")
try:
    I = Iso(formula="C-10")
    print("FAIL: exception not thrown")
except Exception as e:
    print("""exception successfully obtained, message: "{}" -> OK!""".format((str(e))))


print("Checking FASTA + negative formulas... ", end="")
try:
    I = Iso(fasta = "C", formula="C-5")
    print("FAIL: exception not thrown")
except Exception as e:
    print("""exception successfully obtained, message: "{}" -> OK!""".format((str(e))))


print("Checking FASTA + modification... ", end="")
# Selenation of methionine
I = IsoSpecPy.IsoTotalProb(0.999, formula = "C5H9N1O1Se1")
I2 = IsoSpecPy.IsoTotalProb(0.999, fasta = "M", formula = "S-1Se1")
WSD = I.wassersteinDistance(I2)
print(WSD, end="")
assert(math.isclose(I.wassersteinDistance(I2), 0.0))
print(" -> OK!")


print("Checking empiric avg mass... ", end="")
formulas = "C1 P1 H1 H100 P100 N100 O100 C100H100N100 C100O100".split()
isos = [IsoSpecPy.Iso(formula) for formula in formulas]
dists = [IsoSpecPy.IsoThreshold(0.0, formula) for formula in formulas]
diffs = [abs(i.getTheoreticalAverageMass() - d.empiric_average_mass()) for i, d in zip(isos, dists)]
print(max(diffs), end="")
assert max(diffs) < 1e-6
print(" -> OK!")


print("Checking empiric variance... ", end="")
diffs = [abs(i.variance() - d.empiric_variance()) for i, d in zip(isos, dists)]
print(max(diffs), end="")
assert max(diffs) < 1e-6
print(" -> OK!")

print("Checking empiric stddev... ", end="")
diffs = [abs(i.stddev() - d.empiric_stddev()) for i, d in zip(isos, dists)]
print(max(diffs), end="")
assert max(diffs) < 1e-6
print(" -> OK!")
