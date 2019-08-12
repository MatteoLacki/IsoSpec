import IsoSpecPy
from IsoSpecPy.Formulas import *

glu = IsoSpecPy.IsoThreshold(0.0, formula=glucose)
ca = IsoSpecPy.IsoThreshold(0.0, formula=caffeine)

print(ca.wassersteinDistance(glu))

wa = IsoSpecPy.IsoThreshold(0.0, formula=water)
ox = IsoSpecPy.IsoThreshold(0.0, formula=oxygen)
s = wa+ox

s.plot()
print(list(s))
print(len(s), len(wa), len(ox))


