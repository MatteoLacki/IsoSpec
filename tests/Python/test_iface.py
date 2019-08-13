import IsoSpecPy
from IsoSpecPy.Formulas import *

glu = IsoSpecPy.IsoThreshold(0.0, formula=glucose)
ca = IsoSpecPy.IsoThreshold(0.0, formula=caffeine)

print(ca.wassersteinDistance(glu))

wa = IsoSpecPy.IsoThreshold(0.0, formula=water)
ox = IsoSpecPy.IsoThreshold(0.0, formula=oxygen)
s = wa+ox

print(list(s))
print(len(s), len(wa), len(ox))

o = IsoSpecPy.IsoThreshold(0.0, formula="H1")
print(list(o*o))

sur = IsoSpecPy.IsoThreshold(0.0, formula=surcose)
(glu*glu).plot()
print((sur*wa).wassersteinDistance(glu*glu))


