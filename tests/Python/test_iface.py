import IsoSpecPy
from IsoSpecPy.Formulas import *

glu = IsoSpecPy.IsoTotalProb(0.99, formula=glucose)
ca = IsoSpecPy.IsoTotalProb(0.99, formula=caffeine)

print ca.wassersteinDistance(glu)

