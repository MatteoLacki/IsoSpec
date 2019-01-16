import cppyy
from math import exp


cppyy.include("../../IsoSpec++/isoSpec++.h")
cppyy.load_library("../../IsoSpec++/libIsoSpec++.so")



IsoSpec = cppyy.gbl.IsoSpec
std = cppyy.gbl.std

i = IsoSpec.Iso("H2O1")
it = IsoSpec.IsoThresholdGenerator(std.move(i), 0.0001)

while it.advanceToNextConfiguration():
    print it.mass(), exp(it.lprob())
