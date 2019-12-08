import cppyy # NOTE: No import IsoSpecPy or anything like that!
#cppyy.include("unity-build.cpp")
cppyy.include("isoSpec++.h")
cppyy.load_library("libIsoSpec++.so")
i = cppyy.gbl.IsoSpec.Iso("C100H100")
iso = cppyy.gbl.IsoSpec.IsoThresholdGenerator(cppyy.gbl.std.move(i), 0.01)
while iso.advanceToNextConfiguration():
    print iso.prob(), iso.mass()
