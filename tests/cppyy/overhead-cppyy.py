from tqdm import tqdm

import cppyy
import glob

for header in glob.glob("../../IsoSpec++/*.h"):
    if not "mman.h" in header:
        cppyy.include(header)
cppyy.load_library("../../IsoSpec++/libIsoSpec++.so")


t = 0.0
for x in tqdm(xrange(1000000)):
    i = cppyy.gbl.IsoSpec.Iso("C100H100N100O100")
    t += i.getTheoreticalAverageMass()
print t
