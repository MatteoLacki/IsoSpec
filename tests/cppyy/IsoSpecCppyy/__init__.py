import cppyy
import glob
import re

for header in glob.glob("../../../IsoSpec++/*.h"):
    if not "mman.h" in header:
        cppyy.include(header)
cppyy.load_library("../../../IsoSpec++/libIsoSpec++.so")


regex_pattern = re.compile('([A-Z][a-z]?)([0-9]*)')

def IsoParamsFromFormula(formula):
    global regex_pattern

    symbols = []
    atomCounts = []
    for elem, cnt in re.findall(regex_pattern, formula):
        symbols.append(elem)
        atomCounts.append(int(cnt) if cnt is not '' else 1)
    try:
        masses = tuple(PeriodicTbl.symbol_to_masses[s] for s in symbols)
        probs  = tuple(PeriodicTbl.symbol_to_probs[s]  for s in symbols)
    except KeyError:
        raise ValueError("Invalid formula")

    return (atomCounts, masses, probs)


