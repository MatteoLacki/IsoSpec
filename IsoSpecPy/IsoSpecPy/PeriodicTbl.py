from .isoFFI import isoFFI
from collections import defaultdict

try:
    xrange
except NameError:
    xrange = range

number_of_isotopic_entries = isoFFI.clib.NUMBER_OF_ISOTOPIC_ENTRIES

symbol_to_masses = defaultdict(tuple)
symbol_to_probs  = defaultdict(tuple)
symbol_to_atomic_number = {}

for i in xrange(number_of_isotopic_entries):
    symbol = isoFFI.ffi.string(isoFFI.clib.elem_table_symbol[i]).decode("ascii")
    symbol_to_masses[symbol] += (isoFFI.clib.elem_table_mass[i],)
    symbol_to_probs[symbol] += (isoFFI.clib.elem_table_probability[i],)
    symbol_to_atomic_number[symbol] = isoFFI.clib.elem_table_atomicNo[i]

symbol_to_masses = dict(symbol_to_masses)
symbol_to_probs = dict(symbol_to_probs)

# Several derivative convenience dicts...
symbol_to_massprob = dict((key, [zip(symbol_to_masses[key], symbol_to_probs[key])]) for key in symbol_to_probs.keys())

def crossprod(l1, l2):
    return sum(x1*x2 for x1, x2 in zip(l1, l2))

symbol_to_avg_mass = dict((key, crossprod(symbol_to_masses[key], symbol_to_probs[key])) for key in symbol_to_probs.keys())

def maxprod(l1, l2):
    return max(zip(l1, l2), key = lambda x: x[1])[0]

symbol_to_monoisotopic_mass = dict((key, maxprod(symbol_to_masses[key], symbol_to_probs[key])) for key in symbol_to_probs.keys())

