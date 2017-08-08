from IsoSpecPy import isoFFI
from collections import defaultdict


number_of_isotopic_entries = isoFFI.clib.NUMBER_OF_ISOTOPIC_ENTRIES

symbol_to_masses = defaultdict(tuple)
symbol_to_probs  = defaultdict(tuple)

for i in xrange(number_of_isotopic_entries):
    symbol = isoFFI.ffi.string(isoFFI.clib.elem_table_symbol[i])
    symbol_to_masses[symbol] += (isoFFI.clib.elem_table_mass[i],)
    symbol_to_probs[symbol] += (isoFFI.clib.elem_table_probability[i],)

symbol_to_masses = dict(symbol_to_masses)
symbol_to_probs = dict(symbol_to_probs)

