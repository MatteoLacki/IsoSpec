from .isoFFI import isoFFI
from collections import defaultdict


# There IS a real exported C variable with this exact value,
# isospec_number_of_isotopic_entries (element_tables.cpp), which would
# eliminate the hand-duplication below entirely if it were usable -- but
# cffi's dlopen()-based ABI mode (what isoFFI.py uses; there's no
# compile-time cffi build step here) cannot read a plain `extern const
# size_t` global: `NotImplementedError: non-integer constant
# 'isospec_number_of_isotopic_entries' cannot be accessed from a dlopen()
# library` (confirmed by trying it directly). Only #define-style integer
# constants (NUMBER_OF_ISOTOPIC_ENTRIES below) are readable this way, and
# that one is exactly the hand-typed value that must be kept in sync by hand
# with ISOSPEC_NUMBER_OF_ISOTOPIC_ENTRIES in element_tables.h -- forgetting
# to bump it already caused one real bug (the #51/deuterium fix's first
# pass silently missed "D" here even though the compiled library itself
# already had it). No independent runtime cross-check is possible in ABI
# mode, so there's nothing to assert against -- if you add table entries,
# bump NUMBER_OF_ISOTOPIC_ENTRIES in isoFFI.py's cdef by hand, in the same
# commit, or elements will silently go missing here again.
number_of_isotopic_entries = isoFFI.clib.NUMBER_OF_ISOTOPIC_ENTRIES

symbol_to_masses = defaultdict(tuple)
symbol_to_massNo = defaultdict(tuple)
symbol_to_probs  = defaultdict(tuple)
symbol_to_atomic_number = {}

for i in range(number_of_isotopic_entries):
    symbol = isoFFI.ffi.string(isoFFI.clib.elem_table_symbol[i]).decode("ascii")
    symbol_to_masses[symbol] += (isoFFI.clib.elem_table_mass[i],)
    symbol_to_massNo[symbol] += (isoFFI.clib.elem_table_massNo[i],)
    symbol_to_probs[symbol] += (isoFFI.clib.elem_table_probability[i],)
    symbol_to_atomic_number[symbol] = isoFFI.clib.elem_table_atomicNo[i]

symbol_to_masses = dict(symbol_to_masses)
symbol_to_probs = dict(symbol_to_probs)

# Several derivative convenience dicts...
symbol_to_massprob = dict((key, list(zip(symbol_to_masses[key], symbol_to_probs[key]))) for key in symbol_to_probs.keys())

def crossprod(l1, l2):
    return sum(x1*x2 for x1, x2 in zip(l1, l2))

symbol_to_avg_mass = dict((key, crossprod(symbol_to_masses[key], symbol_to_probs[key])) for key in symbol_to_probs.keys())

def maxprod(l1, l2):
    return max(zip(l1, l2), key = lambda x: x[1])[0]

symbol_to_monoisotopic_mass = dict((key, maxprod(symbol_to_masses[key], symbol_to_probs[key])) for key in symbol_to_probs.keys())

