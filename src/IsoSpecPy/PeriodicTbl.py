try:
    from .isoFFI import isoFFI

    _backend = "cffi"
except ImportError:
    # Try loading from nanobind backend
    try:
        from . import _isospec_nb

        isoFFI = None
        _backend = "nanobind"
    except ImportError:
        raise ImportError(
            "Cannot load periodic table data: neither CFFI nor nanobind backend available"
        )

from collections import defaultdict


if _backend == "cffi":
    number_of_isotopic_entries = isoFFI.clib.NUMBER_OF_ISOTOPIC_ENTRIES

    symbol_to_masses = defaultdict(tuple)
    symbol_to_massNo = defaultdict(tuple)
    symbol_to_probs = defaultdict(tuple)
    symbol_to_atomic_number = {}

    for i in range(number_of_isotopic_entries):
        symbol = isoFFI.ffi.string(isoFFI.clib.elem_table_symbol[i]).decode("ascii")
        symbol_to_masses[symbol] += (isoFFI.clib.elem_table_mass[i],)
        symbol_to_massNo[symbol] += (isoFFI.clib.elem_table_massNo[i],)
        symbol_to_probs[symbol] += (isoFFI.clib.elem_table_probability[i],)
        symbol_to_atomic_number[symbol] = isoFFI.clib.elem_table_atomicNo[i]

    symbol_to_masses = dict(symbol_to_masses)
    symbol_to_probs = dict(symbol_to_probs)
else:  # nanobind backend
    number_of_isotopic_entries = _isospec_nb.NUMBER_OF_ISOTOPIC_ENTRIES

    symbol_to_masses = defaultdict(tuple)
    symbol_to_massNo = defaultdict(tuple)
    symbol_to_probs = defaultdict(tuple)
    symbol_to_atomic_number = {}

    for i in range(number_of_isotopic_entries):
        symbol = _isospec_nb.elem_table_symbol(i)
        symbol_to_masses[symbol] += (_isospec_nb.elem_table_mass(i),)
        symbol_to_massNo[symbol] += (_isospec_nb.elem_table_massNo(i),)
        symbol_to_probs[symbol] += (_isospec_nb.elem_table_probability(i),)
        symbol_to_atomic_number[symbol] = _isospec_nb.elem_table_atomicNo(i)

    symbol_to_masses = dict(symbol_to_masses)
    symbol_to_probs = dict(symbol_to_probs)

# Several derivative convenience dicts...
symbol_to_massprob = dict(
    (key, [zip(symbol_to_masses[key], symbol_to_probs[key])])
    for key in symbol_to_probs.keys()
)


def crossprod(l1, l2):
    return sum(x1 * x2 for x1, x2 in zip(l1, l2))


symbol_to_avg_mass = dict(
    (key, crossprod(symbol_to_masses[key], symbol_to_probs[key]))
    for key in symbol_to_probs.keys()
)


def maxprod(l1, l2):
    return max(zip(l1, l2), key=lambda x: x[1])[0]


symbol_to_monoisotopic_mass = dict(
    (key, maxprod(symbol_to_masses[key], symbol_to_probs[key]))
    for key in symbol_to_probs.keys()
)
