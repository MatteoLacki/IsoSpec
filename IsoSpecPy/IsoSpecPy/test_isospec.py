from IsoSpecPy import IsoParamsFromFormula
from IsoSpecPy import isospecify


def test_IsoParamsFromFormula():
    """Test if the formulas are parsed consistently."""
    atom_cnts, iso_masses, iso_probs = IsoParamsFromFormula("C100H202")
    assert atom_cnts == [100, 202]
    # rounding to 8 digits
    iso_masses = [round(x, 8) for Ma in iso_masses for x in Ma]
    assert iso_masses == [12.0, 13.00335484, 1.00782503, 2.01410178]
    iso_probs = [round(x, 8) for Pa in iso_probs for x in Pa]
    assert iso_probs == [0.98921194, 0.01078806, 0.99988429, 0.00011571]


def test_count():
    iso = isospecify("C2H6", 1, get_counts=True)
    for i, x in enumerate(iso):
        pass
    assert i+1 == 3*7
