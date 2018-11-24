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


IsoThresholdGenerator(.99, )


iso = IsoLayered(prob_to_cover=.99, formula="C100H202")
for mass, log_prob in iso:
    print(mass, log_prob)


isothr = IsoThreshold(.01, get_confs=True, formula="C100H202")
for x in isothr:
    print(x)

iso_thr_gen = IsoThresholdGenerator(threshold = .01,
                                    formula   ="C100H202",
                                    get_confs = True)
for x in iso_thr_gen:
    print(x)


iso_lay_gen = IsoLayeredGenerator(prob_to_cover = .99,
                                  formula       ="C100H202",
                                  get_confs     = True)
for x in iso_lay_gen:
    print(x)

iso_ord_gen = IsoOrderedGenerator(formula       ="C100H202",
                                  get_confs     = True)
for x in iso_ord_gen:
    print(x[1] == 0.0)

iso = isospecify("C100H202", .99, get_counts=True)
for x in iso:
    print(x)

