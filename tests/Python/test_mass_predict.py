from IsoSpecPy import IsoThreshold, Iso, ParseFormula, Advanced






def check_estimates(formula, est_fun, thr):
    it = IsoThreshold(thr, formula, absolute = False)
    min_m = min(it.masses)
    max_m = max(it.masses)

    min_e_m, max_e_m = est_fun(formula, thr)

#    assert min_e_m < min_m
#    assert max_e_m > max_e

    return max_e_m, max_m, min_m, min_e_m, max_e_m > max_m >= min_m > min_e_m, (max_e_m - min_e_m) / (max_m - min_m + 0.0001), (max_e_m - min_e_m), (max_m - min_m)



def est_fun_trivial(formula, _):
    i = Iso(formula)
    return i.getLightestPeakMass(), i.getHeaviestPeakMass()


def est_fun_margs(formula, threshold):
    PC = ParseFormula(formula)
    it = IsoThreshold(threshold, formula, get_confs = True)
    marginals = [Iso(elem, cnt) for elem, cnt in zip(PC[0], PC[1])]
    mconfs = [set() for _ in marginals]
    for conf in it.confs:
        for marg_conf, marg_set in zip(conf, mconfs):
            marg_set.add(marg_conf)

    mins  = [min(marginal.conf_mass((conf,)) for conf in cset) for cset, marginal in zip(mconfs, marginals)]
    maxes = [max(marginal.conf_mass((conf,)) for conf in cset) for cset, marginal in zip(mconfs, marginals)]

    return sum(mins) - 1.0, sum(maxes) + 1.0


from IsoSpecPy.Formulas import *

for formula in "bovine_insulin horse_myoglobin surcose water ubiquitin caffeine averagine(100.0) averagine(1000.0) averagine(100000.0) Hg10Sn5 H1000N1000".split():
    try:
        fformula = eval(formula)
    except NameError:
        fformula = formula
    for thr in [0.1, 0.01, 0.001, 0.00001]:
        print(formula, thr, check_estimates(fformula, est_fun_margs, thr))

