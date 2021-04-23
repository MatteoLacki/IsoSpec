import math
import numpy as np

from IsoSpecPy import IsoDistribution

from flows import awsd


def awsd_checked(exp, the_l, exp_ab_cost, th_ab_cost):
    val1 = awsd(exp, the_l, exp_ab_cost, th_ab_cost)
    val2 = exp.abyssalWassersteinDistance(IsoDistribution.LinearCombination(the_l, [1.0]*len(the_l)), (exp_ab_cost + th_ab_cost))
    val3 = exp.abyssalWassersteinDistanceGrad(the_l, [1.0] * (1+len(the_l)), [1.0] * (1+len(the_l)), exp_ab_cost, th_ab_cost)

    if not (math.isclose(val1, val2, rel_tol=1e-05) and math.isclose(val1, val3, rel_tol=1e-05)):
        print("AAA", val1, val2, val3)
        raise Exception()
    return val1
