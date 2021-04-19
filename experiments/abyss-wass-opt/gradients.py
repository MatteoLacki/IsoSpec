import numpy as np

import flows
import parameters


def flow_comp_grad(exps, thes, point, e_ab_c, t_ab_c):
    r_thes = []
    r_thes = [the.scaled(x) for the, x in zip(thes, point)]

    graph_grad = flows.awsd_g(exps, r_thes, e_ab_c, t_ab_c)[1] / np.array(point)

    return graph_grad


def flow_emp_grad(exps, thes, point, e_ab_c, t_ab_c):
    return flows.empiric_gradient(exps, thes, point, e_ab_c, t_ab_c)


def empiric_grad(awsd_fun, exp, the_l, point, exp_ab_cost, th_ab_cost):
    dval = parameters.emp_grad_dval
    res = []
    the_l = [i.copy() for i in the_l]
    point = point.copy()

    def rescaled(p):
        return [iso.scaled(xi) for iso, xi in zip(the_l, p)]

    base = awsd_fun(exp, rescaled(point), exp_ab_cost, th_ab_cost)

    for idx in range(len(the_l)):
        dpoint = point.copy()
        dpoint[idx] += dval
        grpart = (awsd_fun(exp, rescaled(dpoint), exp_ab_cost, th_ab_cost) - base) / dval
        res.append(grpart)

    return np.array(res)

def flow_emp_grad(exps, thes, point, e_ab_c, t_ab_c):
    return empiric_grad(flows.awsd, exps, thes, point, e_ab_c, t_ab_c)
    return flows.empiric_gradient(exps, thes, point, e_ab_c, t_ab_c)


def checked_grad(exps, thes, point, e_ab_c, t_ab_c):
    r_thes = []
    for the, x in zip(thes, point):
        c = the.copy()
        c.scale(x)
        r_thes.append(c)

    graph_grad = flow_comp_grad(exps, thes, point, e_ab_c, t_ab_c)
    emp_grad = flow_emp_grad(exps, thes, point, e_ab_c, t_ab_c)
    if not all(np.isclose(graph_grad, emp_grad, rtol=0.001, atol=0.000001)):
        print(graph_grad)
        print(emp_grad)
        raise Exception()
    return emp_grad


if __name__ == '__main__':
    from test_spectra import *
    import itertools

    grid = np.array([ 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 100.0,])

    # Get off region boundaries,
    # algos are allowed to return grad from different neighbouring regions there
    # (different elements of the subdifferential)
    for x, y, xscale, yscale in itertools.product(grid*1.01, grid*1.01, grid, grid):
        print(checked_grad(EXP.scaled(xscale), [t.scaled(yscale) for t in THEs], [x, y], 0.2, 0.2))

