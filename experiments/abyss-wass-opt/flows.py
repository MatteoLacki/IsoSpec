import networkx
import IsoSpecPy
import math
import sys
import numpy as np
from pprint import pprint

from parameters import int_fact, integerize, flows_loud


def awsd_g(exp, the_l, exp_ab_cost, th_ab_cost):
    exp_ab_cost = int(exp_ab_cost * int_fact)
    th_ab_cost = int(th_ab_cost * int_fact)
    G = networkx.DiGraph()

    exp_sum = sum(x[1] for x in integerize(exp))
    the_sum = sum((sum(x[1] for x in integerize(the_s))) for the_s in the_l)

    for mass, prob in integerize(exp):
        G.add_edge('source', (mass, -1), capacity = prob, weight = 0)
        G.add_edge((mass, -1), 'middle', weight = exp_ab_cost)
#        if exp_sum > the_sum:
#            G.add_edge((mass, -1), 'sink', weight = exp_ab_cost)

    for idx, the_iso in enumerate(the_l):
        for mass, prob in integerize(the_iso):
            G.add_edge((mass, idx), 'sink', capacity = prob, weight = 0)
            G.add_edge('middle', (mass, idx), weight = th_ab_cost)
#            if exp_sum <= the_sum:
#                G.add_edge('source', (mass, idx), weight = th_ab_cost)

    for emass, eprob in integerize(exp):
        for idx, the_iso in enumerate(the_l):
            for tmass, tprob in integerize(the_iso):
                G.add_edge((emass, -1), (tmass, idx), weight = abs(emass - tmass)) #, exp_ab_cost + th_ab_cost)
#                G.add_edge((tmass, idx), (emass, -1), weight = abs(emass - tmass))

    G.add_edge('source', 'middle', capacity = the_sum, weight = 0)
    G.add_edge('middle', 'sink', capacity = exp_sum, weight = 0)

    flow = networkx.max_flow_min_cost(G, 'source', 'sink')

#    print(flow)

    tot_cost = 0
    gradient_exp = 0
    gradient_thr = [0] * len(the_l)
    gradient_exp = 0
    gradient_thr = [0] * len(the_l)

    for start_n, neigh in flow.items():
        for end_n, fl_v in neigh.items():
            # Cost computation
            tot_cost += fl_v * G.edges[start_n, end_n]['weight']
            capacity = G.edges[start_n, end_n].get('capacity', math.inf)
#            print(start_n, end_n, "flow:", fl_v, "weight:", G.edges[start_n, end_n]['weight'], "capacity:", G.edges[start_n, end_n].get('capacity', 'inf'))

            # Grad computation

    def has_spare_flow(exp_peak):
        if exp_peak == 'middle':
            return True
        return flow[exp_peak]['middle'] > 0.0


#    pprint(flow)

    masses2probs = [{mass:prob for mass, prob in integerize(the)} for the in the_l]

    grad = np.zeros(len(the_l), dtype=float)
    for th in G.predecessors('sink'):
        if th != 'middle':
            if flow['middle'][th] > 0.0:
                grad[th[1]] += th_ab_cost * masses2probs[th[1]][th[0]]
            else:
                delta = min(
                                G[exp_p][th]['weight'] if th != 'middle'
                                else math.inf
                                for exp_p in G.predecessors(th)
                            )
                #print("Delta:", delta)

                grad[th[1]] += (delta - exp_ab_cost) * masses2probs[th[1]][th[0]]
                

    return tot_cost/int_fact/int_fact, grad/int_fact/int_fact

def awsd(exp, the_l, exp_ab_cost, th_ab_cost):
    return awsd_g(exp, the_l, exp_ab_cost, th_ab_cost)[0]

def empiric_gradient(exp, the_l, point, exp_ab_cost, th_ab_cost):
    dval = 0.0001
    res = []
    assert dval * int_fact > 1.0
    the_l = [i.copy() for i in the_l]
    point = point.copy()

    def rescaled(p):
        res = []
        for iso, xi in zip(the_l, p):
            iso = iso.copy()
            iso.scale(xi)
            res.append(iso)
        return res
        

    base = awsd(exp, rescaled(point), exp_ab_cost, th_ab_cost)
#    print("base:", base)

    for idx in range(len(the_l)):
        dpoint = point.copy()
        dpoint[idx] += dval
        grpart = (awsd(exp, rescaled(dpoint), exp_ab_cost, th_ab_cost) - base) / dval
#        print("grpart", awsd(exp, n_the_l, exp_ab_cost, th_ab_cost))
        res.append(grpart)

    return np.array(res)


def check_grad(exps, thes, point, e_ab_c, t_ab_c):
    r_thes = []
    for the, x in zip(thes, point):
        c = the.copy()
        c.scale(x)
        r_thes.append(c)
        
    graph_grad = awsd_g(exps, r_thes, e_ab_c, t_ab_c)[1] / np.array(point)
    emp_grad = empiric_gradient(exps, thes, point, e_ab_c, t_ab_c)
    if not all(np.isclose(graph_grad, emp_grad)):
        print(graph_grad)
        print(emp_grad)
        raise Exception()
    return emp_grad

if __name__ == '__main__':
    from test_spectra import *


    print(check_grad(EXP, THEs, [0.3, 1.0], 0.2, 0.2))


