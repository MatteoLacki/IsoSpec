import sys
import IsoSpecPy
from math import exp, log

def sprint(s):
    sys.stdout.write(str(s))
    sys.stdout.flush()


silentish_run = False

if len(sys.argv) < 2:
    try:
        import OldIsoSpecPy
    except ImportError:
        print("This test compares the results of installed IsoSpec with IsoSpec version 1.0.7, installed as OldIsoSpecPy")
        print("You must install it, using:")
        print("pip install OldIsoSpecPy --index-url https://test.pypi.org/simple/")
        print("Or run only the new algo: python " + sys.argv[0] + " True")
        sys.exit(1)
else:
    print("Running only new algo...")
    silentish_run = True

# Correctness tests comparing IsoSpecPy 1.0.7 and HEAD. The trouble here is that due to slightly different way things are calculated
# we get different rounding errors, and therefore slightly different results - and that's OK. The function kinda_like reflects that
# - but is still WAAAY too strict. Often we can justifiably get different configurations, or different counts of configurations...

molecules = "H2O1 C100 P1 P100 C1 H10C10O10N10S5 Se1 Se10 Sn1 Sn4 Sn4C1 C2H6O1 C1000 C520H817N139O147S8 C1H1O2N2Se1Sn1P1 P1C1Sn1 Se5 Sn5 Se50 Sn15 Se2Sn2C2O2N2S2B2He2U2Na2Cl2".split()

parameters = list(map(float, "0.0 0.1 0.5 0.01 0.9 0.99 0.01 0.0001 0.999 0.362 0.852348".split()))

def kinda_like(o1, o2):
    if type(o1) in (list, tuple) and type(o2) in (list, tuple) :
        assert len(o1) == len(o2)
        assert all(kinda_like(oo1, oo2) for oo1, oo2 in zip(o1, o2))
    if type(o1) == type(o2) == float:
        assert o1*o2 >= 0.0 # same sign check
        assert abs(o1*0.99)-0.000001 <= abs(o2) <= abs(o1*1.01)+0.000001
    return True

def sort_confs(confs):
    if len(confs[0]) == 0:
        return confs
    assert len(set(tuple(x) for x in confs[2])) == len(confs[2])
    l = list(zip(*confs))
    l.sort(key = lambda x: -x[1])

    return ([x[0] for x in l], [x[1] for x in l], [x[2] for x in l])

def confs_from_ordered_generator(formula, target_prob):
    ret = ([], [], [])
    prob = 0.0
    for conf in IsoSpecPy.IsoOrderedGenerator(formula=formula, get_confs=True):
        conf = (conf[0], log(conf[1]), conf[2])
        if prob >= target_prob and target_prob < 1.0:
            return ret
        ret[0].append(conf[0])
        prob += exp(conf[1])
        ret[1].append(conf[1])
        ret[2].append([item for sublist in conf[2] for item in sublist])
    return ret

def confs_from_layered_generator(formula, target_prob):
    ret = ([], [], [])
    for conf in IsoSpecPy.IsoLayered(formula=formula, prob_to_cover = target_prob, get_confs=True, get_minimal_pset=True):
        conf = (conf[0], log(conf[1]), conf[2])
        ret[0].append(conf[0])
        ret[1].append(conf[1])
        ret[2].append([item for sublist in conf[2] for item in sublist])

    return sort_confs(ret)

def confs_from_threshold_generator(formula, target_prob):
    ret = ([], [], [])
    for conf in IsoSpecPy.IsoThresholdGenerator(formula=formula, threshold = target_prob, absolute = True, get_confs=True):
        conf = (conf[0], log(conf[1]), conf[2])
        ret[0].append(conf[0])
        ret[1].append(conf[1])
        ret[2].append([item for sublist in conf[2] for item in sublist])

    return sort_confs(ret)

def confs_from_threshold(formula, target_prob):
    ret = ([], [], [])
    for conf in IsoSpecPy.IsoThreshold(formula=formula, threshold = target_prob, absolute = True, get_confs=True):
        conf = (conf[0], log(conf[1]), conf[2])
        ret[0].append(conf[0])
        ret[1].append(conf[1])
        ret[2].append([item for sublist in conf[2] for item in sublist])

    return sort_confs(ret)



is_ok = False
try:
    i = IsoSpecPy.IsoThreshold(0.1, atomCounts = [100], isotopeMasses = [[1.0, 2.0, 3.0]], isotopeProbabilities = [[0.0, 0.6, 0.4]])
    for x in i:
        print(x)
except ValueError:
    is_ok = True
assert is_ok

total_confs = 0

for molecule in molecules:
    for parameter in parameters:
        if not silentish_run:
            sprint("{} {}... ".format(molecule, parameter))
            old_ordered = OldIsoSpecPy.IsoSpecPy.IsoSpec.IsoFromFormula(molecule, parameter, method="ordered").getConfs()
            sprint(len(old_ordered[0]))
        new_ordered = confs_from_ordered_generator(molecule, parameter)
        total_confs += len(new_ordered[0])
        if not silentish_run: assert kinda_like(new_ordered, old_ordered)
        new_layered = confs_from_layered_generator(molecule, parameter)
        total_confs += len(new_layered[0])
        if not silentish_run: assert kinda_like(new_layered, new_ordered)
        if len(new_ordered[1]) > 0:
            new_threshold = exp(new_ordered[1][-1]-0.00000001)
        else:
            new_threshold = 1.1
        new_threshold_res = confs_from_threshold_generator(molecule, new_threshold)
        total_confs += len(new_threshold_res[0])
        if not silentish_run:
            assert kinda_like(new_ordered, new_threshold_res)
            assert kinda_like(new_threshold_res, confs_from_threshold(molecule, new_threshold))

        if parameter > 0:
            if not silentish_run: sprint(" thresholded: ")
            new_thr_r = confs_from_threshold_generator(molecule, parameter)
            if not silentish_run:
                sprint(len(new_thr_r[0]))
                old_thr_r = sort_confs(OldIsoSpecPy.IsoSpecPy.IsoSpec.IsoFromFormula(molecule, parameter, method="threshold_absolute").getConfs())
                assert kinda_like(new_thr_r, old_thr_r)

        if not silentish_run:
            print("... OK!")

sprint("Total confs: " + str(total_confs) + "\n")
