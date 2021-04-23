import random
import numpy as np
import IsoSpecPy


EXP = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [3.0, 7.0])
THE1 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 0.0])

THE2 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [0.0, 1.0])
THE3 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.01], probs = [0.0, 1.0])

THEs = [THE1, THE2]#, THE3]


EXPa = IsoSpecPy.IsoDistribution(masses=[1.0, 5.0], probs=[9.0, 9.0])
THEsa = IsoSpecPy.IsoDistribution(masses=[4.0, 0.0], probs=[5.0, 4.0])
THEsa.sort_by_mass()

std_dim = 2
no_thes = 1

def get_probs(N = std_dim):
    res = [random.randint(1,10) for x in range(N)]
    #res = probs + [N*5 - sum(probs)]
    #print("Probs:", res)
    return res
def get_masses(N = std_dim):
    res = list(np.random.choice(range(N*3), size=N, replace=False))
#        res = [3,3]
    while len(res) != len(set(res)):
        res = [random.randint(0,3+N) for x in range(N)]
    #print("Masses:", res)
    return res

def random_spectrum(N = std_dim):
    return IsoSpecPy.IsoDistribution(masses=get_masses(N), probs = get_probs(N))

def random_point(N = no_thes):
    return [random.randint(1, 1000) / 100.0 for _ in range(N)]

point = random_point()
LEXP = random_spectrum(std_dim)
LTHEs = [random_spectrum(std_dim) for x in range(no_thes)]


def p(spectrum):
    print(list(spectrum.masses), list(spectrum.probs))

def rerandomize():
    LEXP = random_spectrum(std_dim)
    p(LEXP)
    LTHEs = [random_spectrum(std_dim) for x in range(no_thes)]
    print("THEs:")
    for x in LTHEs: p(x)
    
