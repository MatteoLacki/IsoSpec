import IsoSpecPy


EXP = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [3.0, 7.0])
THE1 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [1.0, 0.0])

THE2 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.0], probs = [0.0, 1.0])
THE3 = IsoSpecPy.IsoDistribution(masses=[1.0, 2.01], probs = [0.0, 1.0])

THEs = [THE1, THE2]#, THE3]



def get_probs(N):
    res = [random.randint(1,10) for x in range(N)]
    #res = probs + [N*5 - sum(probs)]
    #print("Probs:", res)
    return res
def get_masses(N):
    res = list(np.random.choice(range(N*3), size=N, replace=False))
#        res = [3,3]
    while len(res) != len(set(res)):
        res = [random.randint(0,3+N) for x in range(N)]
    #print("Masses:", res)
    return res

def random_spectrum(N):
    return IsoSpecPy.IsoDistribution(masses=get_masses(N), probs = get_probs(N))

