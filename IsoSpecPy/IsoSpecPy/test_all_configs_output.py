def binom(n, k):
    """Quickly adapted from https://stackoverflow.com/questions/26560726/python-binomial-coefficient"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    total_ways = 1
    for i in range(min(k, n - k)):
        total_ways = total_ways * (n - i) // (i + 1)
    return total_ways


def max_confs_cnt(formula=""):
    """Get the maximal number of configurations for a given chemical formula."""
    from . import IsoParamsFromFormula
    f = IsoParamsFromFormula(formula)
    if f.atomCount:
        N = 1
        for n, p in zip(f.atomCount, f.prob):
            N *= binom(n+len(p)-1, n) 
        return N
    else:
        return 0


def test_max_confs_cnt():
    assert max_confs_cnt("O100") == 5151
    assert max_confs_cnt("O100N10S6") == 4759524


test_formulas = [   'O100',
                    'O100N10S6',
                    'C100H202',
                    'S10H20'        ]

def test_all_configs_output_cnt():
    """Test if IsoSpecPy output correctly all configurations."""
    from . import IsoThreshold
    global test_formulas
    for f in test_formulas:
        I = IsoThreshold(formula=f, threshold=0.0, absolute=True)
        assert len(I) == max_confs_cnt(f)