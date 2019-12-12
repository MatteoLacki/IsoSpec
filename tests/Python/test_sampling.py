import IsoSpecPy


from random import uniform
from numpy.random import binomial


def _beta_1_b(b):
    '''Returns a random variate from beta(1, b) distribution
    using the inverse CDF method'''
    return 1.0-uniform(0.0, 1.0)**(1.0/b)

def _safe_binom(n, p):
    '''Draw a sample from binomial distribution. Doesn't crash
    if p > 1.0 due to numerical inaccuracies.'''
    if p >= 1.0:
        return n
    return binomial(n, p)





#def _sample_with_replacement_online_impl(population, probabilities, sample_size):
def sample_isospec(formula, count, precision):
    population = IsoSpecPy.IsoLayeredGenerator(formula, t_prob_hint = precision, reorder_marginals = False)
    #population = IsoSpecPy.IsoThresholdGenerator(formula = formula, threshold = -1.0)

    #for x in population:
    #    yield x
    '''Performs sampling with replacement from population argument, with
    associated probabilities from second argument. The probabilities must 
    sum to 1. Yields a stream of tuples: (population_member, times_chosen).
    Accepts generators as first and second argument. May return duplicate
    tuples and tuples with times_chosen == 0.
    '''
    pprob = 0.0
    cprob = 0.0
    accumulated = 0
    iso_iter = population.__iter__()
    while count > 0:
        if accumulated > 0:
            yield (pop_next, accumulated)
            accumulated = 0
        pop_next, prob_next = next(iso_iter)
        pprob += prob_next
        # Beta mode
        while (pprob - cprob) * count / (1.0 - cprob) < 1.0:
            cprob += _beta_1_b(count) * (1.0 - cprob)
            while pprob < cprob:
                if accumulated > 0: 
                    yield (pop_next, accumulated)
                    accumulated = 0
                pop_next, prob_next = next(iso_iter)
                pprob += prob_next
            accumulated += 1
            count -= 1
            if count == 0: break
        if count == 0: break
        # Binomial mode
        nrtaken = _safe_binom(count, (pprob-cprob)/(1.0-cprob))
        accumulated += nrtaken
        count -= nrtaken
        cprob = pprob
    if accumulated > 0:
        yield (pop_next, accumulated)


'''

    population_iter = population.__iter__()
    probabilities_iter = probabilities.__iter__()
    population_next = next(population_iter)
    while sample_size > 0:
        pprob += next(probabilities_iter)
        # Beta mode
        while (pprob - cprob) * sample_size / (1.0 - cprob) < 1.0:
            cprob += _beta_1_b(sample_size) * (1.0 - cprob)
            while pprob < cprob:
                population_next = next(population_iter)
                pprob += next(probabilities_iter)
#            yield (population_next, 1)
            sample_size -= 1
            if sample_size == 0: break
        if sample_size == 0: break
        # Binomial mode
        nrtaken = _safe_binom(sample_size, (pprob-cprob)/(1.0-cprob))
        if nrtaken > 0:
#            yield (population_next, nrtaken)
            sample_size -= nrtaken
        population_next = next(population_iter)
        cprob = pprob
'''

from IsoSpecPy.Formulas import *
from scipy.stats import chisquare
import sys

if __name__ == '__main__':
    test_mol = surcose
    count = 100000000

    print("Starting...")
    X = sorted(x for x in IsoSpecPy.IsoThresholdGenerator(formula=test_mol, threshold=sys.float_info.min, reorder_marginals = False) if x[1] > 0)

    print("No configs: " + str(len(X)))

    Y = dict([(v[0], 0) for v in X])
    #print(Y)

    
    for x in sample_isospec(test_mol, count, 0.9999):
        Y[x[0]] = x[1]

    #print(X)
    #print(Y)

    X = [x[1]*count for x in sorted(X)]
    Y = [x[1] for x in sorted(Y.items())]

    #print(X)
    #print(Y)
    print(chisquare(Y, X))
