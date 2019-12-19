import IsoSpecPy


from random import uniform
from numpy.random import binomial
from math import inf, nan


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


class Sampler:
    def __init__(self, iso, ionic_current, accuracy, beta_bias):
        self.iso = iso.__iter__()
        self.accumulated_prob = 0.0
        self.chasing_prob = 0.0
        self.ionic_current = ionic_current
        self.accuracy = accuracy
        self.beta_bias = beta_bias
        self.accumulated = 0
        self.current_count = nan
    def advance(self):
        if self.ionic_current == 0:
            return False
        while self.chasing_prob >= self.accumulated_prob:
            self.mass, self.cconf_prob = next(self.iso)
            self.accumulated_prob += self.cconf_prob
        prob_diff = self.accumulated_prob - self.chasing_prob
        expected_mols = prob_diff * self.ionic_current
        rem_interval = self.accuracy - self.chasing_prob
        while True:
            print ("Beta appr:", expected_mols / rem_interval)
            if expected_mols / rem_interval < self.beta_bias:
                self.current_count = self.accumulated
                self.ionic_current -= self.accumulated
                self.accumulated = 1
                while self.ionic_current > 0 and self.chasing_prob < self.accumulated_prob:
                    self.chasing_prob += _beta_1_b(self.ionic_current) * rem_interval
                    self.ionic_current -= 1
                    self.current_count += 1
                    rem_interval = self.accuracy - self.chasing_prob
                return self.current_count > 0
            else:
                self.current_count = _safe_binom(self.ionic_current, prob_diff / rem_interval) + self.accumulated
                self.accumulated = 0
                self.chasing_prob = self.accumulated_prob
                if self.current_count > 0:
                    self.ionic_current -= self.current_count
                    return True
                self.mass, self.cconf_prob = next(self.iso)
                self.accumulated_prob += self.cconf_prob
                prob_diff = self.cconf_prob
                expected_mols = prob_diff * self.ionic_current
                rem_interval = self.accuracy - self.chasing_prob


    def current(self):
        return (self.mass, self.current_count)



def sample_isospec2(formula, count, precision):
    population = IsoSpecPy.IsoLayeredGenerator(formula, t_prob_hint = precision, reorder_marginals = False)
    S = Sampler(population, count, precision, 100.0)
    while S.advance():
        yield S.current()



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

    s = 0
    for x in sample_isospec2(test_mol, count, 1.0):
        print(x)
        Y[x[0]] = x[1]
        s += x[1]
    print(s)

    #print(X)
    #print(Y)

    X = [x[1]*count for x in sorted(X)]
    Y = [x[1] for x in sorted(Y.items())]

    #print(X)
    #print(Y)
    print(chisquare(Y, X))
