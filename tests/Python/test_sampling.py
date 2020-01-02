import IsoSpecPy


from random import uniform, choice
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
        if self.ionic_current < 0:
            raise Exception("Ionic current < 0")
        if self.ionic_current <= 0:
            if self.accumulated > 0:
                print("LEFTOVER")
                self.current_count = self.accumulated
                self.accumulated = 0
                return True
            return False
        while self.chasing_prob >= self.accumulated_prob:
            self.mass, self.cconf_prob = next(self.iso)
            self.accumulated_prob += self.cconf_prob
        prob_diff = self.accumulated_prob - self.chasing_prob
        expected_mols = prob_diff * self.ionic_current
        rem_interval = self.accuracy - self.chasing_prob
        while True:
            if expected_mols / rem_interval < self.beta_bias:
                print("BETA")
                self.current_count = self.accumulated
                while self.ionic_current > 0 and self.chasing_prob < self.accumulated_prob:
                    self.chasing_prob += _beta_1_b(self.ionic_current) * rem_interval
                    self.ionic_current -= 1
                    self.current_count += 1
                    rem_interval = self.accuracy - self.chasing_prob
                if self.ionic_current > 0:
                    self.accumulated = 1
                    self.ionic_current -= 1
                def raisee():
                    raise Exception("current count < 0")
                return (True if self.current_count > 0 else raisee())
            else:
                print("BINOM")
                self.current_count = _safe_binom(self.ionic_current, prob_diff / rem_interval)
                self.ionic_current -= self.current_count
                self.current_count += self.accumulated
                self.accumulated = 0
                self.chasing_prob = self.accumulated_prob
                if self.current_count > 0:
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
    S = Sampler(population, count, precision, 1.0)
    while S.advance():
        yield S.current()



class CIIC:
    def __init__(self, iso, no_confs, precision = 0.9999, beta_bias = 1.0):
        self.iso = iso.__iter__()
        self.confs_prob = 0.0
        self.chasing_prob = 0.0
        self.to_sample_left = no_confs
        self.precision = precision
        self.beta_bias = beta_bias
    def next(self):
        while True:
            if self.to_sample_left <= 0:
                return False
            if self.confs_prob < self.chasing_prob:
                # Beta was last
                self.current_count = 1
                self.to_sample_left -= 1
                self.current_conf, self.current_prob = next(self.iso)
                self.confs_prob += self.current_prob
                while self.confs_prob <= self.chasing_prob:
                    self.current_conf, self.current_prob = next(self.iso)
                    self.confs_prob += self.current_prob
                if self.to_sample_left <= 0:
                    assert self.current_count > 0
                    return True
                curr_conf_prob_left = self.confs_prob - self.chasing_prob
            else:
                # Binomial was last
                self.current_count = 0
                self.current_conf, self.current_prob = next(self.iso)
                self.confs_prob += self.current_prob
                curr_conf_prob_left = self.current_prob

            assert self.to_sample_left > 0
            assert self.chasing_prob < self.confs_prob

            prob_left_to_1 = self.precision - self.chasing_prob
            expected_confs = curr_conf_prob_left * self.to_sample_left / prob_left_to_1

            if self.beta_bias < 0.0:
                cond = choice([True, False])
                print("RAND")
            else:
                cond = expected_confs <= self.beta_bias
            if cond:
                # Beta mode: we keep making beta jumps until we leave the current configuration
                self.chasing_prob += _beta_1_b(self.to_sample_left) * prob_left_to_1
                while self.chasing_prob <= self.confs_prob:
                    self.current_count += 1
                    self.to_sample_left -= 1
                    if self.to_sample_left == 0:
                        return True
                    prob_left_to_1 = self.precision - self.chasing_prob
                    self.chasing_prob += _beta_1_b(self.to_sample_left) * prob_left_to_1
                if self.current_count > 0:
                    return True
            else:
                # Binomial mode: a single binomial step
                rbin = _safe_binom(self.to_sample_left, curr_conf_prob_left/prob_left_to_1)
                self.current_count += rbin
                self.to_sample_left -= rbin
                self.chasing_prob = self.confs_prob
                if self.current_count > 0:
                    return True



def sample_ciic(formula, count, precision):
    population = IsoSpecPy.IsoLayeredGenerator(formula, t_prob_hint = precision, reorder_marginals = False)
    S = CIIC(population, count, precision, -1.0)
    while S.next():
        print(S.confs_prob, S.chasing_prob)
        yield (S.current_conf, S.current_count)



from IsoSpecPy.Formulas import *
from scipy.stats import chisquare
import sys

if __name__ == '__main__':
    test_mol = surcose
    count = 10000000

    print("Starting...")
    X = sorted(x for x in IsoSpecPy.IsoThresholdGenerator(formula=test_mol, threshold=sys.float_info.min, reorder_marginals = False) if x[1] > 0)

    print("No configs: " + str(len(X)))

    Y = dict([(v[0], 0) for v in X])
    #print(Y)

    s = 0
    for x in sample_ciic(test_mol, count, 0.999999):
        print(x)
        Y[x[0]] = x[1]
        s += x[1]
    print("S:", s)
    assert s == count

    #print(X)
    #print(Y)

    X = [x[1]*count for x in sorted(X)]
    Y = [x[1] for x in sorted(Y.items())]

    #print(X)
    #print(Y)
    print(chisquare(Y, X))
