import math

from .IsoSpecPy import IsoDistribution




def simple_inverse(fun, prec = 0.01):
    
    def inverse(x):
        start = -1.0
        while fun(start) > x:
            start *= 2.0
        end = 1.0
        while fun(end) < x:
            end *= 2.0

        while (end - start) > prec:
            print(start, end)
            mid = (end + start) * 0.5
            if fun(mid) < x:
                start = mid
            else:
                end = mid

        return (end + start) * 0.5

    return inverse
            


class Distribution(IsoDistribution):
    def __init__(self, cdf, bin_width = 0.01, precision = 0.99, inverse_cdf = None):
        if inverse_cdf is None:
            inverse_cdf = simple_inverse(cdf, prec = bin_width*0.125)

        prec_missing = 0.5 * (1.0 - precision)

        start = inverse_cdf(prec_missing)
        end = inverse_cdf(1.0 - prec_missing)

        half_bin_width = 0.5*bin_width

        bw_inverse = 1.0 / bin_width

        def bw_round(x):
            return math.floor(x*bw_inverse + 0.5) / bw_inverse

        start = bw_round(start)
        end = bw_round(end)

        last_bin_start = start - half_bin_width
        last_cdf = cdf(last_bin_start)
        
        probs = [0.0]
        masses = [last_bin_start - half_bin_width]

        while last_bin_start < end:
            next_bin_start = last_bin_start + bin_width
            next_cdf = cdf(next_bin_start)

            masses.append(last_bin_start + half_bin_width)
            probs.append(next_cdf - last_cdf)

            last_cdf = next_cdf
            last_bin_start = next_bin_start

        sprobs = 1.0/sum(probs)
        probs = [p*sprobs for p in probs]

        probs.append(0.0)
        masses.append(last_bin_start + half_bin_width)

        super(Distribution, self).__init__(masses = masses, probs = probs)

        self.mass_sorted = True



class Gaussian(Distribution):
    def __init__(self, stdev, bin_width = 0.01, precision = 0.99):
        varfactor = 1.0/(stdev * 1.4142135623730951)
        cdf = lambda x: 0.5*(1.0+math.erf(x*varfactor))
        inv_cdf = None
        try:
            # If scipy is available, use its ppf
            from scipy.stats import norm
            inv_cdf = lambda x: norm.ppf(x, loc = 0.0, scale = stdev)
        except ImportError:
            pass # oh well, no scipy, we'll just use binsearch-based inversion

        super(Gaussian, self).__init__(cdf = cdf, bin_width = bin_width, precision = precision, inverse_cdf = inv_cdf)


