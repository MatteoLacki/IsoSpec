/* This file was taken from Boost, with slight modifications
 * boost random/binomial_distribution.hpp header file, at version 1.71
 *
 * Copyright Steven Watanabe 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org for most recent version including documentation.
 *
 *
 */

#include "isoMath.h"
#include <cstdlib>
#include <cmath>


namespace IsoSpec {

typedef double RealType;
typedef ssize_t IntType;


struct binomial_table {
    static const RealType table[10];
};

const RealType binomial_table::table[10] = {
    0.08106146679532726,
    0.04134069595540929,
    0.02767792568499834,
    0.02079067210376509,
    0.01664469118982119,
    0.01387612882307075,
    0.01189670994589177,
    0.01041126526197209,
    0.009255462182712733,
    0.008330563433362871
};


/**
 * The binomial distribution is an integer valued distribution with
 * two parameters, @c t and @c p.  The values of the distribution
 * are within the range [0,t].
 *
 * The distribution function is
 * \f$\displaystyle P(k) = {t \choose k}p^k(1-p)^{t-k}\f$.
 *
 * The algorithm used is the BTRD algorithm described in
 *
 *  @blockquote
 *  "The generation of binomial random variates", Wolfgang Hormann,
 *  Journal of Statistical Computation and Simulation, Volume 46,
 *  Issue 1 & 2 April 1993 , pages 101 - 110
 *  @endblockquote
 */
class binomial_distribution {
    // parameters
    IntType _t;
    RealType _p;

    // common data
    IntType m;

    union {
        // for btrd
        struct {
            RealType r;
            RealType nr;
            RealType npq;
            RealType b;
            RealType a;
            RealType c;
            RealType alpha;
            RealType v_r;
            RealType u_rv_r;
        } btrd;
        // for inversion
        RealType q_n;
    } _u;

public:

    /**
     * Construct a @c binomial_distribution object. @c t and @c p
     * are the parameters of the distribution.
     *
     * Requires: t >=0 && 0 <= p <= 1
     */
    explicit binomial_distribution(IntType t_arg = 1,
                                   RealType p_arg = RealType(0.5))
      : _t(t_arg), _p(p_arg)
    {
        init();
    }


    /**
     * Returns a random variate distributed according to the
     * binomial distribution.
     */
    IntType operator()(std::mt19937& urng = random_gen) const
    {
        if(use_inversion()) {
            if(0.5 < _p) {
                return _t - invert(_t, 1-_p, urng);
            } else {
                return invert(_t, _p, urng);
            }
        } else if(0.5 < _p) {
            return _t - generate(urng);
        } else {
            return generate(urng);
        }
    }

private:

    inline bool use_inversion() const
    {
        // BTRD is safe when np >= 10
        return m < 11;
    }

    // computes the correction factor for the Stirling approximation
    // for log(k!)
    static RealType fc(IntType k)
    {
        if(k < 10) return binomial_table::table[k];
        else {
            RealType ikp1 = RealType(1) / (k + 1);
            return (RealType(1)/12
                 - (RealType(1)/360
                 - (RealType(1)/1260)*(ikp1*ikp1))*(ikp1*ikp1))*ikp1;
        }
    }

    void init()
    {
        using std::sqrt;
        using std::pow;

        RealType p = (0.5 < _p)? (1 - _p) : _p;
        IntType t = _t;

        m = static_cast<IntType>((t+1)*p);

        if(use_inversion()) {
            _u.q_n = pow((1 - p), static_cast<RealType>(t));
        } else {
            _u.btrd.r = p/(1-p);
            _u.btrd.nr = (t+1)*_u.btrd.r;
            _u.btrd.npq = t*p*(1-p);
            RealType sqrt_npq = sqrt(_u.btrd.npq);
            _u.btrd.b = 1.15 + 2.53 * sqrt_npq;
            _u.btrd.a = -0.0873 + 0.0248*_u.btrd.b + 0.01*p;
            _u.btrd.c = t*p + 0.5;
            _u.btrd.alpha = (2.83 + 5.1/_u.btrd.b) * sqrt_npq;
            _u.btrd.v_r = 0.92 - 4.2/_u.btrd.b;
            _u.btrd.u_rv_r = 0.86*_u.btrd.v_r;
        }
    }

    IntType generate(std::mt19937& urng = random_gen) const
    {
        using std::floor;
        using std::abs;
        using std::log;

        while(true) {
            RealType u;
            RealType v = stdunif(urng);
            if(v <= _u.btrd.u_rv_r) {
                u = v/_u.btrd.v_r - 0.43;
                return static_cast<IntType>(floor(
                    (2*_u.btrd.a/(0.5 - abs(u)) + _u.btrd.b)*u + _u.btrd.c));
            }

            if(v >= _u.btrd.v_r) {
                u = stdunif(urng) - 0.5;
            } else {
                u = v/_u.btrd.v_r - 0.93;
                u = ((u < 0)? -0.5 : 0.5) - u;
                v = stdunif(urng) * _u.btrd.v_r;
            }

            RealType us = 0.5 - abs(u);
            IntType k = static_cast<IntType>(floor((2*_u.btrd.a/us + _u.btrd.b)*u + _u.btrd.c));
            if(k < 0 || k > _t) continue;
            v = v*_u.btrd.alpha/(_u.btrd.a/(us*us) + _u.btrd.b);
            RealType km = abs(k - m);
            if(km <= 15) {
                RealType f = 1;
                if(m < k) {
                    IntType i = m;
                    do {
                        ++i;
                        f = f*(_u.btrd.nr/i - _u.btrd.r);
                    } while(i != k);
                } else if(m > k) {
                    IntType i = k;
                    do {
                        ++i;
                        v = v*(_u.btrd.nr/i - _u.btrd.r);
                    } while(i != m);
                }
                if(v <= f) return k;
                else continue;
            } else {
                // final acceptance/rejection
                v = log(v);
                RealType rho =
                    (km/_u.btrd.npq)*(((km/3. + 0.625)*km + 1./6)/_u.btrd.npq + 0.5);
                RealType t = -km*km/(2*_u.btrd.npq);
                if(v < t - rho) return k;
                if(v > t + rho) continue;

                IntType nm = _t - m + 1;
                RealType h = (m + 0.5)*log((m + 1)/(_u.btrd.r*nm))
                           + fc(m) + fc(_t - m);

                IntType nk = _t - k + 1;
                if(v <= h + (_t+1)*log(static_cast<RealType>(nm)/nk)
                          + (k + 0.5)*log(nk*_u.btrd.r/(k+1))
                          - fc(k)
                          - fc(_t - k))
                {
                    return k;
                } else {
                    continue;
                }
            }
        }
    }

    IntType invert(IntType t, RealType p, std::mt19937& urng = random_gen) const
    {
        RealType q = 1 - p;
        RealType s = p / q;
        RealType a = (t + 1) * s;
        RealType r = _u.q_n;
        RealType u = stdunif(urng);
        IntType x = 0;
        while(u > r) {
            u = u - r;
            ++x;
            RealType r1 = ((a/x) - s) * r;
            // If r gets too small then the round-off error
            // becomes a problem.  At this point, p(i) is
            // decreasing exponentially, so if we just call
            // it 0, it's close enough.  Note that the
            // minimum value of q_n is about 1e-7, so we
            // may need to be a little careful to make sure that
            // we don't terminate the first time through the loop
            // for float.  (Hence the test that r is decreasing)
            if(r1 < std::numeric_limits<RealType>::epsilon() && r1 < r) {
                break;
            }
            r = r1;
        }
        return x;
    }

};


}
