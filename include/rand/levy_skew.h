/* Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

#ifndef __IEXP_RAND_LEVY_SKEW__
#define __IEXP_RAND_LEVY_SKEW__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class levysk
{
  public:
    // ========================================
    // generator
    // ========================================

    class rng
    {
      public:
        rng(double c,
            double alpha,
            double beta,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_c(c)
            , m_alpha(alpha)
            , m_beta(beta)
            , m_rng(type, seed)
        {
        }

        rng &seed(unsigned long seed)
        {
            m_rng.seed(seed);
            return *this;
        }

        double next()
        {
            return gsl_ran_levy_skew(m_rng.gsl(), m_c, m_alpha, m_beta);
        }

      private:
        double m_c, m_alpha, m_beta;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            double c,
                            double alpha,
                            double beta,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        levysk::rng r(c, alpha, beta, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, levysk::rng &r)
        -> decltype(x.derived())
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");

        double *data = x.derived().data();
        for (Index i = 0; i < x.size(); ++i) {
            data[i] = r.next();
        }
        return x.derived();
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_LEVY_SKEW__ */
