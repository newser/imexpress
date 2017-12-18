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

#ifndef __IEXP_RAND_RAYLEIGH_TAIL__
#define __IEXP_RAND_RAYLEIGH_TAIL__

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

class raylt_rng
{
  public:
    raylt_rng(double a,
              double sigma,
              rng_type type = DEFAULT_RNG,
              unsigned long seed = 0)
        : m_a(a)
        , m_sigma(sigma)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    double next()
    {
        return gsl_ran_rayleigh_tail(m_rng.gsl(), m_a, m_sigma);
    }

  private:
    double m_a, m_sigma;
    rng m_rng;
};

template <typename T>
inline auto raylt_rand(DenseBase<T> &x,
                       typename T::Scalar a,
                       typename T::Scalar sigma,
                       unsigned long seed = 0,
                       rng_type type = MT19937) -> decltype(x.derived())
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    raylt_rng r(a, sigma, type, seed);

    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        data[i] = r.next();
    }

    return x.derived();
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_RAYLEIGH_TAIL__ */
