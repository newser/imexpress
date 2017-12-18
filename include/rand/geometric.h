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

#ifndef __IEXP_RAND_GEOMETRIC__
#define __IEXP_RAND_GEOMETRIC__

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

class geo_rng
{
  public:
    geo_rng(double p, rng_type type = DEFAULT_RNG, unsigned long seed = 0)
        : m_p(p)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    unsigned int next()
    {
        return gsl_ran_geometric(m_rng.gsl(), m_p);
    }

  private:
    double m_p;
    rng m_rng;
};

template <typename T>
inline auto geo_rand(DenseBase<T> &x,
                     double p,
                     unsigned long seed = 0,
                     rng_type type = DEFAULT_RNG) -> decltype(x.derived())
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    geo_rng r(p, type, seed);

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

#endif /* __IEXP_RAND_GEOMETRIC__ */
