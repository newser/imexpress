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

#ifndef __IEXP_RAND_FLAT__
#define __IEXP_RAND_FLAT__

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

class flat_rng
{
  public:
    flat_rng(double a,
             double b,
             rng::type type = DEFAULT_RNG_TYPE,
             unsigned long seed = 0)
        : m_a(a)
        , m_b(b)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    double next()
    {
        return gsl_ran_flat(m_rng.gsl(), m_a, m_b);
    }

  private:
    double m_a, m_b;
    rng m_rng;
};

template <typename T>
inline auto flat_rand(DenseBase<T> &x,
                      typename T::Scalar a,
                      typename T::Scalar b,
                      unsigned long seed = 0,
                      rng::type type = DEFAULT_RNG_TYPE)
    -> decltype(x.derived())
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    flat_rng r(a, b, type, seed);

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

#endif /* __IEXP_RAND_FLAT__ */
