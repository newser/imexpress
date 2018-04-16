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

#ifndef __IEXP_RAND_DIRICHLET__
#define __IEXP_RAND_DIRICHLET__

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

class drch_rng
{
  public:
    drch_rng(size_t K,
             const double alpha[],
             rng::type type = DEFAULT_RNG_TYPE,
             unsigned long seed = 0)
        : m_K(K)
        , m_alpha(alpha)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double theta[])
    {
        gsl_ran_dirichlet(m_rng.gsl(), m_K, m_alpha, theta);
    }

    size_t K() const
    {
        return m_K;
    }

  private:
    size_t m_K;
    const double *m_alpha;
    rng m_rng;
};

template <typename T>
inline auto drch_rand(DenseBase<T> &x,
                      size_t K,
                      typename T::Scalar *alpha,
                      unsigned long seed = 0,
                      rng::type type = DEFAULT_RNG_TYPE)
    -> decltype(x.derived())
{
    drch_rng r(K, alpha, type, seed);
    return drch_rand(x, r);
}

template <typename T>
inline auto drch_rand(DenseBase<T> &x, drch_rng &r) -> decltype(x.derived())
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");
    eigen_assert((x.size() % r.K()) == 0);

    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); i += r.K()) {
        r.next(&data[i]);
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

#endif /* __IEXP_RAND_DIRICHLET__ */
