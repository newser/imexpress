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

#ifndef __IEXP_RAND_DISCRETE__
#define __IEXP_RAND_DISCRETE__

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

class discrete_rng
{
  public:
    discrete_rng(size_t K,
                 const double P[],
                 rng::type type = DEFAULT_RNG_TYPE,
                 unsigned long seed = 0)
        : m_g(nullptr)
        , m_rng(type, seed)
    {
        m_g = gsl_ran_discrete_preproc(K, P);
        IEXP_NOT_NULLPTR(m_g);
    }

    ~discrete_rng()
    {
        gsl_ran_discrete_free(m_g);
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    size_t next()
    {
        return gsl_ran_discrete(m_rng.gsl(), m_g);
    }

  private:
    discrete_rng(const discrete_rng &) = delete;
    discrete_rng &operator=(const discrete_rng &other) = delete;

    gsl_ran_discrete_t *m_g;
    rng m_rng;
};

template <typename T>
inline auto discrete_rand(DenseBase<T> &x,
                          size_t K,
                          const double P[],
                          unsigned long seed = 0,
                          rng::type type = DEFAULT_RNG_TYPE)
    -> decltype(x.derived())
{
    discrete_rng r(K, P, type, seed);
    return discrete_rand(x, r);
}

template <typename T>
inline auto discrete_rand(DenseBase<T> &x, discrete_rng &r)
    -> decltype(x.derived())
{
    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        data[i] = static_cast<typename T::Scalar>(r.next());
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

#endif /* __IEXP_RAND_DISCRETE__ */
