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

#ifndef __IEXP_RAND_BI_NORMAL__
#define __IEXP_RAND_BI_NORMAL__

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

class bnorm_rng
{
  public:
    bnorm_rng(double sigma_x,
              double sigma_y,
              double rho,
              rng_type type = DEFAULT_RNG,
              unsigned long seed = 0)
        : m_sigma_x(sigma_x)
        , m_sigma_y(sigma_y)
        , m_rho(rho)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double &x, double &y)
    {
        return gsl_ran_bivariate_gaussian(m_rng.gsl(),
                                          m_sigma_x,
                                          m_sigma_y,
                                          m_rho,
                                          &x,
                                          &y);
    }

  private:
    double m_sigma_x, m_sigma_y, m_rho;
    rng m_rng;
};

template <typename T>
inline auto bnorm_rand(DenseBase<T> &x,
                       typename SCALAR(typename T::Scalar) sigma_x,
                       typename SCALAR(typename T::Scalar) sigma_y,
                       typename SCALAR(typename T::Scalar) rho,
                       unsigned long seed = 0,
                       rng_type type = MT19937) -> decltype(x.derived())
{
    static_assert(is_complex<typename T::Scalar>::value,
                  "scalar can only be complex");

    bnorm_rng r(sigma_x, sigma_y, rho, type, seed);

    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        r.next(((double *)&data[i])[0], ((double *)&data[i])[1]);
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

#endif /* __IEXP_RAND_BI_NORMAL__ */
