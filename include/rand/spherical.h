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

#ifndef __IEXP_RAND_SPHERICAL__
#define __IEXP_RAND_SPHERICAL__

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

template <bool trig>
void sph2_impl(const gsl_rng *r, double *x, double *y)
{
    gsl_ran_dir_2d(r, x, y);
}

template <>
void sph2_impl<true>(const gsl_rng *r, double *x, double *y)
{
    gsl_ran_dir_2d_trig_method(r, x, y);
}

template <bool trig = false>
class sph2_rng
{
  public:
    sph2_rng(rng_type type = DEFAULT_RNG, unsigned long seed = 0)
        : m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double &x, double &y)
    {
        return sph2_impl<trig>(m_rng.gsl(), &x, &y);
    }

  private:
    rng m_rng;
};

template <bool trig = false, typename T = void>
inline auto sph2_rand(DenseBase<T> &x,
                      unsigned long seed = 0,
                      rng_type type = MT19937) -> decltype(x.derived())
{
    static_assert(is_complex<typename T::Scalar>::value,
                  "scalar can only be complex");

    sph2_rng<trig> r(type, seed);

    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        r.next(((double *)&data[i])[0], ((double *)&data[i])[1]);
    }

    return x.derived();
}

class sph3_rng
{
  public:
    sph3_rng(rng_type type = DEFAULT_RNG, unsigned long seed = 0)
        : m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double &x, double &y, double &z)
    {
        return gsl_ran_dir_3d(m_rng.gsl(), &x, &y, &z);
    }

  private:
    rng m_rng;
};

class sphn_rng
{
  public:
    sphn_rng(size_t n, rng_type type = DEFAULT_RNG, unsigned long seed = 0)
        : m_n(n)
        , m_rng(type, seed)
    {
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double x[])
    {
        return gsl_ran_dir_nd(m_rng.gsl(), m_n, x);
    }

  private:
    const size_t m_n;
    rng m_rng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_SPHERICAL__ */
