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

#ifndef __IEXP_RAND_BI_GAUSS__
#define __IEXP_RAND_BI_GAUSS__

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

class bgauss
{
  public:
    // ========================================
    // generator
    // ========================================

    class rng
    {
      public:
        rng(double sigma_x,
            double sigma_y,
            double rho,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_sigma_x(sigma_x)
            , m_sigma_y(sigma_y)
            , m_rho(rho)
            , m_rng(type, seed)
        {
        }

        rng &seed(unsigned long seed)
        {
            m_rng.seed(seed);
            return *this;
        }

        rng &next(double &x, double &y)
        {
            gsl_ran_bivariate_gaussian(m_rng.gsl(),
                                       m_sigma_x,
                                       m_sigma_y,
                                       m_rho,
                                       &x,
                                       &y);
            return *this;
        }

      private:
        double m_sigma_x, m_sigma_y, m_rho;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            double sigma_x,
                            double sigma_y,
                            double rho,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        bgauss::rng r(sigma_x, sigma_y, rho, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static auto fill(DenseBase<T> &x, bgauss::rng &r) -> decltype(x.derived())
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");
        eigen_assert((x.size() % 2) == 0);

        typename T::Scalar *data = x.derived().data();
        for (Index i = 0; i < x.size(); i += 2) {
            r.next(data[i], data[i + 1]);
        }
        return x.derived();
    }

    template <typename T>
    static auto fill(DenseBase<T> &x,
                     DenseBase<T> &y,
                     double sigma_x,
                     double sigma_y,
                     double rho,
                     unsigned long seed = 0,
                     rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        bgauss::rng r(sigma_x, sigma_y, rho, type, seed);
        return fill(x, y, r);
    }

    template <typename T>
    static auto fill(DenseBase<T> &x, DenseBase<T> &y, bgauss::rng &r)
        -> decltype(x.derived())
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");
        eigen_assert(x.size() == y.size());

        double *data_x = x.derived().data();
        double *data_y = y.derived().data();
        for (Index i = 0; i < x.size(); ++i) {
            r.next(data_x[i], data_y[i]);
        }
        return x.derived();
    }

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(double sigma_x, double sigma_y, double rho)
            : m_sigma_x(sigma_x)
            , m_sigma_y(sigma_y)
            , m_rho(rho)
        {
        }

        double pdf(double x, double y) const
        {
            return gsl_ran_bivariate_gaussian_pdf(x,
                                                  y,
                                                  m_sigma_x,
                                                  m_sigma_y,
                                                  m_rho);
        }

      private:
        double m_sigma_x, m_sigma_y, m_rho;
    };

  private:
    bgauss() = delete;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_BI_GAUSS__ */
