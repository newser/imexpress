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

#ifndef __IEXP_RAND_LGNORM__
#define __IEXP_RAND_LGNORM__

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

class lgnorm
{
    template <typename T>
    class pdf_functor;

  public:
    // ========================================
    // generator
    // ========================================

    class rng
    {
      public:
        rng(double zeta,
            double sigma,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_zeta(zeta)
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
            return gsl_ran_lognormal(m_rng.gsl(), m_zeta, m_sigma);
        }

      private:
        double m_zeta, m_sigma;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            double zeta,
                            double sigma,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        lgnorm::rng r(zeta, sigma, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, lgnorm::rng &r)
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

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(double zeta, double sigma)
            : m_zeta(zeta)
            , m_sigma(sigma)
        {
        }

        double pdf(double x) const
        {
            return gsl_ran_lognormal_pdf(x, m_zeta, m_sigma);
        }

        double p(double x) const
        {
            return gsl_cdf_lognormal_P(x, m_zeta, m_sigma);
        }

        double invp(double x) const
        {
            return gsl_cdf_lognormal_Pinv(x, m_zeta, m_sigma);
        }

        double q(double x) const
        {
            return gsl_cdf_lognormal_Q(x, m_zeta, m_sigma);
        }

        double invq(double x) const
        {
            return gsl_cdf_lognormal_Qinv(x, m_zeta, m_sigma);
        }

      private:
        double m_zeta, m_sigma;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &x, double zeta, double sigma)
    {
        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(x.rows(),
                                       x.cols(),
                                       pdf_functor<T>(x.derived(),
                                                      zeta,
                                                      sigma));
    }

  private:
    lgnorm() = delete;

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive<T, double>::type;

        pdf_functor(const T &x, double zeta, double sigma)
            : m_x(x)
            , m_dist(zeta, sigma)
        {
        }

        double operator()(Index i, Index j) const
        {
            return m_dist.pdf(m_x(i, j));
        }

      private:
        const T &m_x;
        dist m_dist;
    };
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_LGNORM__ */
