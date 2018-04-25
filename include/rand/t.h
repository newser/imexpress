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

#ifndef __IEXP_RAND_T__
#define __IEXP_RAND_T__

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

class t
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
        rng(double nu,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_nu(nu)
            , m_rng(type, seed)
        {
        }

        void seed(unsigned long seed)
        {
            m_rng.seed(seed);
        }

        double next()
        {
            return gsl_ran_tdist(m_rng.gsl(), m_nu);
        }

      private:
        double m_nu;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            double nu,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        t::rng r(nu, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, t::rng &r) -> decltype(x.derived())
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
        dist(double nu)
            : m_nu(nu)
        {
        }

        double pdf(double x) const
        {
            return gsl_ran_tdist_pdf(x, m_nu);
        }

        double p(double x) const
        {
            return gsl_cdf_tdist_P(x, m_nu);
        }

        double invp(double x) const
        {
            return gsl_cdf_tdist_Pinv(x, m_nu);
        }

        double q(double x) const
        {
            return gsl_cdf_tdist_Q(x, m_nu);
        }

        double invq(double x) const
        {
            return gsl_cdf_tdist_Qinv(x, m_nu);
        }

      private:
        double m_nu;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &x, double nu)
    {
        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(x.rows(),
                                       x.cols(),
                                       pdf_functor<T>(x.derived(), nu));
    }

  private:
    t() = delete;

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive<T, double>::type;

        pdf_functor(const T &x, double nu)
            : m_x(x)
            , m_dist(nu)
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

#endif /* __IEXP_RAND_T__ */
