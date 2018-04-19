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

#ifndef __IEXP_RAND_HYPER_GEOMETRIC__
#define __IEXP_RAND_HYPER_GEOMETRIC__

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

class hgeo
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
        rng(unsigned int n1,
            unsigned int n2,
            unsigned int t,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_n1(n1)
            , m_n2(n2)
            , m_t(t)
            , m_rng(type, seed)
        {
        }

        void seed(unsigned long seed)
        {
            m_rng.seed(seed);
        }

        unsigned int next()
        {
            return gsl_ran_hypergeometric(m_rng.gsl(), m_n1, m_n2, m_t);
        }

      private:
        unsigned int m_n1, m_n2, m_t;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            unsigned int n1,
                            unsigned int n2,
                            unsigned int t,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        hgeo::rng r(n1, n2, t, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, hgeo::rng &r)
        -> decltype(x.derived())
    {
        static_assert(IS_INTEGER(typename T::Scalar),
                      "only support integer scalar");

        typename T::Scalar *data = x.derived().data();
        for (Index i = 0; i < x.size(); ++i) {
            data[i] = (typename T::Scalar)r.next();
        }
        return x.derived();
    }

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(unsigned int n1, unsigned int n2, unsigned int t)
            : m_n1(n1)
            , m_n2(n2)
            , m_t(t)
        {
        }

        double pdf(unsigned int x) const
        {
            return gsl_ran_hypergeometric_pdf(x, m_n1, m_n2, m_t);
        }

        double p(unsigned int x) const
        {
            return gsl_cdf_hypergeometric_P(x, m_n1, m_n2, m_t);
        }

        double q(unsigned int x) const
        {
            return gsl_cdf_hypergeometric_Q(x, m_n1, m_n2, m_t);
        }

      private:
        unsigned int m_n1, m_n2, m_t;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &x, unsigned int n1, unsigned int n2, unsigned int t)
    {
        static_assert(IS_INTEGER(typename T::Scalar),
                      "only support integer scalar");

        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(x.rows(),
                                       x.cols(),
                                       pdf_functor<T>(x.derived(), n1, n2, t));
    }

  private:
    hgeo() = delete;

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive<T, double>::type;

        pdf_functor(const T &x,
                    unsigned int n1,
                    unsigned int n2,
                    unsigned int t)
            : m_x(x)
            , m_dist(n1, n2, t)
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

#endif /* __IEXP_RAND_HYPER_GEOMETRIC__ */
