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

class discrete
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
        rng(size_t K,
            const double P[],
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_g(nullptr)
            , m_rng(type, seed)
        {
            m_g = gsl_ran_discrete_preproc(K, P);
            IEXP_NOT_NULLPTR(m_g);
        }

        ~rng()
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
        rng(const rng &) = delete;
        rng &operator=(const rng &other) = delete;

        gsl_ran_discrete_t *m_g;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            size_t K,
                            const double P[],
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        discrete::rng r(K, P, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, discrete::rng &r)
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
        dist(size_t K, const double P[])
            : m_g(nullptr)
        {
            m_g = gsl_ran_discrete_preproc(K, P);
            IEXP_NOT_NULLPTR(m_g);
        }

        ~dist()
        {
            gsl_ran_discrete_free(m_g);
        }

        double pdf(size_t x) const
        {
            return gsl_ran_discrete_pdf(x, m_g);
        }

      private:
        dist(const dist &) = delete;
        dist &operator=(const dist &other) = delete;

        gsl_ran_discrete_t *m_g;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &x, size_t K, const double P[])
    {
        static_assert(IS_INTEGER(typename T::Scalar),
                      "only support integer scalar");

        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(x.rows(),
                                       x.cols(),
                                       pdf_functor<T>(x.derived(), K, P));
    }

  private:
    discrete() = delete;

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive<T, double>::type;

        pdf_functor(const T &x, size_t K, const double P[])
            : m_x(x)
            , m_dist(new dist(K, P))
        {
        }

        double operator()(Index i, Index j) const
        {
            return m_dist->pdf(m_x(i, j));
        }

      private:
        const T &m_x;
        std::shared_ptr<dist> m_dist;
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

#endif /* __IEXP_RAND_DISCRETE__ */
