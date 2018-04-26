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

class drch
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
        rng(size_t k,
            const double alpha[],
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_k(k)
            , m_alpha(alpha)
            , m_rng(type, seed)
        {
        }

        rng &seed(unsigned long seed)
        {
            m_rng.seed(seed);
            return *this;
        }

        rng &next(double theta[])
        {
            gsl_ran_dirichlet(m_rng.gsl(), m_k, m_alpha, theta);
            return *this;
        }

        size_t k() const
        {
            return m_k;
        }

      private:
        size_t m_k;
        const double *m_alpha;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &theta,
                            size_t k,
                            const double alpha[],
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(theta.derived())
    {
        rng r(k, alpha, type, seed);
        return fill(theta, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &theta, drch::rng &r)
        -> decltype(theta.derived())
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");
        eigen_assert((theta.size() % r.k()) == 0);

        double *data = theta.derived().data();
        for (Index i = 0; i < theta.size(); i += r.k()) {
            r.next(&data[i]);
        }
        return theta.derived();
    }

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(size_t k, const double alpha[])
            : m_k(k)
            , m_alpha(alpha)
        {
        }

        double pdf(const double theta[]) const
        {
            return gsl_ran_dirichlet_pdf(m_k, m_alpha, theta);
        }

        double lnpdf(const double theta[]) const
        {
            return gsl_ran_dirichlet_lnpdf(m_k, m_alpha, theta);
        }

      private:
        size_t m_k;
        const double *m_alpha;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &theta, size_t k, const double alpha[])
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");

        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(M2VNUM_ROW(T, theta),
                                       M2VNUM_COL(T, theta),
                                       pdf_functor<T>(theta.derived(),
                                                      k,
                                                      alpha));
    }

  private:
    drch() = delete;

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive_kpnum<T, double>::type;

        pdf_functor(const T &theta, size_t k, const double alpha[])
            : m_result(new double[M2V_NUM(T, theta)])
        {
            eigen_assert(k == M2V_DIM(T, theta));

            size_t n = M2V_NUM(T, theta);
            typename type_eval<T>::type m_theta(theta.eval());
            dist m_dist(k, alpha);
            for (size_t i = 0; i < n; ++i) {
                m_result.get()[i] = m_dist.pdf(&m_theta.data()[i * k]);
            }
        }

        double operator()(Index i) const
        {
            return m_result.get()[i];
        }

      private:
        std::shared_ptr<double> m_result;
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

#endif /* __IEXP_RAND_DIRICHLET__ */
