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

#ifndef __IEXP_RAND_MUL_NOMIAL__
#define __IEXP_RAND_MUL_NOMIAL__

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

class mnom
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
            const double p[],
            unsigned int N,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_k(k)
            , m_p(p)
            , m_N(N)
            , m_rng(type, seed)
        {
        }

        rng &seed(unsigned long seed)
        {
            m_rng.seed(seed);
            return *this;
        }

        rng &next(unsigned int n[])
        {
            gsl_ran_multinomial(m_rng.gsl(), m_k, m_N, m_p, n);
            return *this;
        }

        size_t k()
        {
            return m_k;
        }

      private:
        size_t m_k;
        const double *m_p;
        unsigned int m_N;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            size_t k,
                            const double p[],
                            unsigned int N,
                            rand::rng::type type = DEFAULT_RNG_TYPE,
                            unsigned long seed = 0) -> decltype(x.derived())
    {
        mnom::rng r(k, p, N, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mnom::rng &r)
        -> decltype(x.derived())
    {
        static_assert(IS_INTEGER(typename T::Scalar),
                      "only support integer scalar");

        return fill(x, r, TYPE_BOOL(TP4(T) == RowMajor)());
    }

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(size_t k, const double p[])
            : m_k(k)
            , m_p(p)
        {
        }

        double pdf(const unsigned int x[]) const
        {
            return gsl_ran_multinomial_pdf(m_k, m_p, x);
        }

        double lnpdf(const unsigned int x[]) const
        {
            return gsl_ran_multinomial_lnpdf(m_k, m_p, x);
        }

        size_t k()
        {
            return m_k;
        }

      private:
        size_t m_k;
        const double *m_p;
    };

    template <typename T>
    static inline CwiseNullaryOp<pdf_functor<T>,
                                 typename pdf_functor<T>::ResultType>
    pdf(const DenseBase<T> &x, size_t k, const double p[])
    {
        static_assert(IS_INTEGER(typename T::Scalar),
                      "only support integer scalar");

        using ResultType = typename pdf_functor<T>::ResultType;
        return ResultType::NullaryExpr(x.IsRowMajor ? x.rows() : 1,
                                       x.IsRowMajor ? 1 : x.cols(),
                                       pdf_functor<T>(x.derived(), k, p));
    }

  private:
    mnom() = delete;

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mnom::rng &r, std::true_type)
        -> decltype(x.derived())
    {
        // row major
        size_t k = r.k();
        eigen_assert(x.cols() == k);

        typename T::Scalar *data = x.derived().data();
        for (Index i = 0; i < x.rows(); ++i) {
            r.next(&data[i * k]);
        }
        return x.derived();
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mnom::rng &r, std::false_type)
        -> decltype(x.derived())
    {
        // column major
        size_t k = r.k();
        eigen_assert(x.rows() == k);

        typename T::Scalar *data = x.derived().data();
        for (Index i = 0; i < x.cols(); ++i) {
            r.next(&data[i * k]);
        }
        return x.derived();
    }

    template <typename T>
    class pdf_functor
    {
      public:
        using ResultType = typename dense_derive_kpnum<T, double>::type;

        pdf_functor(const T &x, size_t k, const double p[])
            : m_result(new double[M2V_NUM(T, x)])
        {
            typename type_eval<T>::type m_x(x.eval());
            dist m_dist(k, p);
            compute(m_x, m_dist, TYPE_BOOL(TP4(T) == RowMajor)());
        }

        double operator()(Index i) const
        {
            return m_result.get()[i];
        }

      private:
        std::shared_ptr<double> m_result;

        void compute(typename type_eval<T>::type &m_x,
                     dist &m_dist,
                     std::true_type)
        {
            // row major
            size_t k = m_dist.k();
            eigen_assert(m_x.cols() == k);
            for (Index i = 0; i < m_x.rows(); ++i) {
                m_result.get()[i] = m_dist.pdf(&m_x.data()[i * k]);
            }
        }

        void compute(typename type_eval<T>::type &m_x,
                     dist &m_dist,
                     std::false_type)
        {
            // column major
            size_t k = m_dist.k();
            eigen_assert(m_x.rows() == k);
            for (Index i = 0; i < m_x.cols(); ++i) {
                m_result.get()[i] = m_dist.pdf(&m_x.data()[i * k]);
            }
        }
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

#endif /* __IEXP_RAND_MUL_NOMIAL__ */
