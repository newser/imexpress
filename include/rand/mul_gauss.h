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

#ifndef __IEXP_RAND_MUL_GAUSS__
#define __IEXP_RAND_MUL_GAUSS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class mgauss
{
    template <typename T>
    class mean_functor;
    template <typename T>
    class cov_functor;
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
            const double mu[],
            const double cov[],
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : m_mu_block{k, const_cast<double *>(mu)}
            , m_mu{k, 1, const_cast<double *>(mu), &m_mu_block, 0}
            , m_rng(type, seed)
        {
            // we have to alloc matrix here, because the cov matrix would be
            // overwritten by cholesky decomposition result
            m_L = gsl_matrix_alloc(k, k);
            IEXP_NOT_NULLPTR(m_L);
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    gsl_matrix_set(m_L, i, j, cov[i * k + j]);
                }
            }
            gsl_linalg_cholesky_decomp1(m_L);
        }

        ~rng()
        {
            gsl_matrix_free(m_L);
        }

        rng &seed(unsigned long seed)
        {
            m_rng.seed(seed);
            return *this;
        }

        rng &next(double x[])
        {
            gsl_block b{m_mu_block.size, x};
            gsl_vector result{m_mu_block.size, 1, x, &b, 0};
            gsl_ran_multivariate_gaussian(m_rng.gsl(), &m_mu, m_L, &result);
            return *this;
        }

        size_t k()
        {
            return m_mu_block.size;
        }

      private:
        // prevent from sharing m_L
        rng(const rng &) = delete;
        rng &operator=(const rng &other) = delete;

        gsl_matrix *m_L;
        gsl_block m_mu_block;
        gsl_vector m_mu;
        rand::rng m_rng;
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            size_t k,
                            const double mu[],
                            const double cov[],
                            rand::rng::type type = DEFAULT_RNG_TYPE,
                            unsigned long seed = 0) -> decltype(x.derived())
    {
        mgauss::rng r(k, mu, cov, type, seed);
        return fill(x, r);
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mgauss::rng &r)
        -> decltype(x.derived())
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");

        return fill(x, r, TYPE_BOOL(TP4(T) == RowMajor)());
    }

    // ========================================
    // distribution
    // ========================================

    class dist
    {
      public:
        dist(size_t k, double mu[], double cov[])
            : m_mu_block{k, mu}
            , m_mu{k, 1, mu, &m_mu_block, 0}
        {
            m_work = gsl_vector_alloc(k);
            IEXP_NOT_NULLPTR(m_work);

            m_L = gsl_matrix_alloc(k, k);
            IEXP_NOT_NULLPTR(m_L);
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    gsl_matrix_set(m_L, i, j, cov[i * k + j]);
                }
            }
            gsl_linalg_cholesky_decomp1(m_L);
        }

        ~dist()
        {
            gsl_vector_free(m_work);
            gsl_matrix_free(m_L);
        }

        double pdf(const double x[]) const
        {
            gsl_block xb{m_mu_block.size, (double *)x};
            gsl_vector xv{m_mu_block.size, 1, (double *)x, &xb, 0};

            // due to m_work, it's not reentrant
            double result;
            gsl_ran_multivariate_gaussian_pdf(&xv, &m_mu, m_L, &result, m_work);
            return result;
        }

        double lnpdf(const double x[]) const
        {
            gsl_block xb{m_mu_block.size, (double *)x};
            gsl_vector xv{m_mu_block.size, 1, (double *)x, &xb, 0};

            // due to m_work, it's not reentrant
            double result;
            gsl_ran_multivariate_gaussian_log_pdf(&xv,
                                                  &m_mu,
                                                  m_L,
                                                  &result,
                                                  m_work);
            return result;
        }

        // n sample, each has k dim
        static void mean(size_t n, size_t k, const double *x, double *mu)
        {
            gsl_block xb{n * k, (double *)x};
            gsl_matrix xm{n, k, k, (double *)x, &xb, 0};

            gsl_block mb{k, mu};
            gsl_vector mv{k, 1, mu, &mb, 0};

            gsl_ran_multivariate_gaussian_mean(&xm, &mv);
        }

        // n sample, each has k dim
        static void cov(size_t n, size_t k, const double *x, double *cov)
        {
            gsl_block xb{n * k, (double *)x};
            gsl_matrix xm{n, k, k, (double *)x, &xb, 0};

            gsl_block cb{k * k, cov};
            gsl_matrix cm{k, k, k, cov, &cb, 0};

            gsl_ran_multivariate_gaussian_vcov(&xm, &cm);
        }

      private:
        dist(const dist &) = delete;
        dist &operator=(const dist &other) = delete;

        gsl_vector *m_work;
        gsl_matrix *m_L;
        gsl_block m_mu_block;
        gsl_vector m_mu;
    };

    // mean
    template <typename T>
    static inline CwiseNullaryOp<mean_functor<T>,
                                 typename mean_functor<T>::ResultType>
    mean(const DenseBase<T> &x)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");

        using ResultType = typename mean_functor<T>::ResultType;
        return ResultType::NullaryExpr((TP4(T) == RowMajor) ? 1 : x.rows(),
                                       (TP4(T) == RowMajor) ? x.cols() : 1,
                                       mean_functor<T>(x.derived()));
    }

    // cov
    template <typename T>
    static inline CwiseNullaryOp<cov_functor<T>,
                                 typename cov_functor<T>::ResultType>
    cov(const DenseBase<T> &x)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");

        using ResultType = typename cov_functor<T>::ResultType;
        return ResultType::NullaryExpr((TP4(T) == RowMajor) ? x.cols()
                                                            : x.rows(),
                                       (TP4(T) == RowMajor) ? x.cols()
                                                            : x.rows(),
                                       cov_functor<T>(x.derived()));
    }

  private:
    mgauss() = delete;

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mgauss::rng &r, std::true_type)
        -> decltype(x.derived())
    {
        // row major
        size_t k = r.k();
        eigen_assert(x.cols() == k);

        double *data = x.derived().data();
        for (Index i = 0; i < x.rows(); ++i) {
            r.next(&data[i * k]);
        }
        return x.derived();
    }

    template <typename T>
    static inline auto fill(DenseBase<T> &x, mgauss::rng &r, std::false_type)
        -> decltype(x.derived())
    {
        // colume major
        size_t k = r.k();
        eigen_assert(x.rows() == k);

        double *data = x.derived().data();
        for (Index i = 0; i < x.cols(); ++i) {
            r.next(&data[i * k]);
        }
        return x.derived();
    }

    template <typename T>
    class mean_functor
    {
      public:
        using Scalar = typename T::Scalar;
        using ResultType = typename dense_derive_kpdim<T>::type;

        mean_functor(const T &x)
            : m_result(new Scalar[M2V_DIM(T, x)])
        {
            typename type_eval<T>::type m_x(x.eval());
            dist::mean(M2V_NUM(T, x),
                       M2V_DIM(T, x),
                       m_x.data(),
                       m_result.get());
        }

        Scalar operator()(Index i) const
        {
            return m_result.get()[i];
        }

      private:
        std::shared_ptr<Scalar> m_result;
    };

    template <typename T>
    class cov_functor
    {
      public:
        using Scalar = typename T::Scalar;
        using ResultType = typename dense_derive<
            T,
            Scalar,
            (TP4(T) == RowMajor) ? T::ColsAtCompileTime : T::RowsAtCompileTime,
            (TP4(T) == RowMajor) ? T::ColsAtCompileTime : T::RowsAtCompileTime,
            (TP4(T) == RowMajor) ? RowMajor : ColMajor,
            (TP4(T) == RowMajor) ? T::MaxColsAtCompileTime
                                 : T::MaxRowsAtCompileTime,
            (TP4(T) == RowMajor) ? T::MaxColsAtCompileTime
                                 : T::MaxRowsAtCompileTime>::type;

        cov_functor(const T &x)
            : m_result((TP4(T) == RowMajor) ? x.cols() : x.rows(),
                       (TP4(T) == RowMajor) ? x.cols() : x.rows())
        {
            typename type_eval<T>::type m_x(x.eval());
            dist::cov((TP4(T) == RowMajor) ? x.rows() : x.cols(),
                      (TP4(T) == RowMajor) ? x.cols() : x.rows(),
                      m_x.data(),
                      m_result.data());
        }

        Scalar operator()(Index i, Index j) const
        {
            return m_result(i, j);
        }

      private:
        buf_rc<Scalar, TP4(T) == RowMajor> m_result;
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

#endif /* __IEXP_RAND_MUL_GAUSS__ */
