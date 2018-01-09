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

#ifndef __IEXP_RANDIST_MUL_GAUSS__
#define __IEXP_RANDIST_MUL_GAUSS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rdist {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class mgauss
{
  public:
    mgauss(size_t k, double mu[], double cov[])
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

    ~mgauss()
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
        if (gsl_ran_multivariate_gaussian_pdf(&xv,
                                              &m_mu,
                                              m_L,
                                              &result,
                                              m_work) != GSL_SUCCESS) {
            RETURN_NAN_OR_THROW(std::runtime_error("mgauss pdf"));
        }
        return result;
    }

    double lnpdf(const double x[]) const
    {
        gsl_block xb{m_mu_block.size, (double *)x};
        gsl_vector xv{m_mu_block.size, 1, (double *)x, &xb, 0};

        // due to m_work, it's not reentrant
        double result;
        if (gsl_ran_multivariate_gaussian_log_pdf(&xv,
                                                  &m_mu,
                                                  m_L,
                                                  &result,
                                                  m_work) != GSL_SUCCESS) {
            RETURN_NAN_OR_THROW(std::runtime_error("mgauss pdf"));
        }
        return result;
    }

    // n sample, each has k dim
    static void mean(const size_t n,
                     const size_t k,
                     const double *x,
                     double *mu)
    {
        gsl_block xb{n * k, (double *)x};
        gsl_matrix xm{n, k, k, (double *)x, &xb, 0};

        gsl_block mb{k, mu};
        gsl_vector mv{k, 1, mu, &mb, 0};

        if (gsl_ran_multivariate_gaussian_mean(&xm, &mv) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("mgauss mean"));
        }
    }

    // n sample, each has k dim
    static void cov(const size_t n,
                    const size_t k,
                    const double *x,
                    double *cov)
    {
        gsl_block xb{n * k, (double *)x};
        gsl_matrix xm{n, k, k, (double *)x, &xb, 0};

        gsl_block cb{k * k, cov};
        gsl_matrix cm{k, k, k, cov, &cb, 0};

        if (gsl_ran_multivariate_gaussian_vcov(&xm, &cm) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("mgauss cov"));
        }
    }

  private:
    mgauss(const mgauss &) = delete;
    mgauss &operator=(const mgauss &other) = delete;

    gsl_vector *m_work;
    gsl_matrix *m_L;
    gsl_block m_mu_block;
    gsl_vector m_mu;
};

// ========================================
// mean
// ========================================

template <typename T>
class mgauss_mean_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::SizeAtCompileTime,
                            1,
                            ColMajor,
                            T::SizeAtCompileTime,
                            1>;

    mgauss_mean_functor(const T &x)
        : m_result(x.IsRowMajor ? x.cols() : x.rows())
    {
        typename type_eval<T>::type m_x(x.eval());
        mgauss::mean(x.IsRowMajor ? x.rows() : x.cols(),
                     x.IsRowMajor ? x.cols() : x.rows(),
                     m_x.data(),
                     m_result.data());
    }

    const typename T::Scalar &operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<mgauss_mean_functor<T>,
                      typename mgauss_mean_functor<T>::ArrayType>
mgauss_mean(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename mgauss_mean_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.IsRowMajor ? x.cols() : x.rows(),
                                  1,
                                  mgauss_mean_functor<T>(x.derived()));
}

// ========================================
// cov
// ========================================

template <typename T>
class mgauss_cov_functor
{
  public:
    using ArrayType = Array<
        typename T::Scalar,
        T::Flags & RowMajorBit ? T::ColsAtCompileTime : T::RowsAtCompileTime,
        T::Flags & RowMajorBit ? T::ColsAtCompileTime : T::RowsAtCompileTime,
        T::Flags & RowMajorBit ? RowMajor : ColMajor,
        T::Flags & RowMajorBit ? T::MaxColsAtCompileTime
                               : T::MaxRowsAtCompileTime,
        T::Flags & RowMajorBit ? T::MaxColsAtCompileTime
                               : T::MaxRowsAtCompileTime>;

    mgauss_cov_functor(const T &x)
        : m_result(x.IsRowMajor ? x.cols() : x.rows(),
                   x.IsRowMajor ? x.cols() : x.rows())
    {
        typename type_eval<T>::type m_x(x.eval());
        mgauss::cov(x.IsRowMajor ? x.rows() : x.cols(),
                    x.IsRowMajor ? x.cols() : x.rows(),
                    m_x.data(),
                    m_result.data());
    }

    const typename T::Scalar &operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<mgauss_cov_functor<T>,
                      typename mgauss_cov_functor<T>::ArrayType>
mgauss_cov(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename mgauss_cov_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.IsRowMajor ? x.cols() : x.rows(),
                                  x.IsRowMajor ? x.cols() : x.rows(),
                                  mgauss_cov_functor<T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_MUL_GAUSS__ */
