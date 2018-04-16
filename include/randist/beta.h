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

#ifndef __IEXP_RANDIST_BETA__
#define __IEXP_RANDIST_BETA__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rdist {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class beta
{
  public:
    beta(double a, double b)
        : m_a(a)
        , m_b(b)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_beta_pdf(x, m_a, m_b);
    }

    double p(double x) const
    {
        return gsl_cdf_beta_P(x, m_a, m_b);
    }

    double invp(double x) const
    {
        return gsl_cdf_beta_Pinv(x, m_a, m_b);
    }

    double q(double x) const
    {
        return gsl_cdf_beta_Q(x, m_a, m_b);
    }

    double invq(double x) const
    {
        return gsl_cdf_beta_Qinv(x, m_a, m_b);
    }

  private:
    double m_a, m_b;
};

template <typename T>
class beta_pdf_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    beta_pdf_functor(const T &x, Scalar a, Scalar b)
        : m_x(x)
        , m_beta(a, b)
    {
    }

    Scalar operator()(Index i, Index j) const
    {
        return m_beta.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    beta m_beta;
};

template <typename T>
inline CwiseNullaryOp<beta_pdf_functor<T>,
                      typename beta_pdf_functor<T>::ResultType>
beta_pdf(const ArrayBase<T> &x, typename T::Scalar a, typename T::Scalar b)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ResultType = typename beta_pdf_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   beta_pdf_functor<T>(x.derived(), a, b));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_BETA__ */
