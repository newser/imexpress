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

#ifndef __IEXP_RANDIST_GAUSS__
#define __IEXP_RANDIST_GAUSS__

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

class gauss
{
  public:
    gauss(double sigma = 1.0)
        : m_sigma(sigma)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_gaussian_pdf(x, m_sigma);
    }

    double p(double x) const
    {
        return gsl_cdf_gaussian_P(x, m_sigma);
    }

    double invp(double x) const
    {
        return gsl_cdf_gaussian_Pinv(x, m_sigma);
    }

    double q(double x) const
    {
        return gsl_cdf_gaussian_Q(x, m_sigma);
    }

    double invq(double x) const
    {
        return gsl_cdf_gaussian_Qinv(x, m_sigma);
    }

  private:
    double m_sigma;
};

template <typename T>
class gauss_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    gauss_pdf_functor(const T &x, typename T::Scalar sigma)
        : m_x(x)
        , m_gauss(sigma)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_gauss.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    gauss m_gauss;
};

template <typename T>
inline CwiseNullaryOp<gauss_pdf_functor<T>,
                      typename gauss_pdf_functor<T>::ArrayType>
gauss_pdf(const ArrayBase<T> &x, typename T::Scalar sigma)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename gauss_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  gauss_pdf_functor<T>(x.derived(), sigma));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_GAUSS__ */
