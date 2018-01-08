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

#ifndef __IEXP_RANDIST_WEIBULL__
#define __IEXP_RANDIST_WEIBULL__

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

class wbl
{
  public:
    wbl(const double a, const double b)
        : m_a(a)
        , m_b(b)
    {
    }

    double pdf(const double x) const
    {
        return gsl_ran_weibull_pdf(x, m_a, m_b);
    }

    double p(const double x) const
    {
        return gsl_cdf_weibull_P(x, m_a, m_b);
    }

    double invp(const double x) const
    {
        return gsl_cdf_weibull_Pinv(x, m_a, m_b);
    }

    double q(const double x) const
    {
        return gsl_cdf_weibull_Q(x, m_a, m_b);
    }

    double invq(const double x) const
    {
        return gsl_cdf_weibull_Qinv(x, m_a, m_b);
    }

  private:
    const double m_a, m_b;
};

template <typename T>
class wbl_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    wbl_pdf_functor(const T &x,
                    const typename T::Scalar a,
                    const typename T::Scalar b)
        : m_x(x)
        , m_wbl(a, b)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_wbl.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    wbl m_wbl;
};

template <typename T>
inline CwiseNullaryOp<wbl_pdf_functor<T>,
                      typename wbl_pdf_functor<T>::ArrayType>
wbl_pdf(const ArrayBase<T> &x,
        const typename T::Scalar a,
        const typename T::Scalar b)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename wbl_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  wbl_pdf_functor<T>(x.derived(), a, b));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_WEIBULL__ */
