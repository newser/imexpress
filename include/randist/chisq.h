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

#ifndef __IEXP_RANDIST_CHISQ__
#define __IEXP_RANDIST_CHISQ__

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

class chisq
{
  public:
    chisq(double nu)
        : m_mu(nu)
    {
    }

    double pdf(double x) const
    {
        return gsl_ran_chisq_pdf(x, m_mu);
    }

    double p(double x) const
    {
        return gsl_cdf_chisq_P(x, m_mu);
    }

    double invp(double x) const
    {
        return gsl_cdf_chisq_Pinv(x, m_mu);
    }

    double q(double x) const
    {
        return gsl_cdf_chisq_Q(x, m_mu);
    }

    double invq(double x) const
    {
        return gsl_cdf_chisq_Qinv(x, m_mu);
    }

  private:
    double m_mu;
};

template <typename T>
class chisq_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    chisq_pdf_functor(const T &x, typename T::Scalar nu)
        : m_x(x)
        , m_chisq(nu)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_chisq.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    chisq m_chisq;
};

template <typename T>
inline CwiseNullaryOp<chisq_pdf_functor<T>,
                      typename chisq_pdf_functor<T>::ArrayType>
chisq_pdf(const ArrayBase<T> &x, typename T::Scalar nu)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename chisq_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  chisq_pdf_functor<T>(x.derived(), nu));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_CHISQ__ */
