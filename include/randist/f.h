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
 * You should have received nu1 copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

#ifndef __IEXP_RANDIST_F__
#define __IEXP_RANDIST_F__

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

class f
{
  public:
    f(const double nu1, const double nu2)
        : m_nu1(nu1)
        , m_nu2(nu2)
    {
    }

    double pdf(const double x) const
    {
        return gsl_ran_fdist_pdf(x, m_nu1, m_nu2);
    }

    double p(const double x) const
    {
        return gsl_cdf_fdist_P(x, m_nu1, m_nu2);
    }

    double invp(const double x) const
    {
        return gsl_cdf_fdist_Pinv(x, m_nu1, m_nu2);
    }

    double q(const double x) const
    {
        return gsl_cdf_fdist_Q(x, m_nu1, m_nu2);
    }

    double invq(const double x) const
    {
        return gsl_cdf_fdist_Qinv(x, m_nu1, m_nu2);
    }

  private:
    const double m_nu1, m_nu2;
};

template <typename T>
class f_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    f_pdf_functor(const T &x,
                  const typename T::Scalar nu1,
                  const typename T::Scalar nu2)
        : m_x(x)
        , m_f(nu1, nu2)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_f.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    f m_f;
};

template <typename T>
inline CwiseNullaryOp<f_pdf_functor<T>, typename f_pdf_functor<T>::ArrayType>
f_pdf(const ArrayBase<T> &x,
      const typename T::Scalar nu1,
      const typename T::Scalar nu2)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename f_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  f_pdf_functor<T>(x.derived(), nu1, nu2));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_F__ */
