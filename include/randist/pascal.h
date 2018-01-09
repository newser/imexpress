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

#ifndef __IEXP_RANDIST_PASCAL__
#define __IEXP_RANDIST_PASCAL__

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

class pascal
{
  public:
    pascal(const double p, const unsigned int n)
        : m_p(p)
        , m_n(n)
    {
    }

    double pdf(const unsigned int x) const
    {
        return gsl_ran_pascal_pdf(x, m_p, m_n);
    }

    double p(const unsigned int x) const
    {
        return gsl_cdf_pascal_P(x, m_p, m_n);
    }

    double q(const unsigned int x) const
    {
        return gsl_cdf_pascal_Q(x, m_p, m_n);
    }

  private:
    const double m_p;
    const unsigned int m_n;
};

template <typename T>
class pascal_pdf_functor
{
  public:
    using ArrayType = Array<double,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    pascal_pdf_functor(const T &x, const double p, const unsigned int n)
        : m_x(x)
        , m_pascal(p, n)
    {
    }

    double operator()(Index i, Index j) const
    {
        return m_pascal.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    pascal m_pascal;
};

template <typename T>
inline CwiseNullaryOp<pascal_pdf_functor<T>,
                      typename pascal_pdf_functor<T>::ArrayType>
pascal_pdf(const ArrayBase<T> &x, const double p, const unsigned int n)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    using ArrayType = typename pascal_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  pascal_pdf_functor<T>(x.derived(), p, n));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_PASCAL__ */
