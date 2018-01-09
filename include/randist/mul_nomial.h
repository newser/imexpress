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

#ifndef __IEXP_RANDIST_MUL_NOMIAL__
#define __IEXP_RANDIST_MUL_NOMIAL__

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

class mnom
{
  public:
    mnom(size_t K, const double p[])
        : m_k(K)
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

  private:
    size_t m_k;
    const double *m_p;
};

#if 0
template <typename T>
class mnom_pdf_functor
{
  public:
    using ArrayType = Array<double,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    mnom_pdf_functor(const T &x, size_t K,
                     const double p[])
        : m_x(x)
        , m_mnom(K, p)
    {
    }

    double operator()(Index i, Index j) const
    {
        return m_mnom.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    mnom m_mnom;
};

template <typename T>
inline CwiseNullaryOp<mnom_pdf_functor<T>,
                      typename mnom_pdf_functor<T>::ArrayType>
    mnom_pdf(const ArrayBase<T> &x, size_t K,
             const double p[])
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                  TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    using ArrayType = typename mnom_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  mnom_pdf_functor<T>(x.derived(), K, p));
}
#endif

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_MUL_NOMIAL__ */
