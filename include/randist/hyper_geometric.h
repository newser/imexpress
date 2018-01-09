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

#ifndef __IEXP_RANDIST_HYPER_GEOMETRIC__
#define __IEXP_RANDIST_HYPER_GEOMETRIC__

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

class hgeo
{
  public:
    hgeo(unsigned int n1, unsigned int n2, unsigned int t)
        : m_n1(n1)
        , m_n2(n2)
        , m_t(t)
    {
    }

    double pdf(unsigned int x) const
    {
        return gsl_ran_hypergeometric_pdf(x, m_n1, m_n2, m_t);
    }

    double p(unsigned int x) const
    {
        return gsl_cdf_hypergeometric_P(x, m_n1, m_n2, m_t);
    }

    double q(unsigned int x) const
    {
        return gsl_cdf_hypergeometric_Q(x, m_n1, m_n2, m_t);
    }

  private:
    unsigned int m_n1, m_n2, m_t;
};

template <typename T>
class hgeo_pdf_functor
{
  public:
    using ArrayType = Array<double,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    hgeo_pdf_functor(const T &x,
                     unsigned int n1,
                     unsigned int n2,
                     unsigned int t)
        : m_x(x)
        , m_hgeo(n1, n2, t)
    {
    }

    double operator()(Index i, Index j) const
    {
        return m_hgeo.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    hgeo m_hgeo;
};

template <typename T>
inline CwiseNullaryOp<hgeo_pdf_functor<T>,
                      typename hgeo_pdf_functor<T>::ArrayType>
hgeo_pdf(const ArrayBase<T> &x,
         unsigned int n1,
         unsigned int n2,
         unsigned int t)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    using ArrayType = typename hgeo_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  hgeo_pdf_functor<T>(x.derived(), n1, n2, t));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_HYPER_GEOMETRIC__ */
