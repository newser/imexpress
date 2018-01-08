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

#ifndef __IEXP_RANDIST_POISSON__
#define __IEXP_RANDIST_POISSON__

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

class poiss
{
  public:
    poiss(const double mu)
        : m_mu(mu)
    {
    }

    double pdf(const double x) const
    {
        return gsl_ran_poisson_pdf(x, m_mu);
    }

    double p(const double x) const
    {
        return gsl_cdf_poisson_P(x, m_mu);
    }

    double q(const double x) const
    {
        return gsl_cdf_poisson_Q(x, m_mu);
    }

  private:
    const double m_mu;
};

template <typename T>
class poiss_pdf_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    poiss_pdf_functor(const T &x, const typename T::Scalar mu)
        : m_x(x)
        , m_poiss(mu)
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_poiss.pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    poiss m_poiss;
};

template <typename T>
inline CwiseNullaryOp<poiss_pdf_functor<T>,
                      typename poiss_pdf_functor<T>::ArrayType>
poiss_pdf(const ArrayBase<T> &x, const typename T::Scalar mu)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "scalar can only be double");

    using ArrayType = typename poiss_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  poiss_pdf_functor<T>(x.derived(), mu));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_POISSON__ */
