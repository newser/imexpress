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

#ifndef __IEXP_RANDIST_DISCRETE__
#define __IEXP_RANDIST_DISCRETE__

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

class discrete
{
  public:
    discrete(const size_t K, const double P[])
        : m_g(nullptr)
    {
        m_g = gsl_ran_discrete_preproc(K, P);
        IEXP_NOT_NULLPTR(m_g);
    }

    ~discrete()
    {
        gsl_ran_discrete_free(m_g);
    }

    double pdf(const size_t x) const
    {
        return gsl_ran_discrete_pdf(x, m_g);
    }

  private:
    discrete(const discrete &) = delete;
    discrete &operator=(const discrete &other) = delete;

    gsl_ran_discrete_t *m_g;
};

template <typename T>
class discrete_pdf_functor
{
  public:
    using ArrayType = Array<double,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    discrete_pdf_functor(const T &x, const size_t K, const double P[])
        : m_x(x)
        , m_discrete(new discrete(K, P))
    {
    }

    typename T::Scalar operator()(Index i, Index j) const
    {
        return m_discrete->pdf(m_x(i, j));
    }

  private:
    const T &m_x;
    std::shared_ptr<discrete> m_discrete;
};

template <typename T>
inline CwiseNullaryOp<discrete_pdf_functor<T>,
                      typename discrete_pdf_functor<T>::ArrayType>
discrete_pdf(const ArrayBase<T> &x, const size_t K, const double P[])
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "scalar can only be int or unsigned int");

    using ArrayType = typename discrete_pdf_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  discrete_pdf_functor<T>(x.derived(), K, P));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RANDIST_DISCRETE__ */
