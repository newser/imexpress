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

#ifndef __IEXP_DCT2__
#define __IEXP_DCT2__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <fft/fftw/plan_cache.h>
#include <fft/fftw/plan_double.h>
#include <fft/fftw/plan_single.h>

IEXP_NS_BEGIN

namespace fft {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <kind k0, kind k1, typename T>
inline void dct2_impl(int n0, int n1, const T *i, T *o)
{
    fftw3::get_plan<k0, k1>(n0, n1, i, o, true)
        .template fwd<k0, k1>(n0, n1, i, o);
}

template <kind k0,
          kind k1,
          typename T,
          bool is_row_major = bool(TP4(T) == RowMajor)>
class dct2_functor;

template <kind k0, kind k1, typename T>
class dct2_functor<k0, k1, T, true>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    dct2_functor(const T &x)
        : m_n0(x.rows())
        , m_n1(x.cols())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        dct2_impl<k0, k1>((int)m_n0, (int)m_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define DCT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return DCT2_RESULT(i, j);
#undef DCT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

template <kind k0, kind k1, typename T>
class dct2_functor<k0, k1, T, false>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    dct2_functor(const T &x)
        : m_n0(x.cols())
        , m_n1(x.rows())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        dct2_impl<k0, k1>((int)m_n0, (int)m_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define DCT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return DCT2_RESULT(j, i);
#undef DCT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

template <kind k0 = DCT_II, kind k1 = DCT_II, typename T = void>
inline CwiseNullaryOp<dct2_functor<k0, k1, T>,
                      typename dct2_functor<k0, k1, T>::ResultType>
dct2(const DenseBase<T> &x)
{
    static_assert(IS_DCT(k0) && IS_DCT(k1), "not dct kind");

    using ResultType = typename dct2_functor<k0, k1, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dct2_functor<k0, k1, T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DCT2__ */
