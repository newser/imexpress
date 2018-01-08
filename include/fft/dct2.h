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
inline void dct2_impl(const int n0, const int n1, const T *i, T *o)
{
    fftw3::get_plan<k0, k1>(n0, n1, i, o, true)
        .template fwd<k0, k1>(n0, n1, i, o);
}

template <kind k0, kind k1, typename T>
class dct2_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            RowMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    dct2_functor(const T &x)
        : m_result(x.rows(), x.cols())
    {
        static_assert(T::Flags & RowMajorBit, "must be row major matrix");

        typename type_eval<T>::type m_x(x.eval());
        dct2_impl<k0, k1>((int)m_x.rows(), (int)m_x.cols(), m_x.data(), m_result.data());
    }

    const typename T::Scalar &operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <kind k0 = DCT_II, kind k1 = DCT_II, typename T = void>
inline CwiseNullaryOp<dct2_functor<k0, k1, T>,
                      typename dct2_functor<k0, k1, T>::ArrayType>
dct2(const ArrayBase<T> &x)
{
    static_assert(IS_DCT(k0) && IS_DCT(k1), "not dct kind");

    using ArrayType = typename dct2_functor<k0, k1, T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
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
