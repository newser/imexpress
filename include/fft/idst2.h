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

#ifndef __IEXP_IDST2__
#define __IEXP_IDST2__

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
inline void idst2_impl(const int n0, const int n1, const T *i, T *o)
{
    fftw3::get_plan<k0, k1>(n0, n1, i, o, false)
        .template inv<k0, k1>(n0, n1, i, o);
}

template <kind k0, kind k1>
inline size_t idst2_scale(const size_t rows, const size_t cols)
{
    const size_t r = k0 == DST_I ? rows + 1 : rows;
    const size_t c = k1 == DST_I ? cols + 1 : cols;
    return (r * c) << 2;
}

template <bool normalize, kind k0, kind k1, typename T>
class idst2_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            RowMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    idst2_functor(const T &x)
        : m_result(x.rows(), x.cols())
    {
        static_assert(T::Flags & RowMajorBit, "must be row major matrix");

        typename type_eval<T>::type m_x(x.eval());
        idst2_impl<k0, k1>(m_x.rows(), m_x.cols(), m_x.data(), m_result.data());

        if (normalize) {
            m_result /= idst2_scale<k0, k1>(m_x.rows(), m_x.cols());
        }
    }

    const typename T::Scalar &operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <bool normalize = false,
          kind k0 = DST_I,
          kind k1 = DST_I,
          typename T = void>
inline CwiseNullaryOp<idst2_functor<normalize, k0, k1, T>,
                      typename idst2_functor<normalize, k0, k1, T>::ArrayType>
idst2(const ArrayBase<T> &x)
{
    static_assert(IS_DST(k0) && IS_DST(k1), "not dst kind");

    using ArrayType = typename idst2_functor<normalize, k0, k1, T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  idst2_functor<normalize, k0, k1, T>(
                                      x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_IDST2__ */
