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

#ifndef __IEXP_IFFT2__
#define __IEXP_IFFT2__

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

template <typename T, typename U>
inline void ifft2_impl(const int n0, const int n1, const T *i, U *o)
{
    fftw3::get_plan(n0, n1, i, o, false).inv(n0, n1, i, o);
}

template <bool normalize, typename T>
class ifft2_functor
{
  public:
    using Scalar = typename TYPE_CHOOSE(IS_COMPLEX(typename T::Scalar),
                                        typename T::Scalar,
                                        std::complex<typename T::Scalar>);
    using ArrayType = Array<Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            RowMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    ifft2_functor(const T &x)
        : m_result(x.rows(), x.cols())
    {
        static_assert(T::Flags & RowMajorBit, "must be row major matrix");

        typename type_eval<T>::type m_x(x.eval());
        ifft2_impl((int)x.rows(), (int)x.cols(), m_x.data(), m_result.data());

        if (normalize) {
            m_result /= (T::Scalar)m_x.size();
        }
    }

    const Scalar &operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <bool normalize = false, typename T = void>
inline CwiseNullaryOp<ifft2_functor<normalize, T>,
                      typename ifft2_functor<normalize, T>::ArrayType>
ifft2(const ArrayBase<T> &x)
{
    using ArrayType = typename ifft2_functor<normalize, T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  ifft2_functor<normalize, T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_IFFT2__ */
