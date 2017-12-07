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

#ifndef __IEXP_FFT2__
#define __IEXP_FFT2__

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
inline void fft2_impl(const int n0, const int n1, const T *i, U *o)
{
    fftw3::get_plan(n0, n1, i, o, true).fwd(n0, n1, i, o);
}

template <typename T>
class fft2_functor
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

    fft2_functor(const T &x)
        : m_in_n1(x.cols())
        , m_out_n1(IS_COMPLEX(typename T::Scalar) ? m_in_n1
                                                  : ((m_in_n1 >> 1) + 1))
        , m_result(x.rows(), m_out_n1)
    {
        static_assert(T::Flags & RowMajorBit, "must be row major matrix");

        typename type_eval<T>::type m_x(x.eval());
        fft2_impl(x.rows(), m_in_n1, m_x.data(), m_result.data());
    }

    const Scalar operator()(Index i, Index j) const
    {
        if (j < m_out_n1) {
            return m_result(i, j);
        } else if (i == 0) {
            return std::conj(m_result(i, m_in_n1 - j));
        } else {
            return std::conj(m_result(m_in_n1 - i, m_in_n1 - j));
        }
    }

  private:
    const Index m_in_n1, m_out_n1;
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<fft2_functor<T>, typename fft2_functor<T>::ArrayType>
fft2(const ArrayBase<T> &x)
{
    using ArrayType = typename fft2_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  fft2_functor<T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_FFT2__ */
