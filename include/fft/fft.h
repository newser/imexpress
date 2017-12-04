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

#ifndef __IEXP_FFT__
#define __IEXP_FFT__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>
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
inline void fft_impl(const int n, const T *i, U *o)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void fft_impl(const int n, const double *i, std::complex<double> *o)
{
    // fftw3::fwd
}

template <typename T>
class fft_functor
{
  public:
    typedef Array<typename std::complex<typename T::Scalar>,
                  Dynamic,
                  1,
                  ColMajor,
                  Dynamic,
                  1>
        ArrayType;

    fft_functor(const T &x)
        : m_in_size(x.size())
        , m_out_size(IS_COMPLEX(typename T::Scalar) ? m_in_size
                                                    : ((m_in_size >> 1) + 1))
        , m_result(m_out_size)
    {
    }

    const typename T::Scalar operator()(Index i) const
    {
        if (i < m_out_size) {
            return m_result[i];
        } else {
            return std::conj(m_result[m_in_size - i]);
        }
    }

  private:
    const Index m_in_size, m_out_size;
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<fft_functor<T>, typename fft_functor<T>::ArrayType> fft(
    const ArrayBase<T> &x)
{
    eigen_assert(IS_VEC(x));

    typedef typename fft_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.size(), fft_functor<T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_FFT__ */
