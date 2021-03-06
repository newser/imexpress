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
inline void fft_impl(int n, const T *i, U *o)
{
    fftw3::get_plan(n, i, o, true).fwd(n, i, o);
}

template <typename T>
class fft_functor
{
  public:
    using Scalar = typename TYPE_CHOOSE(IS_COMPLEX(typename T::Scalar),
                                        typename T::Scalar,
                                        std::complex<typename T::Scalar>);
    using ResultType = typename dense_derive<T, Scalar>::type;

    fft_functor(const T &x)
        : m_in_size(x.size())
        , m_out_size(IS_COMPLEX(typename T::Scalar) ? m_in_size
                                                    : ((m_in_size >> 1) + 1))
        , m_result(new Scalar[m_out_size])
    {
        typename type_eval<T>::type m_x(x.eval());
        fft_impl((int)m_in_size, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i) const
    {
        if (i < m_out_size) {
            return m_result.get()[i];
        } else {
            return std::conj(m_result.get()[m_in_size - i]);
        }
    }

  private:
    Index m_in_size, m_out_size;
    std::shared_ptr<Scalar> m_result;
};

template <typename T>
inline CwiseNullaryOp<fft_functor<T>, typename fft_functor<T>::ResultType> fft(
    const DenseBase<T> &x)
{
    using ResultType = typename fft_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   fft_functor<T>(x.derived()));
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
