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

#ifndef __IEXP_IFFT__
#define __IEXP_IFFT__

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
inline void ifft_impl(const int n, const T *i, U *o)
{
    fftw3::get_plan(n, i, o, false).inv(n, i, o);
}

template <bool normalize, typename T>
class ifft_functor
{
  public:
    using Scalar = typename TYPE_CHOOSE(IS_COMPLEX(typename T::Scalar),
                                        typename T::Scalar,
                                        std::complex<typename T::Scalar>);
    using ResultType = typename dense_derive<T, Scalar>::type;

    ifft_functor(const T &x)
        : m_result(new Scalar[x.size()])
    {
        typename type_eval<T>::type m_x(x.eval());
        ifft_impl((int)m_x.size(), m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            Index n = x.size();
            for (int i = 0; i < n; ++i) {
                p[i] /= n;
            }
        }
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool normalize = false, typename T = void>
inline CwiseNullaryOp<ifft_functor<normalize, T>,
                      typename ifft_functor<normalize, T>::ResultType>
ifft(const DenseBase<T> &x)
{
    // if x is real, x[i](i > 0) should be conjudate with x[size - i]

    using ResultType = typename ifft_functor<normalize, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   ifft_functor<normalize, T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_IFFT__ */
