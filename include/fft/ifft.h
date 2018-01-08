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
    using ArrayType = Array<Scalar, Dynamic, 1, ColMajor, Dynamic, 1>;

    ifft_functor(const T &x)
        : m_result(x.size())
    {
        typename type_eval<T>::type m_x(x.eval());
        ifft_impl((int)m_x.size(), m_x.data(), m_result.data());

        if (normalize) {
            m_result /= (T::Scalar)m_x.size();
        }
    }

    const Scalar &operator()(Index i) const
    {
        return m_result[i];
    }

  private:
    ArrayType m_result;
};

template <bool normalize = false, typename T = void>
inline CwiseNullaryOp<ifft_functor<normalize, T>,
                      typename ifft_functor<normalize, T>::ArrayType>
ifft(const ArrayBase<T> &x)
{
    eigen_assert(IS_VEC(x));
    // if x is real, x[i](i > 0) should be conjudate with x[size - i]

    using ArrayType = typename ifft_functor<normalize, T>::ArrayType;
    return ArrayType::NullaryExpr(x.size(),
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
