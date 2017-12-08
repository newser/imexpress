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

#ifndef __IEXP_IDCT__
#define __IEXP_IDCT__

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

template <kind k, typename T>
inline void idct_impl(const int n, const T *i, T *o)
{
    fftw3::get_plan(n, i, o, false).template inv<k>(n, i, o);
}

template <kind k>
inline size_t idct_scale(const size_t n)
{
    return n << 1;
}

template <>
inline size_t idct_scale<DCT_I>(const size_t n)
{
    return (n - 1) << 1;
}

template <bool normalize, kind k, typename T>
class idct_functor
{
  public:
    using ArrayType =
        Array<typename T::Scalar, Dynamic, 1, ColMajor, Dynamic, 1>;

    idct_functor(const T &x)
        : m_result(x.size())
    {
        typename type_eval<T>::type m_x(x.eval());
        idct_impl<k>(m_x.size(), m_x.data(), m_result.data());

        if (normalize) {
            m_result /= idct_scale<k>(m_x.size());
        }
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result[i];
    }

  private:
    template <kind __k>
    static size_t ss(size_t n)
    {
        return 2 * n;
    }

    ArrayType m_result;
};

template <bool normalize = false, kind k = DCT_II, typename T = void>
inline CwiseNullaryOp<idct_functor<normalize, k, T>,
                      typename idct_functor<normalize, k, T>::ArrayType>
idct(const ArrayBase<T> &x)
{
    static_assert(IS_DCT(k), "not dct kind");
    eigen_assert(IS_VEC(x));

    using ArrayType = typename idct_functor<normalize, k, T>::ArrayType;
    return ArrayType::NullaryExpr(x.size(),
                                  idct_functor<normalize, k, T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_IDCT__ */