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

#ifndef __IEXP_IDST__
#define __IEXP_IDST__

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
inline void idst_impl(int n, const T *i, T *o)
{
    fftw3::get_plan<k>(n, i, o, false).template inv<k>(n, i, o);
}

template <kind k>
inline size_t idst_scale(size_t n)
{
    return n << 1;
}

template <>
inline size_t idst_scale<DST_I>(size_t n)
{
    return (n + 1) << 1;
}

template <bool normalize, kind k, typename T>
class idst_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    idst_functor(const T &x)
        : m_result(new Scalar[x.size()])
    {
        typename type_eval<T>::type m_x(x.eval());
        idst_impl<k>((int)m_x.size(), m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            Index n = x.size();
            Scalar scale = (Scalar)idst_scale<k>(m_x.size());
            for (int i = 0; i < n; ++i) {
                p[i] /= scale;
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

template <bool normalize = false, kind k = DST_I, typename T = void>
inline CwiseNullaryOp<idst_functor<normalize, k, T>,
                      typename idst_functor<normalize, k, T>::ResultType>
idst(const DenseBase<T> &x)
{
    static_assert(IS_DST(k), "not dst kind");

    using ResultType = typename idst_functor<normalize, k, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   idst_functor<normalize, k, T>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_IDST__ */
