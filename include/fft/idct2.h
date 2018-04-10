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

#ifndef __IEXP_IDCT2__
#define __IEXP_IDCT2__

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
inline void idct2_impl(int n0, int n1, const T *i, T *o)
{
    fftw3::get_plan<k0, k1>(n0, n1, i, o, false)
        .template inv<k0, k1>(n0, n1, i, o);
}

template <kind k0, kind k1>
inline size_t idct2_scale(size_t rows, size_t cols)
{
    size_t r = k0 == DCT_I ? rows - 1 : rows;
    size_t c = k1 == DCT_I ? cols - 1 : cols;
    return (r * c) << 2;
}

template <bool normalize,
          kind k0,
          kind k1,
          typename T,
          bool is_row_major = bool(TP4(T) == RowMajor)>
class idct2_functor;

template <bool normalize, kind k0, kind k1, typename T>
class idct2_functor<normalize, k0, k1, T, true>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    idct2_functor(const T &x)
        : m_n0(x.rows())
        , m_n1(x.cols())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        idct2_impl<k0, k1>((int)m_n0, (int)m_n1, m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            size_t n = x.size();
            Scalar scale = idct2_scale<k0, k1>(m_n0, m_n1);
            for (int i = 0; i < n; ++i) {
                p[i] /= scale;
            }
        }
    }

    Scalar operator()(Index i, Index j) const
    {
#define IDCT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return IDCT2_RESULT(i, j);
#undef IDCT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

template <bool normalize, kind k0, kind k1, typename T>
class idct2_functor<normalize, k0, k1, T, false>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    idct2_functor(const T &x)
        : m_n0(x.cols())
        , m_n1(x.rows())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        idct2_impl<k0, k1>((int)m_n0, (int)m_n1, m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            size_t n = x.size();
            Scalar scale = idct2_scale<k0, k1>(m_n0, m_n1);
            for (int i = 0; i < n; ++i) {
                p[i] /= scale;
            }
        }
    }

    Scalar operator()(Index i, Index j) const
    {
#define IDCT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return IDCT2_RESULT(j, i);
#undef IDCT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

template <bool normalize = false,
          kind k0 = DCT_II,
          kind k1 = DCT_II,
          typename T = void>
inline CwiseNullaryOp<idct2_functor<normalize, k0, k1, T>,
                      typename idct2_functor<normalize, k0, k1, T>::ResultType>
idct2(const DenseBase<T> &x)
{
    static_assert(IS_DCT(k0) && IS_DCT(k1), "not dct kind");

    using ResultType = typename idct2_functor<normalize, k0, k1, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   idct2_functor<normalize, k0, k1, T>(
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

#endif /* __IEXP_IDCT2__ */
