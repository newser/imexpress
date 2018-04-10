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
inline void ifft2_impl(int n0, int n1, const T *i, U *o)
{
    fftw3::get_plan(n0, n1, i, o, false).inv(n0, n1, i, o);
}

template <bool normalize,
          typename T,
          bool is_row_major = bool(TP4(T) == RowMajor)>
class ifft2_functor;

// c2c, row major
template <bool normalize, typename T>
class ifft2_functor<normalize, T, true>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T, Scalar>::type;

    ifft2_functor(const T &x)
        : m_n0(x.rows())
        , m_n1(x.cols())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        ifft2_impl((int)m_n0, (int)m_n1, m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            Index n = x.size();
            for (int i = 0; i < n; ++i) {
                p[i] /= n;
            }
        }
    }

    Scalar operator()(Index i, Index j) const
    {
#define IFFT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return IFFT2_RESULT(i, j);
#undef IFFT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

// c2c, col major
template <bool normalize, typename T>
class ifft2_functor<normalize, T, false>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T, Scalar>::type;

    ifft2_functor(const T &x)
        : m_n0(x.cols())
        , m_n1(x.rows())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        ifft2_impl((int)m_n0, (int)m_n1, m_x.data(), m_result.get());

        if (normalize) {
            Scalar *p = m_result.get();
            Index n = x.size();
            for (int i = 0; i < n; ++i) {
                p[i] /= n;
            }
        }
    }

    Scalar operator()(Index i, Index j) const
    {
#define IFFT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return IFFT2_RESULT(j, i);
#undef IFFT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

template <bool normalize = false, typename T = void>
inline CwiseNullaryOp<ifft2_functor<normalize, T>,
                      typename ifft2_functor<normalize, T>::ResultType>
ifft2(const DenseBase<T> &x)
{
    using ResultType = typename ifft2_functor<normalize, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
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
