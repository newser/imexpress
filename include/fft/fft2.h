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

#include <iostream>

IEXP_NS_BEGIN

namespace fft {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T, typename U>
inline void fft2_impl(int n0, int n1, const T *i, U *o)
{
    fftw3::get_plan(n0, n1, i, o, true).fwd(n0, n1, i, o);
}

template <typename T,
          bool is_complex = IS_COMPLEX(typename T::Scalar),
          bool is_row_major = bool(TP4(T) == RowMajor)>
class fft2_functor;

// c2c, row major
template <typename T>
class fft2_functor<T, true, true>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T, Scalar>::type;

    fft2_functor(const T &x)
        : m_n0(x.rows())
        , m_n1(x.cols())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        fft2_impl((int)m_n0, (int)m_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define FFT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return FFT2_RESULT(i, j);
#undef FFT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

// c2c, col major
template <typename T>
class fft2_functor<T, true, false>
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T, Scalar>::type;

    fft2_functor(const T &x)
        : m_n0(x.cols())
        , m_n1(x.rows())
        , m_result(new Scalar[m_n0 * m_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        fft2_impl((int)m_n0, (int)m_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define FFT2_RESULT(i, j) m_result.get()[(i)*m_n1 + (j)]
        return FFT2_RESULT(j, i);
#undef FFT2_RESULT
    }

  private:
    Index m_n0, m_n1;
    std::shared_ptr<Scalar> m_result;
};

// r2c, row major
template <typename T>
class fft2_functor<T, false, true>
{
  public:
    using Scalar = std::complex<typename T::Scalar>;
    using ResultType = typename dense_derive<T, Scalar>::type;

    fft2_functor(const T &x)
        : m_n0(x.rows())
        , m_in_n1(x.cols())
        , m_out_n1((m_in_n1 >> 1) + 1)
        , m_result(new Scalar[m_n0 * m_out_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        fft2_impl((int)m_n0, (int)m_in_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define FFT2_RESULT(i, j) m_result.get()[(i)*m_out_n1 + (j)]
        if (j < m_out_n1) {
            return FFT2_RESULT(i, j);
        } else if (i == 0) {
            return std::conj(FFT2_RESULT(i, m_in_n1 - j));
        } else {
            return std::conj(FFT2_RESULT(m_n0 - i, m_in_n1 - j));
        }
#undef FFT2_RESULT
    }

  private:
    Index m_n0, m_in_n1, m_out_n1;
    std::shared_ptr<Scalar> m_result;
};

// r2c, col major
template <typename T>
class fft2_functor<T, false, false>
{
  public:
    using Scalar = std::complex<typename T::Scalar>;
    using ResultType = typename dense_derive<T, Scalar>::type;

    fft2_functor(const T &x)
        : m_n0(x.cols())
        , m_in_n1(x.rows())
        , m_out_n1((m_in_n1 >> 1) + 1)
        , m_result(new Scalar[m_n0 * m_out_n1])
    {
        typename type_eval<T>::type m_x(x.eval());
        fft2_impl((int)m_n0, (int)m_in_n1, m_x.data(), m_result.get());
    }

    Scalar operator()(Index i, Index j) const
    {
#define FFT2_RESULT(i, j) m_result.get()[(i)*m_out_n1 + (j)]
        if (i < m_out_n1) {
            return FFT2_RESULT(j, i);
        } else if (j == 0) {
            return std::conj(FFT2_RESULT(j, m_in_n1 - i));
        } else {
            return std::conj(FFT2_RESULT(m_n0 - j, m_in_n1 - i));
        }
#undef FFT2_RESULT
    }

  private:
    Index m_n0, m_in_n1, m_out_n1;
    std::shared_ptr<Scalar> m_result;
};

template <typename T>
inline CwiseNullaryOp<fft2_functor<T>, typename fft2_functor<T>::ResultType>
fft2(const DenseBase<T> &x)
{
    using ResultType = typename fft2_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
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
