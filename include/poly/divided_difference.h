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

#ifndef __IEXP_POLY_DIVIDED_DIFFERENCE__
#define __IEXP_POLY_DIVIDED_DIFFERENCE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_poly.h>

IEXP_NS_BEGIN

namespace poly {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// divided difference
// ========================================

template <typename T>
inline void dd_impl(T *dd, const T *xa, const T *ya, int len)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void dd_impl(double *dd, const double *xa, const double *ya, int len)
{
    gsl_poly_dd_init(dd, xa, ya, len);
}

template <typename T>
class dd_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    dd_functor(const T &xa, const T &ya)
        : m_result(new Scalar[xa.size()])
    {
        typename type_eval<T>::type m_xa(xa.eval()), m_ya(ya.eval());
        dd_impl(m_result.get(), m_xa.data(), m_ya.data(), m_xa.size());
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <typename T>
inline CwiseNullaryOp<dd_functor<T>, typename dd_functor<T>::ResultType> dd(
    const DenseBase<T> &xa, const DenseBase<T> &ya)
{
    eigen_assert(MATRIX_SAME_SIZE(xa, ya));

    using ResultType = typename dd_functor<T>::ResultType;
    return ResultType::NullaryExpr(xa.rows(),
                                   xa.cols(),
                                   dd_functor<T>(xa.derived(), ya.derived()));
}

// ========================================
// divided difference evaluation
// ========================================

template <typename T>
inline double dd_eval_impl(const T *dd, const T *xa, int len, T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double dd_eval_impl(const double *dd,
                           const double *xa,
                           int len,
                           double x)
{
    return gsl_poly_dd_eval(dd, xa, len, x);
}

template <typename T>
inline typename T::Scalar dd_eval(const DenseBase<T> &dd,
                                  const DenseBase<T> &xa,
                                  typename T::Scalar x)
{
    eigen_assert(dd.size() == xa.size());

    typename type_eval<T>::type m_dd(dd.eval()), m_xa(xa.eval());
    return dd_eval_impl(m_dd.data(), m_xa.data(), m_dd.size(), x);
}

template <typename T, typename U>
class dd_eval_functor
{
  public:
    using Scalar = typename U::Scalar;
    using ResultType = typename dense_derive<U>::type;

    dd_eval_functor(const DenseBase<T> &dd,
                    const DenseBase<T> &xa,
                    const DenseBase<U> &x)
        : m_result(new Scalar[x.size()])
    {
        typename type_eval<T>::type m_dd(dd.eval()), m_xa(xa.eval());
        ASSIGN_IJ(TP4(U) == RowMajor,
                  m_result.get(),
                  x,
                  dd_eval_impl(m_dd.data(), m_xa.data(), m_dd.size(), x(i, j)));
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <typename T, typename U>
inline CwiseNullaryOp<dd_eval_functor<T, U>,
                      typename dd_eval_functor<T, U>::ResultType>
dd_eval(const DenseBase<T> &dd, const DenseBase<T> &xa, const DenseBase<U> &x)
{
    eigen_assert(MATRIX_SAME_SIZE(dd, xa));

    using ResultType = typename dd_eval_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dd_eval_functor<T, U>(dd.derived(),
                                                         xa.derived(),
                                                         x.derived()));
}

// ========================================
// divided difference to taylor expansion
// ========================================

template <typename T>
inline void dd_taylor_impl(T *c, T xp, const T *dd, const T *xa, int len, T *w)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void dd_taylor_impl(double *c,
                           double xp,
                           const double *dd,
                           const double *xa,
                           int len,
                           double *w)
{
    gsl_poly_dd_taylor(c, xp, dd, xa, len, w);
}

template <typename T>
class dd_taylor_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T>::type;

    dd_taylor_functor(Scalar xp, const T &dd, const T &xa)
        : m_result(new Scalar[dd.size()])
    {
        typename type_eval<T>::type m_dd(dd.eval()), m_xa(xa.eval());
        std::unique_ptr<Scalar> m_w(new Scalar[dd.size()]);
        dd_taylor_impl(m_result.get(),
                       xp,
                       m_dd.data(),
                       m_xa.data(),
                       m_xa.size(),
                       m_w.get());
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <typename T>
inline CwiseNullaryOp<dd_taylor_functor<T>,
                      typename dd_taylor_functor<T>::ResultType>
dd_taylor(typename T::Scalar xp, const DenseBase<T> &dd, const DenseBase<T> &xa)
{
    eigen_assert(MATRIX_SAME_SIZE(dd, xa));

    using ResultType = typename dd_taylor_functor<T>::ResultType;
    return ResultType::NullaryExpr(dd.rows(),
                                   dd.cols(),
                                   dd_taylor_functor<T>(xp,
                                                        dd.derived(),
                                                        xa.derived()));
}

// ========================================
// hermite representation
// ========================================

template <typename T>
inline void dd_hermit_impl(
    T *dd, T *za, const T *xa, const T *ya, const T *dya, int len)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void dd_hermit_impl(double *dd,
                           double *za,
                           const double *xa,
                           const double *ya,
                           const double *dya,
                           int len)
{
    gsl_poly_dd_hermite_init(dd, za, xa, ya, dya, len);
}

template <bool row_form, typename T>
class dd_hermit_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    dd_hermit_functor(const T &xa, const T &ya, const T &dya)
        : m_result(new Scalar[xa.size() << 1])
    {
        std::unique_ptr<Scalar> m_za(new Scalar[xa.size() << 1]);
        typename type_eval<T>::type m_xa(xa.eval()), m_ya(ya.eval()),
            m_dya(dya.eval());
        dd_hermit_impl(m_result.get(),
                       m_za.get(),
                       m_xa.data(),
                       m_ya.data(),
                       m_dya.data(),
                       xa.size());
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<dd_hermit_functor<row_form, T>,
                      typename dd_hermit_functor<row_form, T>::ResultType>
dd_hermit(const DenseBase<T> &xa,
          const DenseBase<T> &ya,
          const DenseBase<T> &dya)
{
    eigen_assert(MATRIX_SAME_SIZE(xa, ya));
    eigen_assert(MATRIX_SAME_SIZE(xa, dya));

    using ResultType = typename dd_hermit_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(row_form ? 1 : (xa.size() << 1),
                                   row_form ? (xa.size() << 1) : 1,
                                   dd_hermit_functor<row_form,
                                                     T>(xa.derived(),
                                                        ya.derived(),
                                                        dya.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_POLY_DIVIDED_DIFFERENCE__ */
