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
inline void dd_impl(T *dd, const T *xa, const T *ya, const int len)
{
    throw std::invalid_argument("todo");
}

template <>
inline void dd_impl(double *dd,
                    const double *xa,
                    const double *ya,
                    const int len)
{
    gsl_poly_dd_init(dd, xa, ya, len);
}

template <typename T>
class dd_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        ArrayType;

    dd_functor(const T &xa, const T &ya)
        : m_result(xa.size(), 1)
    {
        typename Eigen::internal::eval<T>::type m_xa(xa.eval()),
            m_ya(ya.eval());
        dd_impl<typename T::Scalar>(m_result.data(),
                                    m_xa.data(),
                                    m_ya.data(),
                                    m_xa.size());
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<dd_functor<T>, typename dd_functor<T>::ArrayType> dd(
    const Eigen::ArrayBase<T> &xa, const Eigen::ArrayBase<T> &ya)
{
    eigen_assert((xa.derived().cols() == 1) || (xa.derived().rows() == 1));
    eigen_assert((ya.derived().cols() == 1) || (ya.derived().rows() == 1));
    eigen_assert(xa.derived().size() == ya.derived().size());

    typedef typename dd_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(xa.derived().size(),
                                  dd_functor<T>(xa.derived(), ya.derived()));
}

// ========================================
// divided difference evaluation
// ========================================

template <typename T>
inline double dd_eval_impl(const T *dd, const T *xa, const int len, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double dd_eval_impl(const double *dd,
                           const double *xa,
                           const int len,
                           const double x)
{
    return gsl_poly_dd_eval(dd, xa, len, x);
}

template <typename T>
inline typename T::Scalar dd_eval(const Eigen::ArrayBase<T> &dd,
                                  const Eigen::ArrayBase<T> &xa,
                                  const typename T::Scalar x)
{
    eigen_assert((dd.derived().cols() == 1) || (dd.derived().rows() == 1));
    eigen_assert((xa.derived().cols() == 1) || (xa.derived().rows() == 1));
    eigen_assert(dd.derived().size() == xa.derived().size());

    typename Eigen::internal::eval<T>::type m_dd(dd.eval()), m_xa(xa.eval());
    return dd_eval_impl<typename T::Scalar>(m_dd.data(),
                                            m_xa.data(),
                                            m_dd.size(),
                                            x);
}

// ========================================
// divided difference to taylor expansion
// ========================================

template <typename T>
inline void dd_taylor_impl(
    T *c, const T xp, const T *dd, const T *xa, const int len, T *w)
{
    throw std::invalid_argument("todo");
}

template <>
inline void dd_taylor_impl(double *c,
                           const double xp,
                           const double *dd,
                           const double *xa,
                           const int len,
                           double *w)
{
    gsl_poly_dd_taylor(c, xp, dd, xa, len, w);
}

template <typename T>
class dd_taylor_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        ArrayType;

    dd_taylor_functor(typename T::Scalar xp, const T &dd, const T &xa)
        : m_result(dd.size(), 1)
    {
        typename Eigen::internal::eval<T>::type m_dd(dd.eval()),
            m_xa(xa.eval());
        ArrayType m_w(dd.size());
        dd_taylor_impl<typename T::Scalar>(m_result.data(),
                                           xp,
                                           m_dd.data(),
                                           m_xa.data(),
                                           m_xa.size(),
                                           m_w.data());
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<dd_taylor_functor<T>,
                      typename dd_taylor_functor<T>::ArrayType>
dd_taylor(typename T::Scalar xp,
          const Eigen::ArrayBase<T> &dd,
          const Eigen::ArrayBase<T> &xa)
{
    eigen_assert((dd.derived().cols() == 1) || (dd.derived().rows() == 1));
    eigen_assert((xa.derived().cols() == 1) || (xa.derived().rows() == 1));
    eigen_assert(dd.derived().size() == xa.derived().size());

    typedef typename dd_taylor_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(xa.derived().size(),
                                  dd_taylor_functor<T>(xp,
                                                       dd.derived(),
                                                       xa.derived()));
}

// ========================================
// hermite representation
// ========================================

template <typename T>
inline void dd_hermit_impl(
    T *dd, T *za, const T *xa, const T *ya, const T *dya, const int len)
{
    throw std::invalid_argument("todo");
}

template <>
inline void dd_hermit_impl(double *dd,
                           double *za,
                           const double *xa,
                           const double *ya,
                           const double *dya,
                           const int len)
{
    gsl_poly_dd_hermite_init(dd, za, xa, ya, dya, len);
}

template <typename T>
class dd_hermit_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::MaxSizeAtCompileTime,
                  1>
        ArrayType;

    dd_hermit_functor(const T &xa, const T &ya, const T &dya)
        : m_result(xa.size() << 1, 1)
    {
        ArrayType m_za(xa.size() << 1);
        typename Eigen::internal::eval<T>::type m_xa(xa.eval()),
            m_ya(ya.eval()), m_dya(dya.eval());
        dd_hermit_impl<typename T::Scalar>(m_result.data(),
                                           m_za.data(),
                                           m_xa.data(),
                                           m_ya.data(),
                                           m_dya.data(),
                                           xa.size());
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result(i);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<dd_hermit_functor<T>,
                      typename dd_hermit_functor<T>::ArrayType>
dd_hermit(const Eigen::ArrayBase<T> &xa,
          const Eigen::ArrayBase<T> &ya,
          const Eigen::ArrayBase<T> &dya)
{
    eigen_assert((xa.derived().cols() == 1) || (xa.derived().rows() == 1));
    eigen_assert((ya.derived().cols() == 1) || (ya.derived().rows() == 1));
    eigen_assert(xa.derived().size() == ya.derived().size());
    eigen_assert((dya.derived().cols() == 1) || (dya.derived().rows() == 1));
    eigen_assert(xa.derived().size() == dya.derived().size());

    typedef typename dd_hermit_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(xa.derived().size() << 1,
                                  dd_hermit_functor<T>(xa.derived(),
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
