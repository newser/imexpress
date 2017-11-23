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

#ifndef __IEXP_CYLINDRICAL_BESSEL__
#define __IEXP_CYLINDRICAL_BESSEL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// bessel first kind
// ========================================

template <typename T>
inline T cbessel_j_impl(const int n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_j_impl(const int n, const double x)
{
    return gsl_sf_bessel_Jn(n, x);
}

template <typename T>
class cbessel_j_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_j_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_j_impl<typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<cbessel_j_functor<T>,
                      typename cbessel_j_functor<T>::ArrayType>
cbessel_j(const int n, const Eigen::ArrayBase<T> &x)
{
    typedef typename cbessel_j_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_j_functor<T>(n, x.derived()));
}

template <typename T>
inline T cbessel_j_e_impl(const int n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_j_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Jn_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_j"));
}

template <typename T, typename U>
class cbessel_j_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_j_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_j_e_impl<typename T::Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<cbessel_j_e_functor<T, U>,
                      typename cbessel_j_e_functor<T, U>::ArrayType>
cbessel_j(const int n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename cbessel_j_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_j_e_functor<T, U>(n,
                                                            x.derived(),
                                                            e.derived()));
}

// ========================================
// bessel second kind
// ========================================

template <typename T>
inline T cbessel_y_impl(const int n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_y_impl(const int n, const double x)
{
    return gsl_sf_bessel_Yn(n, x);
}

template <typename T>
class cbessel_y_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_y_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_y_impl<typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<cbessel_y_functor<T>,
                      typename cbessel_y_functor<T>::ArrayType>
cbessel_y(const int n, const Eigen::ArrayBase<T> &x)
{
    typedef typename cbessel_y_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_y_functor<T>(n, x.derived()));
}

template <typename T>
inline T cbessel_y_e_impl(const int n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_y_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Yn_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_y"));
}

template <typename T, typename U>
class cbessel_y_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_y_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_y_e_impl<typename T::Scalar>(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<cbessel_y_e_functor<T, U>,
                      typename cbessel_y_e_functor<T, U>::ArrayType>
cbessel_y(const int n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename cbessel_y_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_y_e_functor<T, U>(n,
                                                            x.derived(),
                                                            e.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_CYLINDRICAL_BESSEL__ */
