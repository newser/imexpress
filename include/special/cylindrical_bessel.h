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

template <typename T, typename V>
inline T cbessel_j_impl(const V n, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cbessel_j_impl(const int n, const double x)
{
    return gsl_sf_bessel_Jn(n, x);
}

template <>
inline double cbessel_j_impl(const double n, const double x)
{
    return gsl_sf_bessel_Jnu(n, x);
}

template <typename T, typename V>
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

    cbessel_j_functor(const V n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_j_impl(m_n, m_x(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
};

template <typename T, typename V>
inline CwiseNullaryOp<cbessel_j_functor<T, V>,
                      typename cbessel_j_functor<T, V>::ArrayType>
cbessel_j(const V n, const ArrayBase<T> &x)
{
    typedef typename cbessel_j_functor<T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_j_functor<T, V>(n, x.derived()));
}

template <typename T, typename V>
inline T cbessel_j_e_impl(const V n, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cbessel_j_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Jn_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_j"));
}

template <>
inline double cbessel_j_e_impl(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Jnu_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_j"));
}

template <typename T, typename U, typename V>
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

    cbessel_j_e_functor(const V n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_j_e_impl(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U, typename V>
inline CwiseNullaryOp<cbessel_j_e_functor<T, U, V>,
                      typename cbessel_j_e_functor<T, U, V>::ArrayType>
cbessel_j(const V n, const ArrayBase<T> &x, ArrayBase<U> &e)
{
    typedef typename cbessel_j_e_functor<T, U, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_j_e_functor<T, U, V>(n,
                                                               x.derived(),
                                                               e.derived()));
}

// ========================================
// bessel second kind
// ========================================

template <typename T, typename V>
inline T cbessel_y_impl(const V n, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cbessel_y_impl(const int n, const double x)
{
    return gsl_sf_bessel_Yn(n, x);
}

template <>
inline double cbessel_y_impl(const double n, const double x)
{
    return gsl_sf_bessel_Ynu(n, x);
}

template <typename T, typename V>
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

    cbessel_y_functor(const V n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_y_impl(m_n, m_x(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
};

template <typename T, typename V>
inline CwiseNullaryOp<cbessel_y_functor<T, V>,
                      typename cbessel_y_functor<T, V>::ArrayType>
cbessel_y(const V n, const ArrayBase<T> &x)
{
    typedef typename cbessel_y_functor<T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_y_functor<T, V>(n, x.derived()));
}

template <typename T, typename V>
inline T cbessel_y_e_impl(const V n, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cbessel_y_e_impl(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Yn_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_y"));
}

template <>
inline double cbessel_y_e_impl(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Ynu_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_y"));
}

template <typename T, typename U, typename V>
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

    cbessel_y_e_functor(const V n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_y_e_impl(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U, typename V>
inline CwiseNullaryOp<cbessel_y_e_functor<T, U, V>,
                      typename cbessel_y_e_functor<T, U, V>::ArrayType>
cbessel_y(const V n, const ArrayBase<T> &x, ArrayBase<U> &e)
{
    typedef typename cbessel_y_e_functor<T, U, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_y_e_functor<T, U, V>(n,
                                                               x.derived(),
                                                               e.derived()));
}

// ========================================
// location of x-th 0 of first kind bessel
// ========================================

template <typename V>
inline double cbessel_n0_j_impl(const V n, const unsigned int x)
{
    UNSUPPORTED_TYPE(V);
}

template <>
inline double cbessel_n0_j_impl(const int n, const unsigned int x)
{
    if (n == 0) {
        return gsl_sf_bessel_zero_J0(x);
    } else if (n == 1) {
        return gsl_sf_bessel_zero_J1(x);
    } else {
        RETURN_NAN_OR_THROW(std::invalid_argument("use 'n.0' instead"));
    }
}

template <>
inline double cbessel_n0_j_impl(const double n, const unsigned int x)
{
    return gsl_sf_bessel_zero_Jnu(n, x);
}

template <typename T, typename V>
class cbessel_n0_j_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_n0_j_functor(const V n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return cbessel_n0_j_impl(m_n, (unsigned int)m_x(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
};

template <typename T, typename V>
inline CwiseNullaryOp<cbessel_n0_j_functor<T, V>,
                      typename cbessel_n0_j_functor<T, V>::ArrayType>
cbessel_n0_j(const V n, const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename cbessel_n0_j_functor<T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_n0_j_functor<T, V>(n, x.derived()));
}

template <typename V, typename T>
inline double cbessel_n0_j_e_impl(const V n, const unsigned int x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cbessel_n0_j_e_impl(const int n, const unsigned int x, double &e)
{
    gsl_sf_result r;
    int s;
    if (n == 0) {
        s = gsl_sf_bessel_zero_J0_e(x, &r);
    } else if (n == 1) {
        s = gsl_sf_bessel_zero_J1_e(x, &r);
    } else {
        RETURN_NAN_OR_THROW(std::invalid_argument("use 'n.0' instead"));
    }
    if (s == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_n0_j"));
}

template <>
inline double cbessel_n0_j_e_impl(const double n,
                                  const unsigned int x,
                                  double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_zero_Jnu_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("cbessel_n0_j"));
}

template <typename T, typename U, typename V>
class cbessel_n0_j_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_n0_j_e_functor(const V n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return cbessel_n0_j_e_impl(m_n, m_x(i, j), m_e(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U, typename V>
inline CwiseNullaryOp<cbessel_n0_j_e_functor<T, U, V>,
                      typename cbessel_n0_j_e_functor<T, U, V>::ArrayType>
cbessel_n0_j(const V n, const ArrayBase<T> &x, ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename cbessel_n0_j_e_functor<T, U, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_n0_j_e_functor<T, U, V>(n,
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
