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

#ifndef __IEXP_MODIFIED_SPHERICAL_BESSEL__
#define __IEXP_MODIFIED_SPHERICAL_BESSEL__

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
// modified bessel first kind
// ========================================

template <bool scaled, typename T>
inline T msbessel_i_impl(const int n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double msbessel_i_impl<true, double>(const int n, const double x)
{
    return gsl_sf_bessel_il_scaled(n, x);
}

template <bool scaled, typename T>
class msbessel_i_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    msbessel_i_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return msbessel_i_impl<scaled, typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<msbessel_i_functor<scaled, T>,
                      typename msbessel_i_functor<scaled, T>::ArrayType>
msbessel_i(const int n, const Eigen::ArrayBase<T> &x)
{
    typedef typename msbessel_i_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  msbessel_i_functor<scaled, T>(n,
                                                                x.derived()));
}

template <bool scaled, typename T>
inline T msbessel_i_e_impl(const int n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double msbessel_i_e_impl<true, double>(const int n,
                                              const double x,
                                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_il_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("msbessel_i, scaled"));
}

template <bool scaled, typename T, typename U>
class msbessel_i_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    msbessel_i_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return msbessel_i_e_impl<scaled, typename T::Scalar>(m_n,
                                                             m_x(i, j),
                                                             m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<msbessel_i_e_functor<scaled, T, U>,
                      typename msbessel_i_e_functor<scaled, T, U>::ArrayType>
msbessel_i(const int n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename msbessel_i_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    msbessel_i_e_functor<scaled, T, U>(n,
                                                       x.derived(),
                                                       e.derived()));
}

// ========================================
// modified bessel second kind
// ========================================

template <bool scaled, typename T>
inline T msbessel_k_impl(const int n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double msbessel_k_impl<true, double>(const int n, const double x)
{
    return gsl_sf_bessel_kl_scaled(n, x);
}

template <bool scaled, typename T>
class msbessel_k_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    msbessel_k_functor(const int n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return msbessel_k_impl<scaled, typename T::Scalar>(m_n, m_x(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
};

template <bool scaled = true, typename T = void>
inline CwiseNullaryOp<msbessel_k_functor<scaled, T>,
                      typename msbessel_k_functor<scaled, T>::ArrayType>
msbessel_k(const int n, const Eigen::ArrayBase<T> &x)
{
    typedef typename msbessel_k_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  msbessel_k_functor<scaled, T>(n,
                                                                x.derived()));
}

template <bool scaled, typename T>
inline T msbessel_k_e_impl(const int n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double msbessel_k_e_impl<true, double>(const int n,
                                              const double x,
                                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_kl_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("msbessel_k, scaled"));
}

template <bool scaled, typename T, typename U>
class msbessel_k_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    msbessel_k_e_functor(const int n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return msbessel_k_e_impl<scaled, typename T::Scalar>(m_n,
                                                             m_x(i, j),
                                                             m_e(i, j));
    }

  private:
    const int m_n;
    const T &m_x;
    U &m_e;
};

template <bool scaled = true, typename T = void, typename U = void>
inline CwiseNullaryOp<msbessel_k_e_functor<scaled, T, U>,
                      typename msbessel_k_e_functor<scaled, T, U>::ArrayType>
msbessel_k(const int n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename msbessel_k_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    msbessel_k_e_functor<scaled, T, U>(n,
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

#endif /* __IEXP_MODIFIED_SPHERICAL_BESSEL__ */