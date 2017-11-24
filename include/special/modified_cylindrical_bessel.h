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

#ifndef __IEXP_MODIFIED_CYLINDRICAL_BESSEL__
#define __IEXP_MODIFIED_CYLINDRICAL_BESSEL__

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

template <bool scaled, typename T, typename V>
inline T cbessel_i_impl(const V n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_i_impl<false>(const int n, const double x)
{
    // gsl_sf_bessel_In does not call gsl_sf_bessel_i0...
    if (n == 0) {
        return gsl_sf_bessel_I0(x);
    } else if (n == 1) {
        return gsl_sf_bessel_I1(x);
    } else {
        return gsl_sf_bessel_In(n, x);
    }
}

template <>
inline double cbessel_i_impl<true>(const int n, const double x)
{
    return gsl_sf_bessel_In_scaled(n, x);
}

template <>
inline double cbessel_i_impl<false>(const double n, const double x)
{
    return gsl_sf_bessel_Inu(n, x);
}

template <>
inline double cbessel_i_impl<true>(const double n, const double x)
{
    return gsl_sf_bessel_Inu_scaled(n, x);
}

template <bool scaled, typename T, typename V>
class cbessel_i_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_i_functor(const V n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_i_impl<scaled, typename T::Scalar, V>(m_n, m_x(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
};

template <bool scaled = false, typename T = void, typename V = void>
inline CwiseNullaryOp<cbessel_i_functor<scaled, T, V>,
                      typename cbessel_i_functor<scaled, T, V>::ArrayType>
cbessel_i(const V n, const Eigen::ArrayBase<T> &x)
{
    typedef typename cbessel_i_functor<scaled, T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_i_functor<scaled, T, V>(n,
                                                                  x.derived()));
}

template <bool scaled, typename T, typename V>
inline T cbessel_i_e_impl(const V n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_i_e_impl<false>(const int n, const double x, double &e)
{
    gsl_sf_result r;
    int s;
    // gsl_sf_bessel_In does not call gsl_sf_bessel_i0...
    if (n == 0) {
        s = gsl_sf_bessel_I0_e(x, &r);
    } else if (n == 1) {
        s = gsl_sf_bessel_I1_e(x, &r);
    } else {
        s = gsl_sf_bessel_In_e(n, x, &r);
    }
    if (s == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_i"));
}

template <>
inline double cbessel_i_e_impl<true>(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_In_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_i, scaled"));
}

template <>
inline double cbessel_i_e_impl<false>(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Inu_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_i, fractional"));
}

template <>
inline double cbessel_i_e_impl<true>(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Inu_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_i, scaled, fractional"));
}

template <bool scaled, typename T, typename U, typename V>
class cbessel_i_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_i_e_functor(const V n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_i_e_impl<scaled, typename T::Scalar, V>(m_n,
                                                               m_x(i, j),
                                                               m_e(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
    U &m_e;
};

template <bool scaled = false,
          typename T = void,
          typename U = void,
          typename V = void>
inline CwiseNullaryOp<cbessel_i_e_functor<scaled, T, U, V>,
                      typename cbessel_i_e_functor<scaled, T, U, V>::ArrayType>
cbessel_i(const V n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename cbessel_i_e_functor<scaled, T, U, V>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_i_e_functor<scaled, T, U, V>(n,
                                                         x.derived(),
                                                         e.derived()));
}

// ========================================
// modified bessel second kind
// ========================================

template <bool scaled, typename T, typename V>
inline T cbessel_k_impl(const V n, const T x)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_k_impl<false>(const int n, const double x)
{
    return gsl_sf_bessel_Kn(n, x);
}

template <>
inline double cbessel_k_impl<true>(const int n, const double x)
{
    return gsl_sf_bessel_Kn_scaled(n, x);
}

template <>
inline double cbessel_k_impl<false>(const double n, const double x)
{
    return gsl_sf_bessel_Knu(n, x);
}

template <>
inline double cbessel_k_impl<true>(const double n, const double x)
{
    return gsl_sf_bessel_Knu_scaled(n, x);
}

template <bool scaled, typename T, typename V>
class cbessel_k_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_k_functor(const V n, const T &x)
        : m_n(n)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_k_impl<scaled, typename T::Scalar, V>(m_n, m_x(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
};

template <bool scaled = false, typename T = void, typename V = void>
inline CwiseNullaryOp<cbessel_k_functor<scaled, T, V>,
                      typename cbessel_k_functor<scaled, T, V>::ArrayType>
cbessel_k(const V n, const Eigen::ArrayBase<T> &x)
{
    typedef typename cbessel_k_functor<scaled, T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  cbessel_k_functor<scaled, T, V>(n,
                                                                  x.derived()));
}

template <bool scaled, typename T, typename V>
inline T cbessel_k_e_impl(const V n, const T x, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double cbessel_k_e_impl<false>(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Kn_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_k"));
}

template <>
inline double cbessel_k_e_impl<true>(const int n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Kn_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_k, scaled"));
}

template <>
inline double cbessel_k_e_impl<false>(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Knu_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_k, fractional"));
}

template <>
inline double cbessel_k_e_impl<true>(const double n, const double x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_bessel_Knu_scaled_e(n, x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("cbessel_k, scaled, fractional"));
}

template <bool scaled, typename T, typename U, typename V>
class cbessel_k_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    cbessel_k_e_functor(const V n, const T &x, U &e)
        : m_n(n)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return cbessel_k_e_impl<scaled, typename T::Scalar, V>(m_n,
                                                               m_x(i, j),
                                                               m_e(i, j));
    }

  private:
    const V m_n;
    const T &m_x;
    U &m_e;
};

template <bool scaled = false,
          typename T = void,
          typename U = void,
          typename V = void>
inline CwiseNullaryOp<cbessel_k_e_functor<scaled, T, U, V>,
                      typename cbessel_k_e_functor<scaled, T, U, V>::ArrayType>
cbessel_k(const V n, const Eigen::ArrayBase<T> &x, Eigen::ArrayBase<U> &e)
{
    typedef typename cbessel_k_e_functor<scaled, T, U, V>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    cbessel_k_e_functor<scaled, T, U, V>(n,
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

#endif /* __IEXP_MODIFIED_CYLINDRICAL_BESSEL__ */
