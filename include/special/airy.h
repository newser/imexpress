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

#ifndef __IEXP_SPECIAL_AIRY__
#define __IEXP_SPECIAL_AIRY__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// airy Ai
// ========================================

template <bool scaled, typename T>
inline T airy_Ai_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_impl<false, double>(const double x, const precision p)
{
    return gsl_sf_airy_Ai(x, (int)p);
}

template <>
inline double airy_Ai_impl<true, double>(const double x, const precision p)
{
    return gsl_sf_airy_Ai_scaled(x, (int)p);
}

template <bool scaled, typename T>
class airy_Ai_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_impl<scaled, typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <bool scaled = false, typename T = void>
inline CwiseNullaryOp<airy_Ai_functor<scaled, T>,
                      typename airy_Ai_functor<scaled, T>::ArrayType>
airy_Ai(const ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_functor<scaled, T>(x.derived(), p));
}

template <bool scaled, typename T>
inline T airy_Ai_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_e_impl<false, double>(const double x,
                                            const precision p,
                                            double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Ai"));
}

template <>
inline double airy_Ai_e_impl<true, double>(const double x,
                                           const precision p,
                                           double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Ai_scaled"));
}

template <bool scaled, typename T, typename U>
class airy_Ai_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_e_impl<scaled, typename T::Scalar>(m_x(i, j),
                                                          m_p,
                                                          m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <bool scaled = false, typename T = void, typename U = void>
inline CwiseNullaryOp<airy_Ai_e_functor<scaled, T, U>,
                      typename airy_Ai_e_functor<scaled, T, U>::ArrayType>
airy_Ai(const ArrayBase<T> &x, ArrayBase<U> &e, const precision p = DOUBLE)
{
    typedef typename airy_Ai_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_e_functor<scaled, T, U>(x.derived(),
                                                                  p,
                                                                  e.derived()));
}

// ========================================
// airy Ai derivative
// ========================================

template <bool scaled, typename T>
inline T airy_Ai_deriv_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_impl<false, double>(const double x,
                                                const precision p)
{
    return gsl_sf_airy_Ai_deriv(x, (int)p);
}

template <>
inline double airy_Ai_deriv_impl<true, double>(const double x,
                                               const precision p)
{
    return gsl_sf_airy_Ai_deriv_scaled(x, (int)p);
}

template <bool scaled, typename T>
class airy_Ai_deriv_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_deriv_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_deriv_impl<scaled, typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <bool scaled = false, typename T = void>
inline CwiseNullaryOp<airy_Ai_deriv_functor<scaled, T>,
                      typename airy_Ai_deriv_functor<scaled, T>::ArrayType>
airy_Ai_deriv(const ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_deriv_functor<scaled, T>(x.derived(),
                                                                   p));
}

template <bool scaled, typename T>
inline T airy_Ai_deriv_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_e_impl<false, double>(const double x,
                                                  const precision p,
                                                  double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_deriv_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Ai_deriv"));
}

template <>
inline double airy_Ai_deriv_e_impl<true, double>(const double x,
                                                 const precision p,
                                                 double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_deriv_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Ai_deriv_scaled"));
}

template <bool scaled, typename T, typename U>
class airy_Ai_deriv_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_deriv_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_deriv_e_impl<scaled, typename T::Scalar>(m_x(i, j),
                                                                m_p,
                                                                m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <bool scaled = false, typename T = void, typename U = void>
inline CwiseNullaryOp<airy_Ai_deriv_e_functor<scaled, T, U>,
                      typename airy_Ai_deriv_e_functor<scaled, T, U>::ArrayType>
airy_Ai_deriv(const ArrayBase<T> &x,
              ArrayBase<U> &e,
              const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Ai_deriv_e_functor<scaled, T, U>(x.derived(),
                                                          p,
                                                          e.derived()));
}

// ========================================
// airy Bi
// ========================================

template <bool scaled, typename T>
inline T airy_Bi_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_impl<false, double>(const double x, const precision p)
{
    return gsl_sf_airy_Bi(x, (int)p);
}

template <>
inline double airy_Bi_impl<true, double>(const double x, const precision p)
{
    return gsl_sf_airy_Bi_scaled(x, (int)p);
}

template <bool scaled, typename T>
class airy_Bi_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_impl<scaled, typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <bool scaled = false, typename T = void>
inline CwiseNullaryOp<airy_Bi_functor<scaled, T>,
                      typename airy_Bi_functor<scaled, T>::ArrayType>
airy_Bi(const ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_functor<scaled, T>(x.derived(), p));
}

template <bool scaled, typename T>
inline T airy_Bi_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_e_impl<false, double>(const double x,
                                            const precision p,
                                            double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Bi"));
}

template <>
inline double airy_Bi_e_impl<true, double>(const double x,
                                           const precision p,
                                           double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Bi_scaled"));
}

template <bool scaled, typename T, typename U>
class airy_Bi_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_e_impl<scaled, typename T::Scalar>(m_x(i, j),
                                                          m_p,
                                                          m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <bool scaled = false, typename T = void, typename U = void>
inline CwiseNullaryOp<airy_Bi_e_functor<scaled, T, U>,
                      typename airy_Bi_e_functor<scaled, T, U>::ArrayType>
airy_Bi(const ArrayBase<T> &x, ArrayBase<U> &e, const precision p = DOUBLE)
{
    typedef typename airy_Bi_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_e_functor<scaled, T, U>(x.derived(),
                                                                  p,
                                                                  e.derived()));
}

// ========================================
// airy Bi derivative
// ========================================

template <bool scaled, typename T>
inline T airy_Bi_deriv_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_impl<false, double>(const double x,
                                                const precision p)
{
    return gsl_sf_airy_Bi_deriv(x, (int)p);
}

template <>
inline double airy_Bi_deriv_impl<true, double>(const double x,
                                               const precision p)
{
    return gsl_sf_airy_Bi_deriv_scaled(x, (int)p);
}

template <bool scaled, typename T>
class airy_Bi_deriv_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_deriv_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_deriv_impl<scaled, typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <bool scaled = false, typename T = void>
inline CwiseNullaryOp<airy_Bi_deriv_functor<scaled, T>,
                      typename airy_Bi_deriv_functor<scaled, T>::ArrayType>
airy_Bi_deriv(const ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_functor<scaled, T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_deriv_functor<scaled, T>(x.derived(),
                                                                   p));
}

template <bool scaled, typename T>
inline T airy_Bi_deriv_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_e_impl<false, double>(const double x,
                                                  const precision p,
                                                  double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_deriv_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Bi_deriv"));
}

template <>
inline double airy_Bi_deriv_e_impl<true, double>(const double x,
                                                 const precision p,
                                                 double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_deriv_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_Bi_deriv_scaled"));
}

template <bool scaled, typename T, typename U>
class airy_Bi_deriv_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_deriv_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_deriv_e_impl<scaled, typename T::Scalar>(m_x(i, j),
                                                                m_p,
                                                                m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <bool scaled = false, typename T = void, typename U = void>
inline CwiseNullaryOp<airy_Bi_deriv_e_functor<scaled, T, U>,
                      typename airy_Bi_deriv_e_functor<scaled, T, U>::ArrayType>
airy_Bi_deriv(const ArrayBase<T> &x,
              ArrayBase<U> &e,
              const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_e_functor<scaled, T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Bi_deriv_e_functor<scaled, T, U>(x.derived(),
                                                          p,
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

#endif /* __IEXP_SPECIAL_AIRY__ */
