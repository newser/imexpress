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

#ifndef __IEXP_SPECIAL_AIRY_ZERO__
#define __IEXP_SPECIAL_AIRY_ZERO__

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
// location of n-th 0 of airy Ai
// ========================================

inline double airy_n0_Ai_impl(const unsigned int x)
{
    return gsl_sf_airy_zero_Ai(x);
}

template <typename T>
class airy_n0_Ai_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Ai_functor(const T &x)
        : m_x(x)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Ai_impl((unsigned int)m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Ai_functor<T>,
                      typename airy_n0_Ai_functor<T>::ArrayType>
airy_n0_Ai(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Ai_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Ai_functor<T>(x.derived()));
}

inline double airy_n0_Ai_e_impl(const unsigned int x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_zero_Ai_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_zero_Ai"));
}

template <typename T, typename U>
class airy_n0_Ai_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Ai_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Ai_e_impl(m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Ai_e_functor<T, U>,
                      typename airy_n0_Ai_e_functor<T, U>::ArrayType>
airy_n0_Ai(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Ai_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Ai_e_functor<T, U>(x.derived(),
                                                             e.derived()));
}

// ========================================
// location of n-th 0 of airy Ai derivative
// ========================================

inline double airy_n0_Ai_deriv_impl(const unsigned int x)
{
    return gsl_sf_airy_zero_Ai_deriv(x);
}

template <typename T>
class airy_n0_Ai_deriv_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Ai_deriv_functor(const T &x)
        : m_x(x)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Ai_deriv_impl((unsigned int)m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Ai_deriv_functor<T>,
                      typename airy_n0_Ai_deriv_functor<T>::ArrayType>
airy_n0_Ai_deriv(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Ai_deriv_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Ai_deriv_functor<T>(x.derived()));
}

inline double airy_n0_Ai_deriv_e_impl(const unsigned int x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_zero_Ai_deriv_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_zero_Ai_deriv"));
}

template <typename T, typename U>
class airy_n0_Ai_deriv_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Ai_deriv_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Ai_deriv_e_impl((unsigned int)m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Ai_deriv_e_functor<T, U>,
                      typename airy_n0_Ai_deriv_e_functor<T, U>::ArrayType>
airy_n0_Ai_deriv(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Ai_deriv_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Ai_deriv_e_functor<T,
                                                             U>(x.derived(),
                                                                e.derived()));
}

// ========================================
// location of n-th 0 of airy Bi
// ========================================

inline double airy_n0_Bi_impl(const unsigned int x)
{
    return gsl_sf_airy_zero_Bi(x);
}

template <typename T>
class airy_n0_Bi_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Bi_functor(const T &x)
        : m_x(x)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Bi_impl((unsigned int)m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Bi_functor<T>,
                      typename airy_n0_Bi_functor<T>::ArrayType>
airy_n0_Bi(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Bi_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Bi_functor<T>(x.derived()));
}

inline double airy_n0_Bi_e_impl(const unsigned int x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_zero_Bi_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_zero_Bi"));
}

template <typename T, typename U>
class airy_n0_Bi_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Bi_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Bi_e_impl(m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Bi_e_functor<T, U>,
                      typename airy_n0_Bi_e_functor<T, U>::ArrayType>
airy_n0_Bi(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Bi_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Bi_e_functor<T, U>(x.derived(),
                                                             e.derived()));
}

// ========================================
// location of n-th 0 of airy Bi derivative
// ========================================

inline double airy_n0_Bi_deriv_impl(const unsigned int x)
{
    return gsl_sf_airy_zero_Bi_deriv(x);
}

template <typename T>
class airy_n0_Bi_deriv_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Bi_deriv_functor(const T &x)
        : m_x(x)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Bi_deriv_impl((unsigned int)m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Bi_deriv_functor<T>,
                      typename airy_n0_Bi_deriv_functor<T>::ArrayType>
airy_n0_Bi_deriv(const ArrayBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Bi_deriv_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Bi_deriv_functor<T>(x.derived()));
}

inline double airy_n0_Bi_deriv_e_impl(const unsigned int x, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_zero_Bi_deriv_e(x, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("airy_zero_Bi_deriv"));
}

template <typename T, typename U>
class airy_n0_Bi_deriv_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_n0_Bi_deriv_e_functor(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return airy_n0_Bi_deriv_e_impl((unsigned int)m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Bi_deriv_e_functor<T, U>,
                      typename airy_n0_Bi_deriv_e_functor<T, U>::ArrayType>
airy_n0_Bi_deriv(const ArrayBase<T> &x, ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "x can only be int or uint array");

    typedef typename airy_n0_Bi_deriv_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_n0_Bi_deriv_e_functor<T,
                                                             U>(x.derived(),
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

#endif /* __IEXP_SPECIAL_AIRY_ZERO__ */
