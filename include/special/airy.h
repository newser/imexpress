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

template <typename T>
inline T airy_Ai_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_impl(const double x, const precision p)
{
    return gsl_sf_airy_Ai(x, (int)p);
}

template <typename T>
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
        return airy_Ai_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Ai_functor<T>,
                      typename airy_Ai_functor<T>::ArrayType>
airy_Ai(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Ai_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_e_impl(const double x, const precision p, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
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
        return airy_Ai_e_impl<typename T::Scalar>(m_x(i, j), m_p, m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Ai_e_functor<T, U>,
                      typename airy_Ai_e_functor<T, U>::ArrayType>
airy_Ai(const Eigen::ArrayBase<T> &x,
        Eigen::ArrayBase<U> &e,
        const precision p = DOUBLE)
{
    typedef typename airy_Ai_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_e_functor<T, U>(x.derived(),
                                                          p,
                                                          e.derived()));
}

// ========================================
// airy Ai scaled with exp(+(2/3) x^{3/2})
// ============================￼===========

template <typename T>
inline T airy_Ai_scaled_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_scaled_impl(const double x, const precision p)
{
    return gsl_sf_airy_Ai_scaled(x, (int)p);
}

template <typename T>
class airy_Ai_scaled_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_scaled_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_scaled_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Ai_scaled_functor<T>,
                      typename airy_Ai_scaled_functor<T>::ArrayType>
airy_Ai_scaled(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_scaled_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_scaled_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Ai_scaled_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_scaled_e_impl(const double x,
                                    const precision p,
                                    double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
class airy_Ai_scaled_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_scaled_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_scaled_e_impl<typename T::Scalar>(m_x(i, j),
                                                         m_p,
                                                         m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Ai_scaled_e_functor<T, U>,
                      typename airy_Ai_scaled_e_functor<T, U>::ArrayType>
airy_Ai_scaled(const Eigen::ArrayBase<T> &x,
               Eigen::ArrayBase<U> &e,
               const precision p = DOUBLE)
{
    typedef typename airy_Ai_scaled_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_scaled_e_functor<T, U>(x.derived(),
                                                                 p,
                                                                 e.derived()));
}

// ========================================
// airy Ai derivative
// ========================================

template <typename T>
inline T airy_Ai_deriv_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_impl(const double x, const precision p)
{
    return gsl_sf_airy_Ai_deriv(x, (int)p);
}

template <typename T>
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
        return airy_Ai_deriv_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Ai_deriv_functor<T>,
                      typename airy_Ai_deriv_functor<T>::ArrayType>
airy_Ai_deriv(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_deriv_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Ai_deriv_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_e_impl(const double x, const precision p, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_deriv_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
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
        return airy_Ai_deriv_e_impl<typename T::Scalar>(m_x(i, j),
                                                        m_p,
                                                        m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Ai_deriv_e_functor<T, U>,
                      typename airy_Ai_deriv_e_functor<T, U>::ArrayType>
airy_Ai_deriv(const Eigen::ArrayBase<T> &x,
              Eigen::ArrayBase<U> &e,
              const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_deriv_e_functor<T, U>(x.derived(),
                                                                p,
                                                                e.derived()));
}

// ========================================
// airy Ai derivative scaled with exp(+(2/3) x^{3/2})
// =======================================

template <typename T>
inline T airy_Ai_deriv_scaled_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_scaled_impl(const double x, const precision p)
{
    return gsl_sf_airy_Ai_deriv_scaled(x, (int)p);
}

template <typename T>
class airy_Ai_deriv_scaled_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_deriv_scaled_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_deriv_scaled_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Ai_deriv_scaled_functor<T>,
                      typename airy_Ai_deriv_scaled_functor<T>::ArrayType>
airy_Ai_deriv_scaled(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_scaled_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Ai_deriv_scaled_functor<T>(x.derived(),
                                                                  p));
}

template <typename T>
inline T airy_Ai_deriv_scaled_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Ai_deriv_scaled_e_impl(const double x,
                                          const precision p,
                                          double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Ai_deriv_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
class airy_Ai_deriv_scaled_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Ai_deriv_scaled_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Ai_deriv_scaled_e_impl<typename T::Scalar>(m_x(i, j),
                                                               m_p,
                                                               m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Ai_deriv_scaled_e_functor<T, U>,
                      typename airy_Ai_deriv_scaled_e_functor<T, U>::ArrayType>
airy_Ai_deriv_scaled(const Eigen::ArrayBase<T> &x,
                     Eigen::ArrayBase<U> &e,
                     const precision p = DOUBLE)
{
    typedef typename airy_Ai_deriv_scaled_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Ai_deriv_scaled_e_functor<T, U>(x.derived(),
                                                         p,
                                                         e.derived()));
}

// ========================================
// airy Bi
// ========================================

template <typename T>
inline T airy_Bi_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_impl(const double x, const precision p)
{
    return gsl_sf_airy_Bi(x, (int)p);
}

template <typename T>
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
        return airy_Bi_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Bi_functor<T>,
                      typename airy_Bi_functor<T>::ArrayType>
airy_Bi(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Bi_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_e_impl(const double x, const precision p, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
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
        return airy_Bi_e_impl<typename T::Scalar>(m_x(i, j), m_p, m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Bi_e_functor<T, U>,
                      typename airy_Bi_e_functor<T, U>::ArrayType>
airy_Bi(const Eigen::ArrayBase<T> &x,
        Eigen::ArrayBase<U> &e,
        const precision p = DOUBLE)
{
    typedef typename airy_Bi_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_e_functor<T, U>(x.derived(),
                                                          p,
                                                          e.derived()));
}

// ========================================
// airy Bi scaled with exp(-(2/3) x^{3/2})
// ============================￼===========

template <typename T>
inline T airy_Bi_scaled_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_scaled_impl(const double x, const precision p)
{
    return gsl_sf_airy_Bi_scaled(x, (int)p);
}

template <typename T>
class airy_Bi_scaled_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_scaled_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_scaled_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Bi_scaled_functor<T>,
                      typename airy_Bi_scaled_functor<T>::ArrayType>
airy_Bi_scaled(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_scaled_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_scaled_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Bi_scaled_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_scaled_e_impl(const double x,
                                    const precision p,
                                    double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
class airy_Bi_scaled_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_scaled_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_scaled_e_impl<typename T::Scalar>(m_x(i, j),
                                                         m_p,
                                                         m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Bi_scaled_e_functor<T, U>,
                      typename airy_Bi_scaled_e_functor<T, U>::ArrayType>
airy_Bi_scaled(const Eigen::ArrayBase<T> &x,
               Eigen::ArrayBase<U> &e,
               const precision p = DOUBLE)
{
    typedef typename airy_Bi_scaled_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_scaled_e_functor<T, U>(x.derived(),
                                                                 p,
                                                                 e.derived()));
}

// ========================================
// airy Bi derivative
// ========================================

template <typename T>
inline T airy_Bi_deriv_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_impl(const double x, const precision p)
{
    return gsl_sf_airy_Bi_deriv(x, (int)p);
}

template <typename T>
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
        return airy_Bi_deriv_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Bi_deriv_functor<T>,
                      typename airy_Bi_deriv_functor<T>::ArrayType>
airy_Bi_deriv(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_deriv_functor<T>(x.derived(), p));
}

template <typename T>
inline T airy_Bi_deriv_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_e_impl(const double x, const precision p, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_deriv_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
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
        return airy_Bi_deriv_e_impl<typename T::Scalar>(m_x(i, j),
                                                        m_p,
                                                        m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Bi_deriv_e_functor<T, U>,
                      typename airy_Bi_deriv_e_functor<T, U>::ArrayType>
airy_Bi_deriv(const Eigen::ArrayBase<T> &x,
              Eigen::ArrayBase<U> &e,
              const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_deriv_e_functor<T, U>(x.derived(),
                                                                p,
                                                                e.derived()));
}

// ========================================
// airy Bi derivative scaled with exp(-(2/3) x^{3/2})
// =======================================

template <typename T>
inline T airy_Bi_deriv_scaled_impl(const T x, const precision p)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_scaled_impl(const double x, const precision p)
{
    return gsl_sf_airy_Bi_deriv_scaled(x, (int)p);
}

template <typename T>
class airy_Bi_deriv_scaled_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_deriv_scaled_functor(const T &x, precision p)
        : m_x(x)
        , m_p(p)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_deriv_scaled_impl<typename T::Scalar>(m_x(i, j), m_p);
    }

  private:
    const T &m_x;
    precision m_p;
};

template <typename T>
inline CwiseNullaryOp<airy_Bi_deriv_scaled_functor<T>,
                      typename airy_Bi_deriv_scaled_functor<T>::ArrayType>
airy_Bi_deriv_scaled(const Eigen::ArrayBase<T> &x, const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_scaled_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(x.rows(),
                                  x.cols(),
                                  airy_Bi_deriv_scaled_functor<T>(x.derived(),
                                                                  p));
}

template <typename T>
inline T airy_Bi_deriv_scaled_e_impl(const T x, const precision p, T &e)
{
    throw std::invalid_argument("todo");
}

template <>
inline double airy_Bi_deriv_scaled_e_impl(const double x,
                                          const precision p,
                                          double &e)
{
    gsl_sf_result r;
    if (gsl_sf_airy_Bi_deriv_scaled_e(x, (int)p, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    throw std::invalid_argument("todo");
}

template <typename T, typename U>
class airy_Bi_deriv_scaled_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    airy_Bi_deriv_scaled_e_functor(const T &x, precision p, U &e)
        : m_x(x)
        , m_p(p)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return airy_Bi_deriv_scaled_e_impl<typename T::Scalar>(m_x(i, j),
                                                               m_p,
                                                               m_e(i, j));
    }

  private:
    const T &m_x;
    precision m_p;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_Bi_deriv_scaled_e_functor<T, U>,
                      typename airy_Bi_deriv_scaled_e_functor<T, U>::ArrayType>
airy_Bi_deriv_scaled(const Eigen::ArrayBase<T> &x,
                     Eigen::ArrayBase<U> &e,
                     const precision p = DOUBLE)
{
    typedef typename airy_Bi_deriv_scaled_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Bi_deriv_scaled_e_functor<T, U>(x.derived(),
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
