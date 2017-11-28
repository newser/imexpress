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

#ifndef __IEXP_HYDROGENIC_R__
#define __IEXP_HYDROGENIC_R__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// hydrogenicR_1
// ========================================

template <typename T>
inline T hydroR1_impl(const T z, const T r)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double hydroR1_impl(const double z, const double r)
{
    return gsl_sf_hydrogenicR_1(z, r);
}

template <typename T>
class hydroR1_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    hydroR1_functor(const T &z, const T &r)
        : m_z(z)
        , m_r(r)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return hydroR1_impl(m_z(i, j), m_r(i, j));
    }

  private:
    const T &m_z;
    const T &m_r;
};

template <typename T>
inline CwiseNullaryOp<hydroR1_functor<T>,
                      typename hydroR1_functor<T>::ArrayType>
hydroR(const ArrayBase<T> &z, const ArrayBase<T> &r)
{
    eigen_assert(MATRIX_SAME_SIZE(z, r));

    typedef typename hydroR1_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(z.rows(),
                                  z.cols(),
                                  hydroR1_functor<T>(z.derived(), r.derived()));
}

template <typename T>
inline T hydroR1_e_impl(const T z, const T r, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double hydroR1_e_impl(const double z, const double r, double &e)
{
    gsl_sf_result ret;
    if (gsl_sf_hydrogenicR_1_e(z, r, &ret) == GSL_SUCCESS) {
        e = ret.err;
        return ret.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("hydroR"));
}

template <typename T, typename U>
class hydroR1_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    hydroR1_e_functor(const T &z, const T &r, U &e)
        : m_z(z)
        , m_r(r)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return hydroR1_e_impl(m_z(i, j), m_r(i, j), m_e(i, j));
    }

  private:
    const T &m_z;
    const T &m_r;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<hydroR1_e_functor<T, U>,
                      typename hydroR1_e_functor<T, U>::ArrayType>
hydroR(const ArrayBase<T> &z, const ArrayBase<T> &r, ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(z, r));

    typedef typename hydroR1_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(z.rows(),
                                  z.cols(),
                                  hydroR1_e_functor<T, U>(z.derived(),
                                                          r.derived(),
                                                          e.derived()));
}

// ========================================
// hydrogenicR
// ========================================

template <typename T, typename V>
inline T hydroR_impl(const V n, const V l, const T z, const T r)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double hydroR_impl(const int n,
                          const int l,
                          const double z,
                          const double r)
{
    return gsl_sf_hydrogenicR(n, l, z, r);
}

template <typename T, typename V>
class hydroR_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    hydroR_functor(const V &n, const V &l, const T &z, const T &r)
        : m_n(n)
        , m_l(l)
        , m_z(z)
        , m_r(r)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return hydroR_impl((int)m_n(i, j),
                           (int)m_l(i, j),
                           m_z(i, j),
                           m_r(i, j));
    }

  private:
    const V &m_n;
    const V &m_l;
    const T &m_z;
    const T &m_r;
};

template <typename T, typename V>
inline CwiseNullaryOp<hydroR_functor<T, V>,
                      typename hydroR_functor<T, V>::ArrayType>
hydroR(const ArrayBase<V> &n,
       const ArrayBase<V> &l,
       const ArrayBase<T> &z,
       const ArrayBase<T> &r)
{
    static_assert(TYPE_IS(typename V::Scalar, int) ||
                      TYPE_IS(typename V::Scalar, unsigned int),
                  "n and l can only be int or uint array");
    eigen_assert(MATRIX_SAME_SIZE(n, l));
    eigen_assert(MATRIX_SAME_SIZE(n, z));
    eigen_assert(MATRIX_SAME_SIZE(n, r));

    typedef typename hydroR_functor<T, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(n.rows(),
                                  n.cols(),
                                  hydroR_functor<T, V>(n.derived(),
                                                       l.derived(),
                                                       z.derived(),
                                                       r.derived()));
}

template <typename T, typename V>
inline T hydroR_e_impl(const V n, const V l, const T z, const T r, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double hydroR_e_impl(
    const int n, const int l, const double z, const double r, double &e)
{
    gsl_sf_result ret;
    if (gsl_sf_hydrogenicR_e(n, l, z, r, &ret) == GSL_SUCCESS) {
        e = ret.err;
        return ret.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("hydroR"));
}

template <typename T, typename U, typename V>
class hydroR_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    hydroR_e_functor(const V &n, const V &l, const T &z, const T &r, U &e)
        : m_n(n)
        , m_l(l)
        , m_z(z)
        , m_r(r)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return hydroR_e_impl((int)m_n(i, j),
                             (int)m_l(i, j),
                             m_z(i, j),
                             m_r(i, j),
                             m_e(i, j));
    }

  private:
    const V &m_n;
    const V &m_l;
    const T &m_z;
    const T &m_r;
    U &m_e;
};

template <typename T, typename U, typename V>
inline CwiseNullaryOp<hydroR_e_functor<T, U, V>,
                      typename hydroR_e_functor<T, U, V>::ArrayType>
hydroR(const ArrayBase<V> &n,
       const ArrayBase<V> &l,
       const ArrayBase<T> &z,
       const ArrayBase<T> &r,
       ArrayBase<U> &e)
{
    static_assert(TYPE_IS(typename V::Scalar, int) ||
                      TYPE_IS(typename V::Scalar, unsigned int),
                  "n and l can only be int or uint array");
    eigen_assert(MATRIX_SAME_SIZE(n, l));
    eigen_assert(MATRIX_SAME_SIZE(n, z));
    eigen_assert(MATRIX_SAME_SIZE(n, r));

    typedef typename hydroR_e_functor<T, U, V>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(n.rows(),
                                  n.cols(),
                                  hydroR_e_functor<T, U, V>(n.derived(),
                                                            l.derived(),
                                                            z.derived(),
                                                            r.derived(),
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

#endif /* __IEXP_HYDROGENIC_R__ */
