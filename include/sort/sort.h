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

#ifndef __IEXP_SORT__
#define __IEXP_SORT__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_sort.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// sort
// ========================================

template <typename T>
inline void sort_impl(T *data, const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void sort_impl(double *data, const int n)
{
    gsl_sort(data, 1, n);
}

#define DEFINE_SORT(type, suffix)                                              \
    template <>                                                                \
    inline void sort_impl(type *data, const int n)                             \
    {                                                                          \
        gsl_sort_##suffix(data, 1, n);                                         \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SORT)

template <typename T>
class sort_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sort_functor(const T &v)
        : m_result(v.eval())
    {
        sort_impl(m_result.data(), v.size());
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<sort_functor<T>, typename sort_functor<T>::ArrayType>
sort(const ArrayBase<T> &v)
{
    typedef typename sort_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(v.rows(),
                                  v.cols(),
                                  sort_functor<T>(v.derived()));
}

// ========================================
// sort2
// ========================================

template <typename T>
inline void sort2_impl(T *data1, T *data2, const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void sort2_impl(double *data1, double *data2, const int n)
{
    gsl_sort2(data1, 1, data2, 1, n);
}

#define DEFINE_SORT2(type, suffix)                                             \
    template <>                                                                \
    inline void sort2_impl(type *data1, type *data2, const int n)              \
    {                                                                          \
        gsl_sort2_##suffix(data1, 1, data2, 1, n);                             \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SORT2)

template <typename T>
class sort2_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    sort2_functor(const T &v1, T &v2)
        : m_result(v1.eval())
    {
        sort2_impl(m_result.data(), v2.data(), v1.size());
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<sort2_functor<T>, typename sort2_functor<T>::ArrayType>
sort(const ArrayBase<T> &v1, ArrayBase<T> &v2)
{
    eigen_assert(MATRIX_SAME_SIZE(v1, v2));

    typedef typename sort2_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(v1.rows(),
                                  v1.cols(),
                                  sort2_functor<T>(v1.derived(), v2.derived()));
}

// ========================================
// sort inplace
// ========================================

template <typename T>
void sort_inplace(T &v)
{
    sort_impl(v.data(), v.size());
}

template <typename T>
void sort_inplace(T &v1, T &v2)
{
    eigen_assert(MATRIX_SAME_SIZE(v1, v2));
    sort2_impl(v1.data(), v2.data(), v1.size());
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_SORT__ */
