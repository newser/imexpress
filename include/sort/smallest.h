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

#ifndef __IEXP_SMALLEST__
#define __IEXP_SMALLEST__

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
// smallest
// ========================================

template <typename T>
inline void smallest_impl(T *dest, const int k, const T *src, const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void smallest_impl(double *dest,
                          const int k,
                          const double *src,
                          const int n)
{
    gsl_sort_smallest(dest, k, src, 1, n);
}

#define DEFINE_SMALLEST(type, suffix)                                          \
    template <>                                                                \
    inline void smallest_impl(type *dest,                                      \
                              const int k,                                     \
                              const type *src,                                 \
                              const int n)                                     \
    {                                                                          \
        gsl_sort_##suffix##_smallest(dest, k, src, 1, n);                      \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SMALLEST)

template <typename T>
class smallest_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        ArrayType;

    smallest_functor(const int k, const T &v)
        : m_result(k)
    {
        typename type_eval<T>::type m_v(v.eval());
        smallest_impl(m_result.data(), k, m_v.data(), m_v.size());
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<smallest_functor<T>,
                      typename smallest_functor<T>::ArrayType>
smallest(const int k, const ArrayBase<T> &v)
{
    eigen_assert(IS_VEC(v));

    typedef typename smallest_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(k, 1, smallest_functor<T>(k, v.derived()));
}

// ========================================
// smallest index
// ========================================

template <typename T>
inline void smallest_index_impl(size_t *index,
                                const int k,
                                const T *src,
                                const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void smallest_index_impl(size_t *index,
                                const int k,
                                const double *src,
                                const int n)
{
    gsl_sort_smallest_index(index, k, src, 1, n);
}

#define DEFINE_SMALLEST_INDEX(type, suffix)                                    \
    template <>                                                                \
    inline void smallest_index_impl(size_t *index,                             \
                                    const int k,                               \
                                    const type *src,                           \
                                    const int n)                               \
    {                                                                          \
        gsl_sort_##suffix##_smallest_index(index, k, src, 1, n);               \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SMALLEST_INDEX)

template <typename T>
class smallest_index_functor
{
  public:
    typedef Array<int,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        ArrayType;
    typedef Array<size_t,
                  T::SizeAtCompileTime,
                  1,
                  ColMajor,
                  T::SizeAtCompileTime,
                  1>
        IndexArrayType;

    smallest_index_functor(const int k, const T &v)
        : m_result(k)
    {
        typename type_eval<T>::type m_v(v.eval());
        smallest_index_impl(m_result.data(), k, m_v.data(), m_v.size());
    }

    const int operator()(Index i, Index j) const
    {
        return (int)m_result(i, j);
    }

  private:
    IndexArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<smallest_index_functor<T>,
                      typename smallest_index_functor<T>::ArrayType>
smallest_index(const int k, const ArrayBase<T> &v)
{
    eigen_assert(IS_VEC(v));

    typedef typename smallest_index_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(k,
                                  1,
                                  smallest_index_functor<T>(k, v.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_SMALLEST__ */
