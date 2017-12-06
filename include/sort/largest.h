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

#ifndef __IEXP_LARGEST__
#define __IEXP_LARGEST__

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
// largest
// ========================================

template <typename T>
inline void largest_impl(T *dest, const int k, const T *src, const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void largest_impl(double *dest,
                         const int k,
                         const double *src,
                         const int n)
{
    gsl_sort_largest(dest, k, src, 1, n);
}

#define DEFINE_LARGEST(type, suffix)                                           \
    template <>                                                                \
    inline void largest_impl(type *dest,                                       \
                             const int k,                                      \
                             const type *src,                                  \
                             const int n)                                      \
    {                                                                          \
        gsl_sort_##suffix##_largest(dest, k, src, 1, n);                       \
    }
DEFINE_TYPE_SUFFIX(DEFINE_LARGEST)

template <typename T>
class largest_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::SizeAtCompileTime,
                            1,
                            ColMajor,
                            T::SizeAtCompileTime,
                            1>;

    largest_functor(const int k, const T &v)
        : m_result(k)
    {
        typename type_eval<T>::type m_v(v.eval());
        largest_impl(m_result.data(), k, m_v.data(), m_v.size());
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return m_result(i, j);
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<largest_functor<T>,
                      typename largest_functor<T>::ArrayType>
largest(const int k, const ArrayBase<T> &v)
{
    eigen_assert(IS_VEC(v));

    using ArrayType = typename largest_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(k, 1, largest_functor<T>(k, v.derived()));
}

// ========================================
// largest index
// ========================================

template <typename T>
inline void largest_index_impl(size_t *index,
                               const int k,
                               const T *src,
                               const int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void largest_index_impl(size_t *index,
                               const int k,
                               const double *src,
                               const int n)
{
    gsl_sort_largest_index(index, k, src, 1, n);
}

#define DEFINE_LARGEST_INDEX(type, suffix)                                     \
    template <>                                                                \
    inline void largest_index_impl(size_t *index,                              \
                                   const int k,                                \
                                   const type *src,                            \
                                   const int n)                                \
    {                                                                          \
        gsl_sort_##suffix##_largest_index(index, k, src, 1, n);                \
    }
DEFINE_TYPE_SUFFIX(DEFINE_LARGEST_INDEX)
#undef DEFINE_LARGEST_INDEX

template <typename T>
class largest_index_functor
{
  public:
    using ArrayType =
        Array<int, T::SizeAtCompileTime, 1, ColMajor, T::SizeAtCompileTime, 1>;
    using IndexArrayType = Array<size_t,
                                 T::SizeAtCompileTime,
                                 1,
                                 ColMajor,
                                 T::SizeAtCompileTime,
                                 1>;

    largest_index_functor(const int k, const T &v)
        : m_result(k)
    {
        typename type_eval<T>::type m_v(v.eval());
        largest_index_impl(m_result.data(), k, m_v.data(), m_v.size());
    }

    const int operator()(Index i, Index j) const
    {
        return (int)m_result(i, j);
    }

  private:
    IndexArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<largest_index_functor<T>,
                      typename largest_index_functor<T>::ArrayType>
largest_index(const int k, const ArrayBase<T> &v)
{
    eigen_assert(IS_VEC(v));

    using ArrayType = typename largest_index_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(k,
                                  1,
                                  largest_index_functor<T>(k, v.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_LARGEST__ */
