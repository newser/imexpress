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
inline void largest_impl(T *dest, int k, const T *x, int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void largest_impl(double *dest, int k, const double *x, int n)
{
    gsl_sort_largest(dest, k, x, 1, n);
}

#define DEFINE_LARGEST(type, suffix)                                           \
    template <>                                                                \
    inline void largest_impl(type *dest, int k, const type *x, int n)          \
    {                                                                          \
        gsl_sort_##suffix##_largest(dest, k, x, 1, n);                         \
    }
DEFINE_TYPE_SUFFIX(DEFINE_LARGEST)

template <bool row_form, typename T>
class largest_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1,
                                             0,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    largest_functor(int k, const T &x)
        : m_result(new Scalar[k])
    {
        typename type_eval<T>::type m_x(x.eval());
        largest_impl(m_result.get(), k, m_x.data(), m_x.size());
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<largest_functor<row_form, T>,
                      typename largest_functor<row_form, T>::ResultType>
largest(int k, const DenseBase<T> &x)
{
    using ResultType = typename largest_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(row_form ? 1 : k,
                                   row_form ? k : 1,
                                   largest_functor<row_form, T>(k,
                                                                x.derived()));
}

// ========================================
// largest index
// ========================================

template <typename T>
inline void largest_index_impl(size_t *index, int k, const T *x, int n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void largest_index_impl(size_t *index, int k, const double *x, int n)
{
    gsl_sort_largest_index(index, k, x, 1, n);
}

#define DEFINE_LARGEST_INDEX(type, suffix)                                     \
    template <>                                                                \
    inline void largest_index_impl(size_t *index, int k, const type *x, int n) \
    {                                                                          \
        gsl_sort_##suffix##_largest_index(index, k, x, 1, n);                  \
    }
DEFINE_TYPE_SUFFIX(DEFINE_LARGEST_INDEX)
#undef DEFINE_LARGEST_INDEX

template <bool row_form, typename T, typename U>
class largest_index_functor
{
  public:
    using ResultType = typename dense_derive<T,
                                             U,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1,
                                             0,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    largest_index_functor(int k, const T &x)
        : m_result(new size_t[k])
    {
        typename type_eval<T>::type m_x(x.eval());
        largest_index_impl(m_result.get(), k, m_x.data(), m_x.size());
    }

    U operator()(Index i) const
    {
        // no overhead of casting integer
        return static_cast<U>(m_result.get()[i]);
    }

  private:
    // must be size_t array, not U
    std::shared_ptr<size_t> m_result;
};

template <typename U = size_t, bool row_form = false, typename T = void>
inline CwiseNullaryOp<
    largest_index_functor<row_form, T, U>,
    typename largest_index_functor<row_form, T, U>::ResultType>
largest_index(int k, const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(U), "only support integer index");

    using ResultType =
        typename largest_index_functor<row_form, T, U>::ResultType;
    return ResultType::
        NullaryExpr(row_form ? 1 : k,
                    row_form ? k : 1,
                    largest_index_functor<row_form, T, U>(k, x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_LARGEST__ */
