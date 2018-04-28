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

#ifndef __IEXP_SORT_INDEX__
#define __IEXP_SORT_INDEX__

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

template <typename T>
inline void sort_idx_impl(size_t *index, const T *data, size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void sort_idx_impl(size_t *index, const double *data, size_t n)
{
    gsl_sort_index(index, data, 1, n);
}

#define DEFINE_SORT_INDEX(type, suffix)                                        \
    template <>                                                                \
    inline void sort_idx_impl(size_t *index, const type *data, size_t n)       \
    {                                                                          \
        gsl_sort_##suffix##_index(index, data, 1, n);                          \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SORT_INDEX)
#undef DEFINE_SORT_INDEX

template <typename T, typename U>
class sort_idx_functor
{
  public:
    using ResultType = typename dense_derive<T, U>::type;

    sort_idx_functor(const T &x)
        : m_result(x.rows(), x.cols())
    {
        typename type_eval<T>::type m_v(x.eval());
        sort_idx_impl(m_result.data(), m_v.data(), m_v.size());
    }

    U operator()(Index i, Index j) const
    {
        return static_cast<U>(m_result(i, j));
    }

  private:
    buf_rc<size_t, bool(TP4(T) == RowMajor)> m_result;
};

template <typename U = size_t, typename T = void>
inline CwiseNullaryOp<sort_idx_functor<T, U>,
                      typename sort_idx_functor<T, U>::ResultType>
sort_idx(const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename sort_idx_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sort_idx_functor<T, U>(x.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_SORT_INDEX__ */
