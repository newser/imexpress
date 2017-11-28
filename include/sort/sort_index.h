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

#ifndef __IEXP_INDEX__
#define __IEXP_INDEX__

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
inline void sort_index_impl(size_t *index, const T *data, const size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void sort_index_impl(size_t *index, const double *data, const size_t n)
{
    return gsl_sort_index(index, data, 1, n);
}

#define DEFINE_SORT_INDEX(type, suffix)                                        \
    template <>                                                                \
    inline void sort_index_impl(size_t *index,                                 \
                                const type *data,                              \
                                const size_t n)                                \
    {                                                                          \
        return gsl_sort_##suffix##_index(index, data, 1, n);                   \
    }
DEFINE_TYPE_SUFFIX(DEFINE_SORT_INDEX)

template <typename T>
class sort_index_functor
{
  public:
    typedef Array<int,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;
    typedef Array<size_t,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        IndexArrayType;

    sort_index_functor(const T &v)
        : m_result(v.rows(), v.cols())
    {
        typename type_eval<T>::type m_v(v.eval());
        sort_index_impl(m_result.data(), m_v.data(), m_v.size());
    }

    const int operator()(Index i, Index j) const
    {
        return (int)m_result(i, j);
    }

  private:
    IndexArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<sort_index_functor<T>,
                      typename sort_index_functor<T>::ArrayType>
sort_index(const ArrayBase<T> &v)
{
    typedef typename sort_index_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(v.rows(),
                                  v.cols(),
                                  sort_index_functor<T>(v.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_INDEX__ */
