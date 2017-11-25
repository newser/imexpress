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

#ifndef __IEXP_DEBYE__
#define __IEXP_DEBYE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_debye.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

#define DEFINE_DEBYE(n)                                                        \
    template <typename T>                                                      \
    inline T debye##n##_impl(const T x)                                        \
    {                                                                          \
        throw std::invalid_argument("todo");                                   \
    }                                                                          \
                                                                               \
    template <>                                                                \
    inline double debye##n##_impl(const double x)                              \
    {                                                                          \
        return gsl_sf_debye_##n(x);                                            \
    }                                                                          \
                                                                               \
    template <typename T>                                                      \
    class debye##n##_functor                                                   \
    {                                                                          \
      public:                                                                  \
        typedef Array<typename T::Scalar,                                      \
                      T::RowsAtCompileTime,                                    \
                      T::ColsAtCompileTime,                                    \
                      T::Flags & RowMajorBit ? RowMajor : ColMajor,            \
                      T::MaxRowsAtCompileTime,                                 \
                      T::MaxColsAtCompileTime>                                 \
            ArrayType;                                                         \
                                                                               \
        debye##n##_functor(const T &x)                                         \
            : m_x(x)                                                           \
        {                                                                      \
        }                                                                      \
                                                                               \
        const typename T::Scalar operator()(Index i, Index j) const            \
        {                                                                      \
            return debye##n##_impl(m_x(i, j));                                 \
        }                                                                      \
                                                                               \
      private:                                                                 \
        const T &m_x;                                                          \
    };                                                                         \
                                                                               \
    template <typename T>                                                      \
    inline CwiseNullaryOp<debye##n##_functor<T>,                               \
                          typename debye##n##_functor<T>::ArrayType>           \
        debye##n(const ArrayBase<T> &x)                                        \
    {                                                                          \
        typedef typename debye##n##_functor<T>::ArrayType ArrayType;           \
        return ArrayType::NullaryExpr(x.rows(),                                \
                                      x.cols(),                                \
                                      debye##n##_functor<T>(x.derived()));     \
    }                                                                          \
                                                                               \
    template <typename T>                                                      \
    inline T debye##n##_e_impl(const T x, T &e)                                \
    {                                                                          \
        throw std::invalid_argument("todo");                                   \
    }                                                                          \
                                                                               \
    template <>                                                                \
    inline double debye##n##_e_impl(const double x, double &e)                 \
    {                                                                          \
        gsl_sf_result r;                                                       \
        if (gsl_sf_debye_##n##_e(x, &r) == GSL_SUCCESS) {                      \
            e = r.err;                                                         \
            return r.val;                                                      \
        }                                                                      \
        THROW_OR_RETURN_NAN(std::runtime_error("debye##n##"));                 \
    }                                                                          \
                                                                               \
    template <typename T, typename U>                                          \
    class debye##n##_e_functor                                                 \
    {                                                                          \
      public:                                                                  \
        typedef Array<typename T::Scalar,                                      \
                      T::RowsAtCompileTime,                                    \
                      T::ColsAtCompileTime,                                    \
                      T::Flags & RowMajorBit ? RowMajor : ColMajor,            \
                      T::MaxRowsAtCompileTime,                                 \
                      T::MaxColsAtCompileTime>                                 \
            ArrayType;                                                         \
                                                                               \
        debye##n##_e_functor(const T &x, U &e)                                 \
            : m_x(x)                                                           \
            , m_e(e)                                                           \
        {                                                                      \
        }                                                                      \
                                                                               \
        const typename T::Scalar operator()(Index i, Index j) const            \
        {                                                                      \
            return debye##n##_e_impl(m_x(i, j), m_e(i, j));                    \
        }                                                                      \
                                                                               \
      private:                                                                 \
        const T &m_x;                                                          \
        U &m_e;                                                                \
    };                                                                         \
                                                                               \
    template <typename T, typename U>                                          \
    inline CwiseNullaryOp<debye##n##_e_functor<T, U>,                          \
                          typename debye##n##_e_functor<T, U>::ArrayType>      \
        debye##n(const ArrayBase<T> &x, ArrayBase<U> &e)                       \
    {                                                                          \
        typedef typename debye##n##_e_functor<T, U>::ArrayType ArrayType;      \
        return ArrayType::NullaryExpr(x.rows(),                                \
                                      x.cols(),                                \
                                      debye##n##_e_functor<T,                  \
                                                           U>(x.derived(),     \
                                                              e.derived()));   \
    }

DEFINE_DEBYE(1)
DEFINE_DEBYE(2)
DEFINE_DEBYE(3)
DEFINE_DEBYE(4)
DEFINE_DEBYE(5)
DEFINE_DEBYE(6)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DEBYE__ */
