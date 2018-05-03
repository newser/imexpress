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

#ifndef __IEXP_FUNCTOR_M2VDIM__
#define __IEXP_FUNCTOR_M2VDIM__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

// matrix to vector, keep dimension
#define M2VDIM_ROW(type, x) ((TP4(T) == RowMajor) ? 1 : (x).rows())

#define M2VDIM_COL(type, x) ((TP4(T) == RowMajor) ? (x).cols() : 1)

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T,
          typename _Scalar = TP1(T),
          int _Rows = ((TP4(T) == RowMajor) ? 1 : T::RowsAtCompileTime),
          int _Cols = ((TP4(T) == RowMajor) ? T::ColsAtCompileTime : 1),
          int _Options = TP4(T) & ~RowMajor,
          int _MaxRows = ((TP4(T) == RowMajor) ? 1 : T::MaxRowsAtCompileTime),
          int _MaxCols = ((TP4(T) == RowMajor) ? T::MaxColsAtCompileTime : 1)>
struct dense_derive_m2vdim
{
    using matrix = Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using array = Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using type = typename std::conditional<IS_MATRIX(T), matrix, array>::type;
};

// ========================================
// array parameter
// ========================================

template <typename Derived, typename T, typename R>
class functor_m2vdim_va
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vdim<T, R>::type;

    functor_m2vdim_va(const T &x)
        : m_result(new R[M2V_DIM(T, x)])
    {
        // it has to eval whole x and then pass to m2vdim_va_impl
        size_t num = M2V_NUM(T, x);
        size_t dim = M2V_DIM(T, x);
        typename type_eval<T>::type m_x(x.eval());
        derived().m2vdim_va_impl(m_x.data(), num, dim, m_result.get());
    }

    R operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<R> m_result;
};

template <typename Derived, typename T, typename U, typename R>
class functor_m2vdim_va_e
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vdim<T, R>::type;

    functor_m2vdim_va_e(const T &x, U &e)
        : m_result(new R[M2V_DIM(T, x)])
    {
        size_t num = M2V_NUM(T, x);
        size_t dim = M2V_DIM(T, x);
        typename type_eval<T>::type m_x(x.eval());
        eigen_assert(e.size() == num);
        derived().m2vdim_va_impl(m_x.data(),
                                 num,
                                 dim,
                                 m_result.get(),
                                 e.data());
    }

    R operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<R> m_result;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FUNCTOR_M2VDIM__ */
