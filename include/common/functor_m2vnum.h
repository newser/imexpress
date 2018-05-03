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

#ifndef __IEXP_FUNCTOR_M2VNUM__
#define __IEXP_FUNCTOR_M2VNUM__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

// matrix to vector, keep count
#define M2VNUM_ROW(type, x) ((TP4(T) == RowMajor) ? (x).rows() : 1)

#define M2VNUM_COL(type, x) ((TP4(T) == RowMajor) ? 1 : (x).cols())

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T,
          typename _Scalar = TP1(T),
          int _Rows = ((TP4(T) == RowMajor) ? T::RowsAtCompileTime : 1),
          int _Cols = ((TP4(T) == RowMajor) ? 1 : T::ColsAtCompileTime),
          int _Options = TP4(T) & ~RowMajor,
          int _MaxRows = ((TP4(T) == RowMajor) ? T::MaxRowsAtCompileTime : 1),
          int _MaxCols = ((TP4(T) == RowMajor) ? 1 : T::MaxColsAtCompileTime)>
struct dense_derive_m2vnum
{
    using matrix = Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using array = Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using type = typename std::conditional<IS_MATRIX(T), matrix, array>::type;
};

// ========================================
// array parameter
// ========================================

template <typename Derived, typename T, typename R>
class functor_m2vnum_va
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vnum<T, R>::type;

    functor_m2vnum_va(const T &x)
        : m_result(new R[M2V_NUM(T, x)])
    {
        // it has to eval whole x and asseble a row or a colume into an array
        size_t num = M2V_NUM(T, x);
        size_t dim = M2V_DIM(T, x);
        typename type_eval<T>::type m_x(x.eval());
        for (size_t i = 0; i < num; ++i) {
            m_result.get()[i] = derived().m2vnum_va_impl(&m_x.data()[i * dim]);
        }
    }

    R operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<R> m_result;
};

template <typename Derived, typename T, typename U, typename R>
class functor_m2vnum_va_e
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vnum<T, R>::type;

    functor_m2vnum_va_e(const T &x, U &e)
        : m_result(new R[M2V_NUM(T, x)])
    {
        size_t num = M2V_NUM(T, x);
        size_t dim = M2V_DIM(T, x);
        typename type_eval<T>::type m_x(x.eval());
        eigen_assert(e.size() == num);
        for (size_t i = 0; i < num; ++i) {
            m_result.get()[i] =
                derived().m2vnum_e_impl(&m_x.data()[i * dim], e.data()[i]);
        }
    }

    R operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<R> m_result;
};

// ========================================
// fixed parameter
// ========================================

template <typename Derived, typename T, typename R, int dim>
class functor_m2vnum_fixed
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vnum<T, R>::type;

    functor_m2vnum_fixed(const T &x)
        : m_x(x)
    {
        static_assert(dim > 1, "dimension must be greater than 1");
        eigen_assert(M2V_DIM(T, x) == dim);
    }

    R operator()(Index i) const
    {
        return derived().compute(i, TYPE_BOOL(TP4(T) == RowMajor)());
    }

  protected:
    const T &m_x;
};

template <typename Derived, typename T, typename U, typename R, int dim>
class functor_m2vnum_fixed_e
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive_m2vnum<T, R>::type;

    functor_m2vnum_fixed_e(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
        static_assert(dim > 1, "dimension must be greater than 1");
        eigen_assert(M2V_DIM(T, x) == dim);
        eigen_assert(M2V_NUM(T, x) == e.size());
    }

    R operator()(Index i) const
    {
        return derived().compute(i,
                                 m_e.data()[i],
                                 TYPE_BOOL(TP4(T) == RowMajor)());
    }

  protected:
    const T &m_x;
    U &m_e;
};

// clang-format off
#define DEFINE_FUNCTOR_M2VNUM_ROWMAJOR(dim)                                    \
    /* dim parameters */                                                       \
    template <typename Derived, typename T, typename R>                        \
    class functor_m2vnum_##dim##d                                              \
        : public functor_m2vnum_fixed<functor_m2vnum_##dim##d<Derived, T, R>,  \
                                      T,                                       \
                                      R,                                       \
                                      dim>                                     \
    {                                                                          \
        DEFINE_DERIVED                                                         \
        using Base =                                                           \
            functor_m2vnum_fixed<functor_m2vnum_##dim##d<Derived, T, R>,       \
                                 T,                                            \
                                 R,                                            \
                                 dim>;                                         \
                                                                               \
      public:                                                                  \
        using ResultType = typename dense_derive_m2vnum<T, R>::type;           \
                                                                               \
        functor_m2vnum_##dim##d(const T &x)                                    \
            : Base(x)                                                          \
        {                                                                      \
        }                                                                      \
                                                                               \
        R compute(Index i, std::true_type) const                               \
        {                                                                      \
            /* row major */                                                    \
            return derived().m2vnum_impl(
#define DEFINE_FUNCTOR_M2VNUM_COLMAJOR(dim)                                    \
                                         );                                    \
        }                                                                      \
                                                                               \
        R compute(Index i, std::false_type) const                              \
        {                                                                      \
            /* col major */                                                    \
            return derived().m2vnum_impl(
#define DEFINE_FUNCTOR_M2VNUM_END(dim)                                         \
                                         );                                    \
        }                                                                      \
    };

#define DEFINE_FUNCTOR_M2VNUM_E_ROWMAJOR(dim)                                  \
    template <typename Derived, typename T, typename U, typename R>            \
    class functor_m2vnum_##dim##d_e                                            \
        : public functor_m2vnum_fixed_e<                                       \
              functor_m2vnum_##dim##d_e<Derived, T, U, R>,                     \
              T,                                                               \
              U,                                                               \
              R,                                                               \
              dim>                                                             \
    {                                                                          \
        DEFINE_DERIVED                                                         \
        using Base = functor_m2vnum_fixed_e<                                   \
            functor_m2vnum_##dim##d_e<Derived, T, U, R>,                       \
            T,                                                                 \
            U,                                                                 \
            R,                                                                 \
            dim>;                                                              \
                                                                               \
      public:                                                                  \
        using ResultType = typename dense_derive_m2vnum<T, R>::type;           \
                                                                               \
        functor_m2vnum_##dim##d_e(const T &x, U &e)                            \
            : Base(x, e)                                                       \
        {                                                                      \
        }                                                                      \
                                                                               \
        R compute(Index i, typename U::Scalar &e, std::true_type) const        \
        {                                                                      \
            /* row major */                                                    \
            return derived().m2vnum_e_impl(
#define DEFINE_FUNCTOR_M2VNUM_E_COLMAJOR(dim)                                  \
                                           , e);                               \
        }                                                                      \
                                                                               \
        R compute(Index i, typename U::Scalar &e, std::false_type) const       \
        {                                                                      \
            /* col major */                                                    \
            return derived().m2vnum_e_impl(
#define DEFINE_FUNCTOR_M2VNUM_E_END(dim)                                       \
                                           , e);                               \
        }                                                                      \
    };

// 2 parameters
DEFINE_FUNCTOR_M2VNUM_ROWMAJOR(2)
    Base::m_x(i, 0), Base::m_x(i, 1)
DEFINE_FUNCTOR_M2VNUM_COLMAJOR(2)
    Base::m_x(0, i), Base::m_x(1, i)
DEFINE_FUNCTOR_M2VNUM_END(2)

DEFINE_FUNCTOR_M2VNUM_E_ROWMAJOR(2)
    Base::m_x(i, 0), Base::m_x(i, 1)
DEFINE_FUNCTOR_M2VNUM_E_COLMAJOR(2)
    Base::m_x(0, i), Base::m_x(1, i)
DEFINE_FUNCTOR_M2VNUM_E_END(2)

// 3 parameters
DEFINE_FUNCTOR_M2VNUM_ROWMAJOR(3)
    Base::m_x(i, 0), Base::m_x(i, 1), Base::m_x(i, 2)
DEFINE_FUNCTOR_M2VNUM_COLMAJOR(3)
    Base::m_x(0, i), Base::m_x(1, i), Base::m_x(2, i)
DEFINE_FUNCTOR_M2VNUM_END(3)

DEFINE_FUNCTOR_M2VNUM_E_ROWMAJOR(3)
    Base::m_x(i, 0), Base::m_x(i, 1), Base::m_x(i, 2)
DEFINE_FUNCTOR_M2VNUM_E_COLMAJOR(3)
    Base::m_x(0, i), Base::m_x(1, i), Base::m_x(2, i)
DEFINE_FUNCTOR_M2VNUM_E_END(3)

// 4 parameters
DEFINE_FUNCTOR_M2VNUM_ROWMAJOR(4)
    Base::m_x(i, 0), Base::m_x(i, 1), Base::m_x(i, 2), Base::m_x(i, 3)
DEFINE_FUNCTOR_M2VNUM_COLMAJOR(4)
    Base::m_x(0, i), Base::m_x(1, i), Base::m_x(2, i), Base::m_x(3, i)
DEFINE_FUNCTOR_M2VNUM_END(4)

DEFINE_FUNCTOR_M2VNUM_E_ROWMAJOR(4)
    Base::m_x(i, 0), Base::m_x(i, 1), Base::m_x(i, 2), Base::m_x(i, 3)
DEFINE_FUNCTOR_M2VNUM_E_COLMAJOR(4)
    Base::m_x(0, i), Base::m_x(1, i), Base::m_x(2, i), Base::m_x(3, i)
DEFINE_FUNCTOR_M2VNUM_E_END(4)

    // clang-format on

    ////////////////////////////////////////////////////////////
    // global variants
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // interface declaration
    ////////////////////////////////////////////////////////////

    IEXP_NS_END

#endif /* __IEXP_FUNCTOR_M2VNUM__ */
