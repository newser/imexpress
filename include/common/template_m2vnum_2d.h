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

#ifndef __IEXP_TEMPLATE_M2VNUM_2D__
#define __IEXP_TEMPLATE_M2VNUM_2D__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

#define DEFINE_TEMPLATE_M2VNUM_2D(name, return_t, scalar_t, p1t, p2t, gf)      \
    template <typename T>                                                      \
    class name##_functor                                                       \
        : public functor_m2vnum_2d<name##_functor<T>, T, return_t>             \
    {                                                                          \
      public:                                                                  \
        name##_functor(const T &x)                                             \
            : functor_m2vnum_2d<name##_functor<T>, T, return_t>(x)             \
        {                                                                      \
        }                                                                      \
                                                                               \
        return_t m2vnum_impl(scalar_t x1, scalar_t x2) const                   \
        {                                                                      \
            return gf((p1t)x1, (p2t)x2);                                       \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <typename T>                                                      \
    inline CwiseNullaryOp<name##_functor<T>,                                   \
                          typename name##_functor<T>::ResultType>              \
    name(const DenseBase<T> &x)                                                \
    {                                                                          \
        static_assert(TYPE_IS(typename T::Scalar, scalar_t),                   \
                      "only support " #scalar_t " scalar");                    \
                                                                               \
        using ResultType = typename name##_functor<T>::ResultType;             \
        return ResultType::NullaryExpr(M2VNUM_ROW(T, x),                       \
                                       M2VNUM_COL(T, x),                       \
                                       name##_functor<T>(x.derived()));        \
    }

#define DEFINE_TEMPLATE_M2VNUM_2D_E(name, return_t, scalar_t, p1t, p2t, gfe)   \
    template <typename T, typename U>                                          \
    class name##_e_functor                                                     \
        : public functor_m2vnum_2d_e<name##_e_functor<T, U>, T, U, return_t>   \
    {                                                                          \
      public:                                                                  \
        name##_e_functor(const T &x, U &e)                                     \
            : functor_m2vnum_2d_e<name##_e_functor<T, U>, T, U, return_t>(x,   \
                                                                          e)   \
        {                                                                      \
        }                                                                      \
                                                                               \
        return_t m2vnum_e_impl(scalar_t x1, scalar_t x2, double &e) const      \
        {                                                                      \
            gsl_sf_result r;                                                   \
            gfe((p1t)x1, (p2t)x2, &r);                                         \
            e = r.err;                                                         \
            return r.val;                                                      \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <typename T, typename U>                                          \
    inline CwiseNullaryOp<name##_e_functor<T, U>,                              \
                          typename name##_e_functor<T, U>::ResultType>         \
    name(const DenseBase<T> &x, DenseBase<U> &e)                               \
    {                                                                          \
        static_assert(TYPE_IS(typename T::Scalar, scalar_t),                   \
                      "only support " #scalar_t " scalar");                    \
                                                                               \
        using ResultType = typename name##_e_functor<T, U>::ResultType;        \
        return ResultType::NullaryExpr(M2VNUM_ROW(T, x),                       \
                                       M2VNUM_COL(T, x),                       \
                                       name##_e_functor<T, U>(x.derived(),     \
                                                              e.derived()));   \
    }

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

#endif /* __IEXP_TEMPLATE_M2VNUM_2D__ */
