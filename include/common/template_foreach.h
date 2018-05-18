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

#ifndef __IEXP_TEMPLATE_FOREACH__
#define __IEXP_TEMPLATE_FOREACH__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

#define DEFINE_TEMPLATE_FOREACH(name, rt, p1t, gf)                             \
    template <typename T>                                                      \
    class name##_functor : public functor_foreach<name##_functor<T>, T, rt>    \
    {                                                                          \
      public:                                                                  \
        name##_functor(const T &x)                                             \
            : functor_foreach<name##_functor<T>, T, rt>(x)                     \
        {                                                                      \
        }                                                                      \
                                                                               \
        rt foreach_impl(p1t x) const                                           \
        {                                                                      \
            return gf(x);                                                      \
        }                                                                      \
    };                                                                         \
                                                                               \
    template <typename T>                                                      \
    inline CwiseNullaryOp<name##_functor<T>,                                   \
                          typename name##_functor<T>::ResultType>              \
    name(const DenseBase<T> &x)                                                \
    {                                                                          \
        using ResultType = typename name##_functor<T>::ResultType;             \
        return ResultType::NullaryExpr(x.rows(),                               \
                                       x.cols(),                               \
                                       name##_functor<T>(x.derived()));        \
    }

#define DEFINE_TEMPLATE_FOREACH_E(name, rt, p1t, gfe)                          \
    template <typename T, typename U>                                          \
    class name##_e_functor                                                     \
        : public functor_foreach_e<name##_e_functor<T, U>, T, U, rt>           \
    {                                                                          \
      public:                                                                  \
        name##_e_functor(const T &x, U &e)                                     \
            : functor_foreach_e<name##_e_functor<T, U>, T, U, rt>(x, e)        \
        {                                                                      \
        }                                                                      \
                                                                               \
        rt foreach_e_impl(p1t x, p1t &e) const                                 \
        {                                                                      \
            gsl_sf_result r;                                                   \
            gfe(x, &r);                                                        \
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
        using ResultType = typename name##_e_functor<T, U>::ResultType;        \
        return ResultType::NullaryExpr(x.rows(),                               \
                                       x.cols(),                               \
                                       name##_e_functor<T, U>(x.derived(),     \
                                                              e.derived()));   \
    }

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

#endif /* __IEXP_TEMPLATE_FOREACH__ */
