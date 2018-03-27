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

#ifndef __IEXP_INTEGRAL_FUNCTION__
#define __IEXP_INTEGRAL_FUNCTION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

IEXP_NS_BEGIN

namespace integral {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class unary_func
{
  public:
    using type = std::function<T(T)>;

    static T s_func(T x, void *param)
    {
        return ((type *)param)->operator()(x);
    }

    unary_func(const type &f)
        : m_fn(f)
        , m_gsl_fn{s_func,
                   const_cast<void *>(reinterpret_cast<const void *>(&m_fn))}
    {
        static_assert(TYPE_IS(T, double), "only support double now");
    }

    const gsl_function *gsl()
    {
        return &m_gsl_fn;
    }

  private:
    // m_fn can not be reference, msvc requires copying it
    const type m_fn;
    const gsl_function m_gsl_fn;
};

template <typename T>
class monte_func
{
  public:
    using type_1d = std::function<T(T)>;
    using type_nd = std::function<T(T *, size_t)>;

    static T s_func_1d(T *x, size_t dim, void *param)
    {
        return ((type_1d *)param)->operator()(x[0]);
    }

    static T s_func_nd(T *x, size_t dim, void *param)
    {
        return ((type_nd *)param)->operator()(x, dim);
    }

    monte_func(const type_1d &f)
        : m_fn(f)
        , m_gsl_fn{s_func_1d,
                   1,
                   const_cast<void *>(
                       reinterpret_cast<const void *>(&m_fn._1d))}
    {
        static_assert(TYPE_IS(T, double), "only support double now");
    }

    monte_func(size_t dim, const type_nd &f)
        : m_fn(f)
        , m_gsl_fn{s_func_nd,
                   dim,
                   const_cast<void *>(
                       reinterpret_cast<const void *>(&m_fn._nd))}
    {
        static_assert(TYPE_IS(T, double), "only support double now");

        // dim must be large than 1, or destructor does not work correctly
        eigen_assert(dim > 1);
    }

    ~monte_func()
    {
        if (m_gsl_fn.dim == 1) {
            m_fn._1d.~type_1d();
        } else {
            m_fn._nd.~type_nd();
        }
    }

    const gsl_monte_function *gsl()
    {
        return &m_gsl_fn;
    }

  private:
    // m_fn can not be reference, msvc requires copying it
    union fn
    {
        fn(const type_1d &f)
            : _1d(f)
        {
        }
        fn(const type_nd &f)
            : _nd(f)
        {
        }
        ~fn()
        {
        }

        const type_1d _1d;
        const type_nd _nd;
    } m_fn;
    const gsl_monte_function m_gsl_fn;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_FUNCTION__ */
