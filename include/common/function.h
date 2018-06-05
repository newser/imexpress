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

#ifndef __IEXP_FUNCTION__
#define __IEXP_FUNCTION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>

IEXP_NS_BEGIN

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

    T operator()(T x)
    {
        return m_fn(x);
    }

    const gsl_function *gsl() const
    {
        return &m_gsl_fn;
    }

  private:
    // m_fn can not be reference, msvc requires copying it
    const type m_fn;
    const gsl_function m_gsl_fn;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FUNCTION__ */
