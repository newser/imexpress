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

#ifndef __IEXP_FUNCTOR_FOREACH__
#define __IEXP_FUNCTOR_FOREACH__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename Derived, typename T, typename R>
class functor_foreach
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive<T, R>::type;

    functor_foreach(const T &x)
        : m_x(x)
    {
    }

    R operator()(Index i, Index j) const
    {
        return derived().foreach_impl(m_x(i, j));
    }

  private:
    const T &m_x;
};

template <typename Derived, typename T, typename U, typename R>
class functor_foreach_e
{
    DEFINE_DERIVED

  public:
    using ResultType = typename dense_derive<T, R>::type;

    functor_foreach_e(const T &x, U &e)
        : m_x(x)
        , m_e(e)
    {
        eigen_assert(m_x.size() == m_e.size());
    }

    R operator()(Index i, Index j) const
    {
        return derived().foreach_e_impl(m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_x;
    U &m_e;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FUNCTOR_FOREACH__ */
