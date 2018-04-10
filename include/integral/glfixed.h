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

#ifndef __IEXP_INTEGRAL_GLFIXED__
#define __IEXP_INTEGRAL_GLFIXED__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>
#include <integral/type.h>

#include <gsl/gsl_integration.h>

IEXP_NS_BEGIN

namespace integral {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class glfixed_t
{
  public:
    glfixed_t(size_t n)
        : m_n(n)
        , m_table(nullptr)
    {
    }

    ~glfixed_t()
    {
        if (m_table != nullptr) {
            gsl_integration_glfixed_table_free(m_table);
        }
    }

    T operator()(const typename unary_func<T>::type &fn, T a, T b)
    {
        UNSUPPORTED_TYPE(T);
    }

    T n() const
    {
        return m_n;
    }

    glfixed_t &n(T n)
    {
        if (m_n != n) {
            m_n = n;
            if (m_table != nullptr) {
                gsl_integration_glfixed_table_free(m_table);
                m_table = gsl_integration_glfixed_table_alloc(m_n);
                IEXP_NOT_NULLPTR(m_table);
            }
        }
        return *this;
    }

  private:
    glfixed_t(const glfixed_t &) = delete;
    glfixed_t &operator=(const glfixed_t &) = delete;

    size_t m_n;
    gsl_integration_glfixed_table *m_table;
};

template <>
double glfixed_t<double>::operator()(
    const typename unary_func<double>::type &fn, double a, double b)
{
    if (m_table == nullptr) {
        m_table = gsl_integration_glfixed_table_alloc(m_n);
        IEXP_NOT_NULLPTR(m_table);
    }

    unary_func<double> m_fn(fn);
    return gsl_integration_glfixed(m_fn.gsl(), a, b, m_table);
}

typedef glfixed_t<double> glfixed;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_GLFIXED__ */
