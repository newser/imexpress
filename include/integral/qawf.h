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

#ifndef __IEXP_INTEGRAL_QAWF__
#define __IEXP_INTEGRAL_QAWF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>
#include <integral/qawo.h>
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
class qawf_t
{
  public:
    qawf_t(const typename unary_func<T>::type &fn, T epsabs, size_t limit)
        : m_fn(fn)
        , m_epsabs(epsabs)
        , m_limit(limit)
        , m_workspace(nullptr)
        , m_cycle_workspace(nullptr)
    {
    }

    ~qawf_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }

        if (m_cycle_workspace != nullptr) {
            gsl_integration_workspace_free(m_cycle_workspace);
        }
    }

    int operator()(qawo_table_t<T> &w,
                   const T a,
                   T *result,
                   T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T &epsabs()
    {
        return m_epsabs;
    }

  private:
    qawf_t(const qawf_t &) = delete;
    qawf_t &operator=(const qawf_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs;
    size_t m_limit;
    gsl_integration_workspace *m_workspace, *m_cycle_workspace;
};

template <>
int qawf_t<double>::operator()(qawo_table &w,
                               const double a,
                               double *result,
                               double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    if (m_cycle_workspace == nullptr) {
        m_cycle_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qawf((gsl_function *)m_fn.gsl(),
                                a,
                                m_epsabs,
                                m_limit,
                                m_workspace,
                                m_cycle_workspace,
                                w.get(),
                                result,
                                abserr != nullptr ? abserr : &__abserr);
}

typedef qawf_t<double> qawf;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAWF__ */
