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
    qawf_t(T epsabs, size_t limit)
        : m_epsabs(epsabs)
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

    T operator()(const typename unary_func<T>::type &fn,
                 qawo_table_t<T> &w,
                 T a,
                 void *opaque = nullptr,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T epsabs()
    {
        return m_epsabs;
    }

    qawf_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qawf_t &limit(T l)
    {
        if (m_limit != l) {
            m_limit = l;
            if (m_workspace != nullptr) {
                gsl_integration_workspace_free(m_workspace);
                m_workspace = gsl_integration_workspace_alloc(m_limit);
                IEXP_NOT_NULLPTR(m_workspace);
            }
        }
        return *this;
    }

  private:
    qawf_t(const qawf_t &) = delete;
    qawf_t &operator=(const qawf_t &) = delete;

    T m_epsabs;
    size_t m_limit;
    gsl_integration_workspace *m_workspace, *m_cycle_workspace;
};

template <>
double qawf_t<double>::operator()(const typename unary_func<double>::type &fn,
                                  qawo_table &w,
                                  double a,
                                  void *opaque,
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

    unary_func<double> m_fn(fn, opaque);
    double r, e;
    gsl_integration_qawf((gsl_function *)m_fn.gsl(),
                         a,
                         m_epsabs,
                         m_limit,
                         m_workspace,
                         m_cycle_workspace,
                         w.gsl(),
                         &r,
                         abserr != nullptr ? abserr : &e);
    return r;
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
