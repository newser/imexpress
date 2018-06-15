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

#ifndef __IEXP_INTEGRAL_QAWS__
#define __IEXP_INTEGRAL_QAWS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>

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
class qaws_table_t
{
  public:
    qaws_table_t(T alpha, T beta, int mu, int nu)
        : m_alpha(alpha)
        , m_beta(beta)
        , m_mu(mu)
        , m_nu(nu)
        , m_table(nullptr)
    {
    }

    ~qaws_table_t()
    {
        if (m_table != nullptr) {
            gsl_integration_qaws_table_free(m_table);
        }
    }

    T alpha()
    {
        return m_alpha;
    }

    qaws_table_t &alpha(T e)
    {
        m_alpha = e;
        return *this;
    }

    T beta()
    {
        return m_beta;
    }

    qaws_table_t &beta(T e)
    {
        m_beta = e;
        return *this;
    }

    T mu()
    {
        return m_mu;
    }

    qaws_table_t &mu(T e)
    {
        m_mu = e;
        return *this;
    }

    T nu()
    {
        return m_nu;
    }

    qaws_table_t &nu(T e)
    {
        m_nu = e;
        return *this;
    }

    void reload()
    {
        if (m_table != nullptr) {
            gsl_integration_qaws_table_free(m_table);
            m_table = nullptr;
        }

        gsl();
    }

    gsl_integration_qaws_table *gsl()
    {
        if (m_table == nullptr) {
            m_table =
                gsl_integration_qaws_table_alloc(m_alpha, m_beta, m_mu, m_nu);
            IEXP_NOT_NULLPTR(m_table);
        }
        return m_table;
    }

  private:
    T m_alpha, m_beta;
    int m_mu, m_nu;
    gsl_integration_qaws_table *m_table;
};

typedef qaws_table_t<double> qaws_table;

template <typename T>
class qaws_t
{
  public:
    qaws_t(T epsabs, T epsrel, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qaws_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    T operator()(const typename unary_func<T>::type &fn,
                 qaws_table_t<T> &w,
                 T a,
                 T b,
                 void *opaque = nullptr,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T epsabs() const
    {
        return m_epsabs;
    }

    qaws_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qaws_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qaws_t &limit(T l)
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
    qaws_t(const qaws_t &) = delete;
    qaws_t &operator=(const qaws_t &) = delete;

    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qaws_t<double>::operator()(const typename unary_func<double>::type &fn,
                                  qaws_table_t<double> &w,
                                  double a,
                                  double b,
                                  void *opaque,
                                  double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    unary_func<double> m_fn(fn, opaque);
    double r, e;
    gsl_integration_qaws((gsl_function *)m_fn.gsl(),
                         a,
                         b,
                         w.gsl(),
                         m_epsabs,
                         m_epsrel,
                         m_limit,
                         m_workspace,
                         &r,
                         abserr != nullptr ? abserr : &e);
    return r;
}

typedef qaws_t<double> qaws;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAWS__ */
