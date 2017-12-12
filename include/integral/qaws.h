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
class qaws_table_t
{
  public:
    qaws_table_t(const T alpha, const T beta, const int mu, const int nu)
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

    gsl_integration_qaws_table *get()
    {
        if (m_table == nullptr) {
            m_table =
                gsl_integration_qaws_table_alloc(m_alpha, m_beta, m_mu, m_nu);
            IEXP_NOT_NULLPTR(m_table);
        }
        return m_table;
    }

    void set(const T alpha, const T beta, const int mu, const int nu)
    {
        if (m_table == nullptr) {
            m_table = gsl_integration_qaws_table_alloc(alpha, beta, mu, nu);
            IEXP_NOT_NULLPTR(m_table);
        } else {
            eigen_assert(
                gsl_integration_qaws_table_set(m_table, alpha, beta, mu, nu) ==
                GSL_SUCCESS);
        }
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
    qaws_t(const typename unary_func<T>::type &fn,
           T epsabs,
           T epsrel,
           size_t limit)
        : m_fn(fn)
        , m_epsabs(epsabs)
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

    int operator()(qaws_table_t<T> &w,
                   const T a,
                   const T b,
                   T *result,
                   T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T &epsabs()
    {
        return m_epsabs;
    }

    T &epsrel()
    {
        return m_epsrel;
    }

  private:
    qaws_t(const qaws_t &) = delete;
    qaws_t &operator=(const qaws_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qaws_t<double>::operator()(qaws_table_t<double> &w,
                               const double a,
                               const double b,
                               double *result,
                               double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qaws((gsl_function *)m_fn.gsl(),
                                a,
                                b,
                                w.get(),
                                m_epsabs,
                                m_epsrel,
                                m_limit,
                                m_workspace,
                                result,
                                abserr != nullptr ? abserr : &__abserr);
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
