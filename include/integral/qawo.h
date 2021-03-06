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

#ifndef __IEXP_INTEGRAL_QAWO__
#define __IEXP_INTEGRAL_QAWO__

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
class qawo_table_t
{
  public:
    qawo_table_t(T omega, T length, bool sine, size_t n)
        : m_omega(omega)
        , m_length(length)
        , m_sine(sine)
        , m_n(n)
        , m_table(nullptr)
    {
    }

    ~qawo_table_t()
    {
        if (m_table != nullptr) {
            gsl_integration_qawo_table_free(m_table);
        }
    }

    T omega() const
    {
        return m_omega;
    }

    qawo_table_t &omega(T o)
    {
        m_omega = o;
        return *this;
    }

    T length() const
    {
        return m_length;
    }

    qawo_table_t &length(T len)
    {
        m_length = len;
        return *this;
    }

    bool sine() const
    {
        return m_sine;
    }

    qawo_table_t &sine(bool s)
    {
        m_sine = s;
        return *this;
    }

    size_t n() const
    {
        return m_n;
    }

    qawo_table_t &n(size_t s)
    {
        m_n = s;
        return *this;
    }

    void reload()
    {
        if (m_table != nullptr) {
            gsl_integration_qawo_table_free(m_table);
            m_table = nullptr;
        }

        gsl();
    }

    gsl_integration_qawo_table *gsl()
    {
        if (m_table == nullptr) {
            m_table =
                gsl_integration_qawo_table_alloc(m_omega,
                                                 m_length,
                                                 m_sine ? GSL_INTEG_SINE
                                                        : GSL_INTEG_COSINE,
                                                 m_n);
            IEXP_NOT_NULLPTR(m_table);
        }
        return m_table;
    }

  private:
    T m_omega, m_length;
    bool m_sine;
    size_t m_n;
    gsl_integration_qawo_table *m_table;
};

typedef qawo_table_t<double> qawo_table;

template <typename T>
class qawo_t
{
  public:
    qawo_t(T epsabs, T epsrel, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qawo_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
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

    T epsabs() const
    {
        return m_epsabs;
    }

    qawo_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qawo_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qawo_t &limit(T l)
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
    qawo_t(const qawo_t &) = delete;
    qawo_t &operator=(const qawo_t &) = delete;

    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qawo_t<double>::operator()(const typename unary_func<double>::type &fn,
                                  qawo_table_t<double> &w,
                                  double a,
                                  void *opaque,
                                  double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    unary_func<double> m_fn(fn, opaque);
    double r, e;
    gsl_integration_qawo((gsl_function *)m_fn.gsl(),
                         a,
                         m_epsabs,
                         m_epsrel,
                         m_limit,
                         m_workspace,
                         w.gsl(),
                         &r,
                         abserr != nullptr ? abserr : &e);
    return r;
}

typedef qawo_t<double> qawo;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAWO__ */
