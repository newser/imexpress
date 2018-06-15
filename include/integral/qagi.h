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

#ifndef __IEXP_INTEGRAL_QAGI__
#define __IEXP_INTEGRAL_QAGI__

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

// ========================================
// qagi
// ========================================

template <typename T>
class qagi_t
{
  public:
    qagi_t(T epsabs, T epsrel, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qagi_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    T operator()(const typename unary_func<T>::type &fn,
                 void *opaque = nullptr,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T epsabs() const
    {
        return m_epsabs;
    }

    qagi_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qagi_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qagi_t &limit(T l)
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
    qagi_t(const qagi_t &) = delete;
    qagi_t &operator=(const qagi_t &) = delete;

    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qagi_t<double>::operator()(const typename unary_func<double>::type &fn,
                                  void *opaque,
                                  double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    unary_func<double> m_fn(fn, opaque);
    double r, e;
    gsl_integration_qagi((gsl_function *)m_fn.gsl(),
                         m_epsabs,
                         m_epsrel,
                         m_limit,
                         m_workspace,
                         &r,
                         abserr != nullptr ? abserr : &e);
    return r;
}

typedef qagi_t<double> qagi;

// ========================================
// qagiu
// ========================================

template <typename T>
class qagiu_t
{
  public:
    qagiu_t(T epsabs, T epsrel, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qagiu_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    double operator()(const typename unary_func<T>::type &fn,
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

    qagiu_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qagiu_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qagiu_t &limit(T l)
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

  protected:
    qagiu_t(const qagiu_t &) = delete;
    qagiu_t &operator=(const qagiu_t &) = delete;

    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qagiu_t<double>::operator()(const typename unary_func<double>::type &fn,
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
    gsl_integration_qagiu((gsl_function *)m_fn.gsl(),
                          a,
                          m_epsabs,
                          m_epsrel,
                          m_limit,
                          m_workspace,
                          &r,
                          abserr != nullptr ? abserr : &e);
    return r;
}

typedef qagiu_t<double> qagiu;

// ========================================
// qagil
// ========================================

template <typename T>
class qagil_t
{
  public:
    qagil_t(T epsabs, T epsrel, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qagil_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    double operator()(const typename unary_func<T>::type &fn,
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

    qagil_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qagil_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qagil_t &limit(T l)
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
    qagil_t(const qagil_t &) = delete;
    qagil_t &operator=(const qagil_t &) = delete;

    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qagil_t<double>::operator()(const typename unary_func<double>::type &fn,
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
    gsl_integration_qagil((gsl_function *)m_fn.gsl(),
                          b,
                          m_epsabs,
                          m_epsrel,
                          m_limit,
                          m_workspace,
                          &r,
                          abserr != nullptr ? abserr : &e);
    return r;
}

typedef qagil_t<double> qagil;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAGI__ */
