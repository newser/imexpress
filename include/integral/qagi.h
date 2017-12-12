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

// ========================================
// qagi
// ========================================

template <typename T>
class qagi_t
{
  public:
    qagi_t(const typename unary_func<T>::type &fn,
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

    ~qagi_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    int operator()(T *result, T *abserr = nullptr)
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
    qagi_t(const qagi_t &) = delete;
    qagi_t &operator=(const qagi_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qagi_t<double>::operator()(double *result, double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qagi((gsl_function *)m_fn.gsl(),
                                m_epsabs,
                                m_epsrel,
                                m_limit,
                                m_workspace,
                                result,
                                abserr != nullptr ? abserr : &__abserr);
}

typedef qagi_t<double> qagi;

// ========================================
// qagiu
// ========================================

template <typename T>
class qagiu_t
{
  public:
    qagiu_t(const typename unary_func<T>::type &fn,
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

    ~qagiu_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    int operator()(const T a, T *result, T *abserr = nullptr)
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
    qagiu_t(const qagiu_t &) = delete;
    qagiu_t &operator=(const qagiu_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qagiu_t<double>::operator()(const double a, double *result, double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qagiu((gsl_function *)m_fn.gsl(),
                                 a,
                                 m_epsabs,
                                 m_epsrel,
                                 m_limit,
                                 m_workspace,
                                 result,
                                 abserr != nullptr ? abserr : &__abserr);
}

typedef qagiu_t<double> qagiu;

// ========================================
// qagil
// ========================================

template <typename T>
class qagil_t
{
  public:
    qagil_t(const typename unary_func<T>::type &fn,
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

    ~qagil_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    int operator()(const T b, T *result, T *abserr = nullptr)
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
    qagil_t(const qagil_t &) = delete;
    qagil_t &operator=(const qagil_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qagil_t<double>::operator()(const double b, double *result, double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qagil((gsl_function *)m_fn.gsl(),
                                 b,
                                 m_epsabs,
                                 m_epsrel,
                                 m_limit,
                                 m_workspace,
                                 result,
                                 abserr != nullptr ? abserr : &__abserr);
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
