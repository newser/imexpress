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
class qawo_table_t
{
  public:
    qawo_table_t(const T omega, const T length, const bool sine, const size_t n)
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

    gsl_integration_qawo_table *get()
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

    void set(const T omega, const T length, const bool sine)
    {
        if (m_table != nullptr) {
            eigen_assert(
                gsl_integration_qawo_table_set(m_table,
                                               omega,
                                               length,
                                               sine ? GSL_INTEG_SINE
                                                    : GSL_INTEG_COSINE) ==
                GSL_SUCCESS);
        }

        m_omega = omega;
        m_length = length;
        m_sine = sine;
    }

    void set_length(const T length)
    {
        if (m_table != nullptr) {
            eigen_assert(
                gsl_integration_qawo_table_set_length(m_table, length) ==
                GSL_SUCCESS);
        }

        m_length = length;
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
    qawo_t(const typename unary_func<T>::type &fn,
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

    ~qawo_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
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

    T &epsrel()
    {
        return m_epsrel;
    }

  private:
    qawo_t(const qawo_t &) = delete;
    qawo_t &operator=(const qawo_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qawo_t<double>::operator()(qawo_table_t<double> &w,
                               const double a,
                               double *result,
                               double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qawo((gsl_function *)m_fn.gsl(),
                                a,
                                m_epsabs,
                                m_epsrel,
                                m_limit,
                                m_workspace,
                                w.get(),
                                result,
                                abserr != nullptr ? abserr : &__abserr);
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
