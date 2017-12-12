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

#ifndef __IEXP_INTEGRAL_QAGS__
#define __IEXP_INTEGRAL_QAGS__

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
class qags_t
{
  public:
    qags_t(const typename unary_func<T>::type &fn,
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

    ~qags_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    int operator()(const T a, const T b, T *result, T *abserr = nullptr)
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
    qags_t(const qags_t &) = delete;
    qags_t &operator=(const qags_t &) = delete;

    unary_func<T> m_fn;
    T m_epsabs, m_epsrel;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
int qags_t<double>::operator()(const double a,
                               const double b,
                               double *result,
                               double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    double __abserr;
    return gsl_integration_qags(m_fn.gsl(),
                                a,
                                b,
                                m_epsabs,
                                m_epsrel,
                                m_limit,
                                m_workspace,
                                result,
                                abserr != nullptr ? abserr : &__abserr);
}

typedef qags_t<double> qags;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAGS__ */
