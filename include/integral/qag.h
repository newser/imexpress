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

#ifndef __IEXP_INTEGRAL_QAG__
#define __IEXP_INTEGRAL_QAG__

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
class qag_t
{
  public:
    enum class key
    {
        GAUSS15 = GSL_INTEG_GAUSS15,
        GAUSS21 = GSL_INTEG_GAUSS21,
        GAUSS31 = GSL_INTEG_GAUSS31,
        GAUSS41 = GSL_INTEG_GAUSS41,
        GAUSS51 = GSL_INTEG_GAUSS51,
        GAUSS61 = GSL_INTEG_GAUSS61,
    };

    qag_t(T epsabs, T epsrel, key k, size_t limit)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
        , m_key(k)
        , m_limit(limit)
        , m_workspace(nullptr)
    {
    }

    ~qag_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_workspace_free(m_workspace);
        }
    }

    T operator()(const typename unary_func<T>::type &fn,
                 T a,
                 T b,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T epsabs() const
    {
        return m_epsabs;
    }

    qag_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qag_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

    key key() const
    {
        return m_key;
    }

    qag_t &key(enum key k)
    {
        m_key = k;
        return *this;
    }

    T limit() const
    {
        return m_limit;
    }

    qag_t &limit(T l)
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
    qag_t(const qag_t &) = delete;
    qag_t &operator=(const qag_t &) = delete;

    T m_epsabs, m_epsrel;
    enum key m_key;
    size_t m_limit;
    gsl_integration_workspace *m_workspace;
};

template <>
double qag_t<double>::operator()(const typename unary_func<double>::type &fn,
                                 double a,
                                 double b,
                                 double *abserr)
{
    if (m_workspace == nullptr) {
        m_workspace = gsl_integration_workspace_alloc(m_limit);
        IEXP_NOT_NULLPTR(m_workspace);
    }

    unary_func<double> m_fn(fn);
    double r, e;
    gsl_integration_qag(m_fn.gsl(),
                        a,
                        b,
                        m_epsabs,
                        m_epsrel,
                        m_limit,
                        static_cast<int>(m_key),
                        m_workspace,
                        &r,
                        abserr != nullptr ? abserr : &e);
    return r;
}

typedef qag_t<double> qag;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QAG__ */
