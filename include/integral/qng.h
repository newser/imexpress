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

#ifndef __IEXP_INTEGRAL_QNG__
#define __IEXP_INTEGRAL_QNG__

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
class qng_t
{
  public:
    qng_t(T epsabs, T epsrel)
        : m_epsabs(epsabs)
        , m_epsrel(epsrel)
    {
    }

    T operator()(const typename unary_func<T>::type &fn,
                 T a,
                 T b,
                 T *abserr = nullptr,
                 size_t *neval = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T epsabs() const
    {
        return m_epsabs;
    }

    qng_t &epsabs(T e)
    {
        m_epsabs = e;
        return *this;
    }

    T epsrel() const
    {
        return m_epsrel;
    }

    qng_t &epsrel(T e)
    {
        m_epsrel = e;
        return *this;
    }

  private:
    qng_t(const qng_t &) = delete;
    qng_t &operator=(const qng_t &) = delete;

    T m_epsabs, m_epsrel;
};

template <>
double qng_t<double>::operator()(const typename unary_func<double>::type &fn,
                                 double a,
                                 double b,
                                 double *abserr,
                                 size_t *neval)
{
    unary_func<double> m_fn(fn);
    double r, e;
    size_t n;
    gsl_integration_qng(m_fn.gsl(),
                        a,
                        b,
                        m_epsabs,
                        m_epsrel,
                        &r,
                        abserr != nullptr ? abserr : &e,
                        neval != nullptr ? neval : &n);
    return r;
}

typedef qng_t<double> qng;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_QNG__ */
