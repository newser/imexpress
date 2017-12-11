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
class qng_t
{
  public:
    qng_t(const typename unary_func<T>::type &fn, T epsabs, T epsrel)
        : m_fn(fn)
        , m_epsabs(epsabs)
        , m_epsrel(epsrel)
    {
    }

    int operator()(const T a,
                   const T b,
                   T *result,
                   T *abserr = nullptr,
                   size_t *neval = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

  private:
    unary_func<T> m_fn;
    const T m_epsabs, m_epsrel;
};

template <>
int qng_t<double>::operator()(const double a,
                              const double b,
                              double *result,
                              double *abserr,
                              size_t *neval)
{
    double __abserr;
    size_t __neval;
    return gsl_integration_qng(m_fn.gsl(),
                               a,
                               b,
                               m_epsabs,
                               m_epsrel,
                               result,
                               abserr != nullptr ? abserr : &__abserr,
                               neval != nullptr ? neval : &__neval);
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
