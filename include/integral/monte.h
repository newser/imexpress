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

#ifndef __IEXP_INTEGRAL_MONTE__
#define __IEXP_INTEGRAL_MONTE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>
#include <integral/type.h>
#include <rand/rng.h>

#include <gsl/gsl_monte_plain.h>

IEXP_NS_BEGIN

namespace integral {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class monte_t
{
  public:
    monte_t(size_t dim,
            rand::rng_type type = rand::MT19937,
            unsigned long seed = 0)
        : m_state(nullptr)
        , m_rng(type, seed)
    {
        m_state = gsl_monte_plain_alloc(dim);
        IEXP_NOT_NULLPTR(m_state);
    }

    monte_t(rand::rng_type type = rand::MT19937, unsigned long seed = 0)
        : m_state(nullptr)
        , m_rng(type, seed)
    {
        m_state = gsl_monte_plain_alloc(1);
        IEXP_NOT_NULLPTR(m_state);
    }

    ~monte_t()
    {
        gsl_monte_plain_free(m_state);
    }

    T operator()(const typename monte_func<T>::type_1d &fn,
                 const T a,
                 const T b,
                 size_t calls,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T operator()(const typename monte_func<T>::type_nd &fn,
                 const T a[],
                 const T b[],
                 size_t calls,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

  private:
    monte_t(const monte_t &) = delete;
    monte_t &operator=(const monte_t &) = delete;

    gsl_monte_plain_state *m_state;
    rand::rng m_rng;
};


template <>
double monte_t<double>::operator()(
    const typename monte_func<double>::type_1d &fn,
    const double a,
    const double b,
    size_t calls,
    double *abserr)
{
    gsl_monte_plain_init(m_state);

    monte_func<double> m_fn(fn);
    double result, __abserr;
    gsl_monte_plain_integrate(m_fn.gsl(),
                              &a,
                              &b,
                              1,
                              calls,
                              m_rng.gsl(),
                              m_state,
                              &result,
                              abserr != nullptr ? abserr : &__abserr);
    return result;
}

template <>
double monte_t<double>::operator()(
    const typename monte_func<double>::type_nd &fn,
    const double a[],
    const double b[],
    size_t calls,
    double *abserr)
{
    gsl_monte_plain_init(m_state);

    monte_func<double> m_fn(m_state->dim, fn);
    double result, __abserr;
    gsl_monte_plain_integrate(m_fn.gsl(),
                              a,
                              b,
                              m_state->dim,
                              calls,
                              m_rng.gsl(),
                              m_state,
                              &result,
                              abserr != nullptr ? abserr : &__abserr);
    return result;
}

typedef monte_t<double> monte;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_MONTE__ */
