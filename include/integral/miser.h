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

#ifndef __IEXP_INTEGRAL_MISER__
#define __IEXP_INTEGRAL_MISER__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>
#include <integral/type.h>
#include <rand/rng.h>

#include <gsl/gsl_monte_miser.h>

IEXP_NS_BEGIN

namespace integral {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class miser_t
{
  public:
    miser_t(size_t dim,
            rand::rng_type type = rand::MT19937,
            unsigned long seed = 0)
        : m_state(nullptr)
        , m_rng(type, seed)
    {
        m_state = gsl_monte_miser_alloc(dim);
        IEXP_NOT_NULLPTR(m_state);
    }

    ~miser_t()
    {
        gsl_monte_miser_free(m_state);
    }

    T operator()(const typename monte_func<T>::type_1d &fn,
                 T a,
                 T b,
                 size_t calls,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    T operator()(const typename monte_func<T>::type_nd &fn,
                 const T a[],
                 const T b[],
                 size_t dim,
                 size_t calls,
                 T *abserr = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    double estimate_frac()
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        return p.estimate_frac;
    }

    miser_t &estimate_frac(double frac)
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        p.estimate_frac = frac;
        gsl_monte_miser_params_set(m_state, &p);

        return *this;
    }

    size_t min_calls()
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        return p.min_calls;
    }

    miser_t &min_calls(size_t calls)
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        p.min_calls = calls;
        gsl_monte_miser_params_set(m_state, &p);

        return *this;
    }

    size_t min_calls_per_bisection()
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        return p.min_calls_per_bisection;
    }

    miser_t &min_calls_per_bisection(size_t calls)
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        p.min_calls_per_bisection = calls;
        gsl_monte_miser_params_set(m_state, &p);

        return *this;
    }

    double alpha()
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        return p.alpha;
    }

    miser_t &alpha(double alpha)
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        p.alpha = alpha;
        gsl_monte_miser_params_set(m_state, &p);

        return *this;
    }

    double dither()
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        return p.dither;
    }

    miser_t &dither(double dither)
    {
        gsl_monte_miser_params p;
        gsl_monte_miser_params_get(m_state, &p);
        p.dither = dither;
        gsl_monte_miser_params_set(m_state, &p);

        return *this;
    }

  private:
    miser_t(const miser_t &) = delete;
    miser_t &operator=(const miser_t &) = delete;

    gsl_monte_miser_state *m_state;
    rand::rng m_rng;
};

template <>
double miser_t<double>::operator()(
    const typename monte_func<double>::type_1d &fn,
    double a,
    double b,
    size_t calls,
    double *abserr)
{
    gsl_monte_miser_init(m_state);

    monte_func<double> m_fn(fn);
    double result, __abserr;
    gsl_monte_miser_integrate(const_cast<gsl_monte_function *>(m_fn.gsl()),
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
double miser_t<double>::operator()(
    const typename monte_func<double>::type_nd &fn,
    const double a[],
    const double b[],
    size_t dim,
    size_t calls,
    double *abserr)
{
    eigen_assert(dim == m_state->dim);

    gsl_monte_miser_init(m_state);

    monte_func<double> m_fn(dim, fn);
    double result, __abserr;
    gsl_monte_miser_integrate(const_cast<gsl_monte_function *>(m_fn.gsl()),
                              a,
                              b,
                              dim,
                              calls,
                              m_rng.gsl(),
                              m_state,
                              &result,
                              abserr != nullptr ? abserr : &__abserr);
    return result;
}

typedef miser_t<double> miser;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_MISER__ */
