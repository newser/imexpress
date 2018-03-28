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

#ifndef __IEXP_INTEGRAL_VEGAS__
#define __IEXP_INTEGRAL_VEGAS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <integral/function.h>
#include <integral/type.h>
#include <rand/rng.h>

#include <gsl/gsl_monte_vegas.h>

IEXP_NS_BEGIN

namespace integral {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

enum vegas_mode
{
    IMPORTANCE = GSL_VEGAS_MODE_IMPORTANCE,
    IMPORTANCE_ONLY = GSL_VEGAS_MODE_IMPORTANCE_ONLY,
    STRATIFIED = GSL_VEGAS_MODE_STRATIFIED,
};

template <typename T>
class vegas_t
{
  public:
    vegas_t(size_t dim,
            rand::rng_type type = rand::MT19937,
            unsigned long seed = 0)
        : m_state(nullptr)
        , m_rng(type, seed)
    {
        m_state = gsl_monte_vegas_alloc(dim);
        IEXP_NOT_NULLPTR(m_state);
    }

    ~vegas_t()
    {
        gsl_monte_vegas_free(m_state);
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

    double chisq()
    {
        return gsl_monte_vegas_chisq(m_state);
    }

    double alpha()
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        return p.alpha;
    }

    vegas_t &alpha(double v)
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        p.alpha = v;
        gsl_monte_vegas_params_set(m_state, &p);

        return *this;
    }

    size_t iterations()
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        return p.iterations;
    }

    vegas_t &iterations(size_t iterations)
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        p.iterations = iterations;
        gsl_monte_vegas_params_set(m_state, &p);

        return *this;
    }

    int stage()
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        return p.stage;
    }

    vegas_t &stage(int stage)
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        p.stage = stage;
        gsl_monte_vegas_params_set(m_state, &p);

        return *this;
    }

    vegas_mode mode()
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        return (vegas_mode)p.mode;
    }

    vegas_t &mode(vegas_mode mode)
    {
        gsl_monte_vegas_params p;
        gsl_monte_vegas_params_get(m_state, &p);
        p.mode = mode;
        gsl_monte_vegas_params_set(m_state, &p);

        return *this;
    }

  private:
    vegas_t(const vegas_t &) = delete;
    vegas_t &operator=(const vegas_t &) = delete;

    gsl_monte_vegas_state *m_state;
    rand::rng m_rng;
};

template <>
double vegas_t<double>::operator()(
    const typename monte_func<double>::type_1d &fn,
    double a,
    double b,
    size_t calls,
    double *abserr)
{
    gsl_monte_vegas_init(m_state);

    monte_func<double> m_fn(fn);
    double result, __abserr;
    gsl_monte_vegas_integrate(const_cast<gsl_monte_function *>(m_fn.gsl()),
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
double vegas_t<double>::operator()(
    const typename monte_func<double>::type_nd &fn,
    const double a[],
    const double b[],
    size_t dim,
    size_t calls,
    double *abserr)
{
    eigen_assert(dim == m_state->dim);

    gsl_monte_vegas_init(m_state);

    monte_func<double> m_fn(dim, fn);
    double result, __abserr;
    gsl_monte_vegas_integrate(const_cast<gsl_monte_function *>(m_fn.gsl()),
                              const_cast<double *>(a),
                              const_cast<double *>(b),
                              dim,
                              calls,
                              m_rng.gsl(),
                              m_state,
                              &result,
                              abserr != nullptr ? abserr : &__abserr);
    return result;
}

typedef vegas_t<double> vegas;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_VEGAS__ */
