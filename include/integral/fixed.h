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

#ifndef __IEXP_INTEGRAL_FIXED__
#define __IEXP_INTEGRAL_FIXED__

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
class fixed_t
{
  public:
    enum class type
    {
        LEGENDRE,
        CHEBYSHEV,
        GEGENBAUER,
        JACOBI,
        LAGUERRE,
        HERMITE,
        EXPONENTIAL,
        RATIONAL,
        CHEBYSHEV2
    };

    fixed_t(type t, size_t n, T a, T b, T alpha, T beta)
        : m_type(t)
        , m_n(n)
        , m_a(a)
        , m_b(b)
        , m_alpha(alpha)
        , m_beta(beta)
        , m_workspace(nullptr)
    {
    }

    ~fixed_t()
    {
        if (m_workspace != nullptr) {
            gsl_integration_fixed_free(m_workspace);
        }
    }

    T operator()(const typename unary_func<T>::type &fn, void *opaque = nullptr)
    {
        UNSUPPORTED_TYPE(T);
    }

    size_t n()
    {
        init();
        return gsl_integration_fixed_n(m_workspace);
    }

    T a() const
    {
        return m_a;
    }

    T b() const
    {
        return m_a;
    }

    T alpha() const
    {
        return m_alpha;
    }

    T beta() const
    {
        return m_beta;
    }

    const T *x()
    {
        UNSUPPORTED_TYPE(T);
    }

    const T *w()
    {
        UNSUPPORTED_TYPE(T);
    }

  private:
    fixed_t(const fixed_t &) = delete;
    fixed_t &operator=(const fixed_t &) = delete;

    type m_type;
    size_t m_n;
    T m_a, m_b, m_alpha, m_beta;
    gsl_integration_fixed_workspace *m_workspace;

    static const gsl_integration_fixed_type *gsl_type(type t)
    {
        switch (t) {
            case type::LEGENDRE: {
                return gsl_integration_fixed_legendre;
            } break;
            case type::CHEBYSHEV: {
                return gsl_integration_fixed_chebyshev;
            } break;
            case type::GEGENBAUER: {
                return gsl_integration_fixed_gegenbauer;
            } break;
            case type::JACOBI: {
                return gsl_integration_fixed_jacobi;
            } break;
            case type::LAGUERRE: {
                return gsl_integration_fixed_laguerre;
            } break;
            case type::HERMITE: {
                return gsl_integration_fixed_hermite;
            } break;
            case type::EXPONENTIAL: {
                return gsl_integration_fixed_exponential;
            } break;
            case type::RATIONAL: {
                return gsl_integration_fixed_rational;
            } break;
            case type::CHEBYSHEV2: {
                return gsl_integration_fixed_chebyshev2;
            } break;
            default: {
                throw std::invalid_argument("invalid type");
            }
        }
    }

    void init()
    {
        if (m_workspace == nullptr) {
            m_workspace = gsl_integration_fixed_alloc(gsl_type(m_type),
                                                      m_n,
                                                      m_a,
                                                      m_b,
                                                      m_alpha,
                                                      m_beta);
            IEXP_NOT_NULLPTR(m_workspace);
        }
    }
};

template <>
double fixed_t<double>::operator()(const typename unary_func<double>::type &fn,
                                   void *opaque)
{
    init();

    unary_func<double> m_fn(fn, opaque);
    double r;
    gsl_integration_fixed(m_fn.gsl(), &r, m_workspace);
    return r;
}

template <>
const double *fixed_t<double>::x()
{
    init();
    return gsl_integration_fixed_nodes(m_workspace);
}

template <>
const double *fixed_t<double>::w()
{
    init();
    return gsl_integration_fixed_weights(m_workspace);
}

typedef fixed_t<double> fixed;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_FIXED__ */
