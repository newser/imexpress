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

#ifndef __IEXP_FMIN_FMIN__
#define __IEXP_FMIN_FMIN__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class fmin
{
  public:
    enum class type
    {
        GOLDENSECTION,
        BRENT,
        QUAD_GOLDEN,
    };

    fmin(const unary_func<double>::type &fn,
         double lower,
         double upper,
         double guess,
         type t = type::BRENT)
        : m_fn(fn)
        , m_fmin(gsl_min_fminimizer_alloc(s_fmin_map[(int)t]))
    {
        IEXP_NOT_NULLPTR(m_fmin);

        eigen_assert(lower <= upper);
        gsl_min_fminimizer_set(m_fmin,
                               const_cast<gsl_function *>(m_fn.gsl()),
                               guess,
                               lower,
                               upper);
    }

    ~fmin()
    {
        gsl_min_fminimizer_free(m_fmin);
    }

    const char *name()
    {
        return gsl_min_fminimizer_name(m_fmin);
    }

    double find(double epsabs, double epsrel, int max_iter = 100)
    {
        int status, i = 0;
        do {
            gsl_min_fminimizer_iterate(m_fmin);
            status = gsl_min_test_interval(gsl_min_fminimizer_x_lower(m_fmin),
                                           gsl_min_fminimizer_x_upper(m_fmin),
                                           epsabs,
                                           epsrel);
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return gsl_min_fminimizer_x_minimum(m_fmin);
    }

  private:
    static const gsl_min_fminimizer_type *s_fmin_map[];

    fmin(const fmin &) = delete;
    fmin &operator=(const fmin &other) = delete;

    const unary_func<double> m_fn;
    gsl_min_fminimizer *m_fmin;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FMIN_FMIN__ */
