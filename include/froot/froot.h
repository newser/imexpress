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

#ifndef __IEXP_FROOT_FROOT__
#define __IEXP_FROOT_FROOT__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class froot
{
  public:
    enum class type
    {
        BISECTION,
        FALSEPOS,
        BRENT,
    };

    froot(const unary_func<double>::type &fn,
          double lower,
          double upper,
          type t = type::BRENT)
        : m_fn(fn)
        , m_fsolver(gsl_root_fsolver_alloc(s_fsolver_map[(int)t]))
    {
        IEXP_NOT_NULLPTR(m_fsolver);

        eigen_assert(lower <= upper);
        gsl_root_fsolver_set(m_fsolver,
                             const_cast<gsl_function *>(m_fn.gsl()),
                             lower,
                             upper);
    }

    ~froot()
    {
        gsl_root_fsolver_free(m_fsolver);
    }

    const char *name()
    {
        return gsl_root_fsolver_name(m_fsolver);
    }

    template <bool test_interval = true>
    double find(double epsabs, double epsrel, int max_iter = 100)
    {
        return find(epsabs, epsrel, max_iter, TYPE_BOOL(test_interval)());
    }

    double find(double epsabs, int max_iter = 100)
    {
        double x;
        int status, i = 0;
        do {
            gsl_root_fsolver_iterate(m_fsolver);
            x = gsl_root_fsolver_root(m_fsolver);
            status = gsl_root_test_residual(m_fn(x), epsabs);
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return x;
    }

  private:
    static const gsl_root_fsolver_type *s_fsolver_map[];

    froot(const froot &) = delete;
    froot &operator=(const froot &other) = delete;

    double find(double epsabs, double epsrel, int max_iter, std::true_type)
    {
        int status, i = 0;
        do {
            gsl_root_fsolver_iterate(m_fsolver);
            status = gsl_root_test_interval(gsl_root_fsolver_x_lower(m_fsolver),
                                            gsl_root_fsolver_x_upper(m_fsolver),
                                            epsabs,
                                            epsrel);
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return gsl_root_fsolver_root(m_fsolver);
    }

    double find(double epsabs, double epsrel, int max_iter, std::false_type)
    {
        double x0, x1;
        gsl_root_fsolver_iterate(m_fsolver);
        x0 = gsl_root_fsolver_root(m_fsolver);

        int status, i = 0;
        do {
            gsl_root_fsolver_iterate(m_fsolver);
            x1 = gsl_root_fsolver_root(m_fsolver);
            status = gsl_root_test_delta(x1, x0, epsabs, epsrel);
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return x1;
    }

    unary_func<double> m_fn;
    gsl_root_fsolver *m_fsolver;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FROOT_FROOT__ */
