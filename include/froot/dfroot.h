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

#ifndef __IEXP_FROOT_DFROOT__
#define __IEXP_FROOT_DFROOT__

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

class dfroot
{
    class fdf_func
    {
      public:
        using f_type = std::function<double(double)>;
        using fdf_type = std::function<void(double, double &, double &)>;

        static double s_f(double x, void *param)
        {
            return ((fdf_func *)param)->f(x);
        }

        static double s_df(double x, void *param)
        {
            return ((fdf_func *)param)->df(x);
        }

        static void s_fdf(double x, void *param, double *f, double *df)
        {
            return ((fdf_func *)param)->fdf(x, *f, *df);
        }

        fdf_func(const f_type &f,
                 const f_type &df,
                 const fdf_type &fdf,
                 void *param)
            : m_f(f)
            , m_df(df)
            , m_fdf(fdf)
            , m_gsl_fdf{s_f, s_df, s_fdf, param}
        {
        }

        const gsl_function_fdf *gsl() const
        {
            return &m_gsl_fdf;
        }

        double f(double x)
        {
            return m_f(x);
        }

        double df(double x)
        {
            return m_df(x);
        }

        void fdf(double x, double &f, double &df)
        {
            return m_fdf(x, f, df);
        }

      private:
        const f_type m_f, m_df;
        const fdf_type m_fdf;
        const gsl_function_fdf m_gsl_fdf;
    };

  public:
    enum class type
    {
        NEWTON,
        SECANT,
        STEFFENSON,
    };

    dfroot(const fdf_func::f_type &f,
           const fdf_func::f_type &df,
           const fdf_func::fdf_type &fdf,
           double x0,
           type t = type::STEFFENSON)
        : m_fdf(f, df, fdf, this)
        , m_fsolver(gsl_root_fdfsolver_alloc(s_fsolver_map[(int)t]))
    {
        IEXP_NOT_NULLPTR(m_fsolver);

        gsl_root_fdfsolver_set(m_fsolver,
                               const_cast<gsl_function_fdf *>(m_fdf.gsl()),
                               x0);
    }

    ~dfroot()
    {
        gsl_root_fdfsolver_free(m_fsolver);
    }

    const char *name()
    {
        return gsl_root_fdfsolver_name(m_fsolver);
    }

    double find(double epsabs, double epsrel, int max_iter = 100)
    {
        double x0, x1;
        x0 = gsl_root_fdfsolver_root(m_fsolver);

        int status, i = 0;
        do {
            gsl_root_fdfsolver_iterate(m_fsolver);
            x1 = gsl_root_fdfsolver_root(m_fsolver);
            status = gsl_root_test_delta(x1, x0, epsabs, epsrel);
            x0 = x1;
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return x1;
    }

    double find(double epsabs, int max_iter = 100)
    {
        double x;
        int status, i = 0;
        do {
            gsl_root_fdfsolver_iterate(m_fsolver);
            x = gsl_root_fdfsolver_root(m_fsolver);
            status = gsl_root_test_residual(m_fdf.f(x), epsabs);
        } while ((status == GSL_CONTINUE) && (i++ < max_iter));
        return x;
    }

  private:
    static const gsl_root_fdfsolver_type *s_fsolver_map[];

    dfroot(const dfroot &) = delete;
    dfroot &operator=(const dfroot &other) = delete;

    fdf_func m_fdf;
    gsl_root_fdfsolver *m_fsolver;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FROOT_DFROOT__ */
