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

#ifndef __IEXP_FMIN_MDFMIN__
#define __IEXP_FMIN_MDFMIN__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/gslvec.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// f_type definition
////////////////////////////////////////////////////////////

class mdfmin
{
    class fdf_func
    {
      public:
        using f_type =
            std::function<double(Map<const VectorXd> &x, void *opaque)>;
        using df_type = std::function<void(Map<const VectorXd> &x,
                                           Map<VectorXd> &gradient,
                                           void *opaque)>;
        using fdf_type = std::function<double(Map<const VectorXd> &x,
                                              Map<VectorXd> &gradient,
                                              void *opaque)>;

        static double s_f(const gsl_vector *x, void *opaque)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            return ((fdf_func *)opaque)->f(mapped_x);
        }

        static void s_df(const gsl_vector *x, void *opaque, gsl_vector *g)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<VectorXd> mapped_g(g->data, g->size);
            ((fdf_func *)opaque)->df(mapped_x, mapped_g);
        }

        static void s_fdf(const gsl_vector *x,
                          void *opaque,
                          double *f,
                          gsl_vector *g)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<VectorXd> mapped_g(g->data, g->size);
            *f = ((fdf_func *)opaque)->fdf(mapped_x, mapped_g);
        }

        fdf_func(const f_type &f,
                 const df_type &df,
                 const fdf_type &fdf,
                 size_t n,
                 void *opaque)
            : m_f(f)
            , m_df(df)
            , m_fdf(fdf)
            , m_gsl_mf{s_f, s_df, s_fdf, n, this}
            , m_opaque(opaque)
        {
        }

        const gsl_multimin_function_fdf *gsl() const
        {
            return &m_gsl_mf;
        }

        double f(Map<const VectorXd> &x) const
        {
            return m_f(x, m_opaque);
        }

        void df(Map<const VectorXd> &x, Map<VectorXd> &gradient) const
        {
            return m_df(x, gradient, m_opaque);
        }

        double fdf(Map<const VectorXd> &x, Map<VectorXd> &gradient) const
        {
            return m_fdf(x, gradient, m_opaque);
        }

      private:
        const f_type m_f;
        const df_type m_df;
        const fdf_type m_fdf;
        const gsl_multimin_function_fdf m_gsl_mf;
        void *m_opaque;
    };

  public:
    enum class type
    {
        CONJUGATE_FR,
        CONJUGATE_PR,
        VECTOR_BFGS,
        VECTOR_BFGS2,
        STEEPEST_DESCENT,
    };

    template <typename T>
    mdfmin(const fdf_func::f_type &f,
           const fdf_func::df_type &df,
           const fdf_func::fdf_type &fdf,
           const DenseBase<T> &x,
           double step,
           double tolerance,
           void *opaque,
           type t = type::CONJUGATE_FR)
        : m_fdf(f, df, fdf, x.size(), opaque)
        , m_fmin(gsl_multimin_fdfminimizer_alloc(s_fmin_map[(int)t], x.size()))
    {
        IEXP_NOT_NULLPTR(m_fmin);

        gslvec<typename T::Scalar> gv(x, false);
        gsl_multimin_fdfminimizer_set(m_fmin,
                                      const_cast<gsl_multimin_function_fdf *>(
                                          m_fdf.gsl()),
                                      gv.gsl_vector(),
                                      step,
                                      tolerance);
    }

    ~mdfmin()
    {
        gsl_multimin_fdfminimizer_free(m_fmin);
    }

    const char *name()
    {
        return gsl_multimin_fdfminimizer_name(m_fmin);
    }

    Map<VectorXd> find(double epsabs, int max_iter = 100)
    {
        int i = 0;
        do {
            gsl_multimin_fdfminimizer_iterate(m_fmin);
        } while ((gsl_multimin_test_gradient(gsl_multimin_fdfminimizer_gradient(
                                                 m_fmin),
                                             epsabs) == GSL_CONTINUE) &&
                 (i++ < max_iter));

        return Map<VectorXd>(m_fmin->x->data, m_fmin->x->size);
    }

  private:
    static const gsl_multimin_fdfminimizer_type *s_fmin_map[];

    mdfmin(const mdfmin &) = delete;
    mdfmin &operator=(const mdfmin &other) = delete;

    const fdf_func m_fdf;
    gsl_multimin_fdfminimizer *m_fmin;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FMIN_MDFMIN__ */
