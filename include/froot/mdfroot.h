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

#ifndef __IEXP_FROOT_MDFROOT__
#define __IEXP_FROOT_MDFROOT__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/gslvec.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// f_type definition
////////////////////////////////////////////////////////////

class mdfroot
{
    class fdf_func
    {
      public:
        using f_type = std::function<
            bool(Map<const VectorXd> &x, Map<VectorXd> &f, void *opaque)>;
        using df_type = std::function<
            bool(Map<const VectorXd> &x, Map<RowMatrixXd> &jac, void *opaque)>;
        using fdf_type = std::function<bool(Map<const VectorXd> &x,
                                            Map<VectorXd> &f,
                                            Map<RowMatrixXd> &jac,
                                            void *opaque)>;

        static int s_f(const gsl_vector *x, void *opaque, gsl_vector *f)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<VectorXd> mapped_f(f->data, f->size);
            return ((fdf_func *)opaque)->f(mapped_x, mapped_f) ? GSL_SUCCESS
                                                               : GSL_FAILURE;
        }

        static int s_df(const gsl_vector *x, void *opaque, gsl_matrix *jac)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<RowMatrixXd> mapped_jac(jac->data, jac->size1, jac->size2);
            return ((fdf_func *)opaque)->df(mapped_x, mapped_jac) ? GSL_SUCCESS
                                                                  : GSL_FAILURE;
        }

        static int s_fdf(const gsl_vector *x,
                         void *opaque,
                         gsl_vector *f,
                         gsl_matrix *jac)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<VectorXd> mapped_f(f->data, f->size);
            Map<RowMatrixXd> mapped_jac(jac->data, jac->size1, jac->size2);
            return ((fdf_func *)opaque)->fdf(mapped_x, mapped_f, mapped_jac)
                       ? GSL_SUCCESS
                       : GSL_FAILURE;
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

        const gsl_multiroot_function_fdf *gsl() const
        {
            return &m_gsl_mf;
        }

        bool f(Map<const VectorXd> &x, Map<VectorXd> &f) const
        {
            return m_f(x, f, m_opaque);
        }

        bool df(Map<const VectorXd> &x, Map<RowMatrixXd> &jac) const
        {
            return m_df(x, jac, m_opaque);
        }

        bool fdf(Map<const VectorXd> &x,
                 Map<VectorXd> &f,
                 Map<RowMatrixXd> &jac) const
        {
            return m_fdf(x, f, jac, m_opaque);
        }

      private:
        const f_type m_f;
        const df_type m_df;
        const fdf_type m_fdf;
        const gsl_multiroot_function_fdf m_gsl_mf;
        void *m_opaque;
    };

  public:
    enum class type
    {
        HYBRIDSJ,
        HYBRIDJ,
        NEWTON,
        GNEWTON,
    };

    template <typename T>
    mdfroot(const fdf_func::f_type &f,
            const fdf_func::df_type &df,
            const fdf_func::fdf_type &fdf,
            const DenseBase<T> &x,
            void *opaque,
            type t = type::HYBRIDSJ)
        : m_fdf(f, df, fdf, x.size(), opaque)
        , m_fsolver(
              gsl_multiroot_fdfsolver_alloc(s_fsolver_map[(int)t], x.size()))
    {
        IEXP_NOT_NULLPTR(m_fsolver);

        gslvec<typename T::Scalar> gv(x, false);
        gsl_multiroot_fdfsolver_set(m_fsolver,
                                    const_cast<gsl_multiroot_function_fdf *>(
                                        m_fdf.gsl()),
                                    gv.gsl_vector());
    }

    ~mdfroot()
    {
        gsl_multiroot_fdfsolver_free(m_fsolver);
    }

    const char *name()
    {
        return gsl_multiroot_fdfsolver_name(m_fsolver);
    }

    Map<VectorXd> find(double epsabs, double epsrel, int max_iter = 100)
    {
        int i = 0;
        do {
            gsl_multiroot_fdfsolver_iterate(m_fsolver);
        } while ((gsl_multiroot_test_delta(m_fsolver->dx,
                                           m_fsolver->x,
                                           epsabs,
                                           epsrel) == GSL_CONTINUE) &&
                 (i++ < max_iter));

        return Map<VectorXd>(m_fsolver->x->data, m_fsolver->x->size);
    }

    Map<VectorXd> find(double epsabs, int max_iter = 100)
    {
        int i = 0;
        do {
            gsl_multiroot_fdfsolver_iterate(m_fsolver);
        } while ((gsl_multiroot_test_residual(m_fsolver->f, epsabs) ==
                  GSL_CONTINUE) &&
                 (i++ < max_iter));

        return Map<VectorXd>(m_fsolver->x->data, m_fsolver->x->size);
    }

  private:
    static const gsl_multiroot_fdfsolver_type *s_fsolver_map[];

    mdfroot(const mdfroot &) = delete;
    mdfroot &operator=(const mdfroot &other) = delete;

    const fdf_func m_fdf;
    gsl_multiroot_fdfsolver *m_fsolver;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FROOT_MDFROOT__ */
