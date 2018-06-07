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

#ifndef __IEXP_FROOT_MFROOT__
#define __IEXP_FROOT_MFROOT__

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
// type definition
////////////////////////////////////////////////////////////

class mfroot
{
    class f_func
    {
      public:
        using type =
            std::function<bool(Map<const VectorXd> &x, Map<VectorXd> &f)>;

        static int s_f(const gsl_vector *x, void *param, gsl_vector *f)
        {
            eigen_assert(x->size == f->size);

            Map<const VectorXd> mapped_x(x->data, x->size);
            Map<VectorXd> mapped_f(f->data, f->size);
            return ((f_func *)param)->f(mapped_x, mapped_f) ? GSL_SUCCESS
                                                            : GSL_FAILURE;
        }

        f_func(const type &f, size_t n, void *param)
            : m_f(f)
            , m_gsl_mf{s_f, n, param}
        {
        }

        const gsl_multiroot_function *gsl() const
        {
            return &m_gsl_mf;
        }

        bool f(Map<const VectorXd> &x, Map<VectorXd> &f)
        {
            return m_f(x, f);
        }

      private:
        const type m_f;
        const gsl_multiroot_function m_gsl_mf;
    };

  public:
    enum class type
    {
        HYBRIDS,
        HYBRID,
        DNEWTON,
        BROYDEN,
    };

    template <typename T>
    mfroot(const f_func::type &f, const DenseBase<T> &x, type t = type::HYBRIDS)
        : m_fn(f, x.size(), const_cast<f_func *>(&m_fn))
        , m_fsolver(
              gsl_multiroot_fsolver_alloc(s_fsolver_map[(int)t], x.size()))
    {
        IEXP_NOT_NULLPTR(m_fsolver);

        gslvec<typename T::Scalar> gv(x, false);
        gsl_multiroot_fsolver_set(m_fsolver,
                                  const_cast<gsl_multiroot_function *>(
                                      m_fn.gsl()),
                                  gv.gsl_vector());
    }

    ~mfroot()
    {
        gsl_multiroot_fsolver_free(m_fsolver);
    }

    const char *name()
    {
        return gsl_multiroot_fsolver_name(m_fsolver);
    }

    Map<VectorXd> find(double epsabs, double epsrel, int max_iter = 100)
    {
        int i = 0;
        do {
            gsl_multiroot_fsolver_iterate(m_fsolver);
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
            gsl_multiroot_fsolver_iterate(m_fsolver);
        } while ((gsl_multiroot_test_residual(m_fsolver->f, epsabs) ==
                  GSL_CONTINUE) &&
                 (i++ < max_iter));

        return Map<VectorXd>(m_fsolver->x->data, m_fsolver->x->size);
    }

  private:
    static const gsl_multiroot_fsolver_type *s_fsolver_map[];

    mfroot(const mfroot &) = delete;
    mfroot &operator=(const mfroot &other) = delete;

    const f_func m_fn;
    gsl_multiroot_fsolver *m_fsolver;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FROOT_MFROOT__ */
