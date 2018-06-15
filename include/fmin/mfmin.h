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

#ifndef __IEXP_FMIN_MFMIN__
#define __IEXP_FMIN_MFMIN__

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
// type definition
////////////////////////////////////////////////////////////

class mfmin
{
    class f_func
    {
      public:
        using type = std::function<double(Map<const VectorXd> &x)>;

        static double s_f(const gsl_vector *x, void *param)
        {
            Map<const VectorXd> mapped_x(x->data, x->size);
            return ((f_func *)param)->f(mapped_x);
        }

        f_func(const type &f, size_t n, void *param)
            : m_f(f)
            , m_gsl_mf{s_f, n, param}
        {
        }

        const gsl_multimin_function *gsl() const
        {
            return &m_gsl_mf;
        }

        double f(Map<const VectorXd> &x)
        {
            return m_f(x);
        }

      private:
        const type m_f;
        const gsl_multimin_function m_gsl_mf;
    };

  public:
    enum class type
    {
        NMSIMPLEX,
        NMSIMPLEX2,
        NMSIMPLEX2RAND,
    };

    template <typename T, typename U>
    mfmin(const f_func::type &f,
          const DenseBase<T> &x,
          const DenseBase<U> &step,
          type t = type::NMSIMPLEX2RAND)
        : m_fn(f, x.size(), const_cast<f_func *>(&m_fn))
        , m_fmin(gsl_multimin_fminimizer_alloc(s_mfmin_map[(int)t], x.size()))
    {
        IEXP_NOT_NULLPTR(m_fmin);

        static_assert(TYPE_IS(typename T::Scalar, typename U::Scalar),
                      "must be same scalar");
        eigen_assert(x.size() == step.size());

        gslvec<typename T::Scalar> gx(x, false);
        gslvec<typename U::Scalar> gs(x, false);
        gsl_multimin_fminimizer_set(m_fmin,
                                    const_cast<gsl_multimin_function *>(
                                        m_fn.gsl()),
                                    gx.gsl_vector(),
                                    gs.gsl_vector());
    }

    ~mfmin()
    {
        gsl_multimin_fminimizer_free(m_fmin);
    }

    const char *name()
    {
        return gsl_multimin_fminimizer_name(m_fmin);
    }

    Map<VectorXd> find(double epsabs, int max_iter = 100)
    {
        int i = 0;
        do {
            gsl_multimin_fminimizer_iterate(m_fmin);
        } while ((gsl_multimin_test_size(gsl_multimin_fminimizer_size(m_fmin),
                                         epsabs) == GSL_CONTINUE) &&
                 (i++ < max_iter));

        return Map<VectorXd>(m_fmin->x->data, m_fmin->x->size);
    }

  private:
    static const gsl_multimin_fminimizer_type *s_mfmin_map[];

    mfmin(const mfmin &) = delete;
    mfmin &operator=(const mfmin &other) = delete;

    const f_func m_fn;
    gsl_multimin_fminimizer *m_fmin;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_FMIN_MFMIN__ */
