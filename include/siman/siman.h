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

#ifndef __IEXP_INTEGRAL_SIMAN__
#define __IEXP_INTEGRAL_SIMAN__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_siman.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class siman_test;

template <typename C>
class siman
{
    friend struct state;
    friend class siman_test;

  public:
    using f_energy = std::function<double(const C &)>;
    using f_step = std::function<void(rand::rng &, C &, double)>;
    using f_metric = std::function<double(const C &, const C &)>;
    using f_print = std::function<void(const C &)>;

    using f_copy = std::function<void(const C &, C &)>;
    using f_copy_construct = std::function<C &(const C &)>;

    siman(rand::rng::type type = DEFAULT_RNG_TYPE, unsigned long seed = 0)
        : m_rng(type, seed)
    {
    }

    siman &energy(const f_energy &e)
    {
        m_energy = e;
        return *this;
    }

    siman &step(const f_step &s)
    {
        m_step = s;
        return *this;
    }

    siman &metric(const f_metric &m)
    {
        m_metric = m;
        return *this;
    }

    siman &print(const f_print &p)
    {
        m_print = p;
        return *this;
    }

    siman &step_size(double step_size)
    {
        m_step_size = step_size;
        return *this;
    }

    siman &cooling(double k, double t_initial, double mu_t, double t_min)
    {
        m_k = k;
        m_t_initial = t_initial;
        m_mu_t = mu_t;
        m_t_min = t_min;
        return *this;
    }

    siman &n_tries(int n_tries)
    {
        m_n_tries = n_tries;
        return *this;
    }

    siman &iters_fixed_T(int iters_fixed_T)
    {
        m_iters_fixed_T = iters_fixed_T;
        return *this;
    }

    void solve(C &init)
    {
        state st(*this, &init);
        gsl_siman_params_t params{m_n_tries,
                                  m_iters_fixed_T,
                                  m_step_size,
                                  m_k,
                                  m_t_initial,
                                  m_mu_t,
                                  m_t_min};
        gsl_siman_solve(m_rng.gsl(),
                        &st,
                        s_Efunc,
                        s_step,
                        s_metric,
                        m_print != nullptr ? s_print : nullptr,
                        s_copy,
                        s_copy_construct,
                        s_destroy,
                        0,
                        params);

        // need not destroy passed param
        st.m_cfg = nullptr;
    }

  private:
    rand::rng m_rng;
    f_energy m_energy{nullptr};
    f_step m_step{nullptr};
    f_metric m_metric{nullptr};
    f_print m_print{nullptr};
    f_copy m_copy{nullptr};
    f_copy_construct m_copy_construct{nullptr};
    double m_step_size{0};
    double m_k{0};
    double m_t_initial{0};
    double m_mu_t{0};
    double m_t_min{0};
    int m_n_tries{0};
    int m_iters_fixed_T{0};

  private:
    struct state
    {
        state(siman &s, C *cfg)
            : m_siman(s)
            , m_cfg(cfg)
        {
        }

        ~state()
        {
            delete m_cfg;
        }

        siman &m_siman;
        C *m_cfg;
    };

    static double s_Efunc(void *xp)
    {
        const state *st = (state *)xp;
        return st->m_siman.m_energy(*st->m_cfg);
    }

    static void s_step(const gsl_rng *r, void *xp, double step_size)
    {
        const state *st = (state *)xp;
        return st->m_siman.m_step(st->m_siman.m_rng, *st->m_cfg, step_size);
    }

    static double s_metric(void *xp, void *yp)
    {
        const state *x = (state *)xp;
        const state *y = (state *)yp;
        return x->m_siman.m_metric(*x->m_cfg, *y->m_cfg);
    }

    static void s_print(void *xp)
    {
        const state *st = (state *)xp;
        st->m_siman.m_print(*st->m_cfg);
    }

    static void s_copy(void *source, void *dest)
    {
        const state *s = (state *)source;
        const state *d = (state *)dest;
        *d->m_cfg = *s->m_cfg;
    }

    static void *s_copy_construct(void *xp)
    {
        const state *st = (state *)xp;
        return new state(st->m_siman, new C(*st->m_cfg));
    }

    static void s_destroy(void *xp)
    {
        state *st = (state *)xp;
        delete st;
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_SIMAN__ */
