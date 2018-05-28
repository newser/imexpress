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

#ifndef __IEXP_DAE_ODE__
#define __IEXP_DAE_ODE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <dae/def.h>
#include <dae/error.h>
#include <dae/function.h>
#include <dae/sunvec.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename Derived>
class ode
{
    DEFINE_DERIVED

  public:
    template <typename T>
    Derived &tolerance(double rel, const DenseBase<T> &abs)
    {
        sunvec_serial sv_abs(abs);
        cv_check(CVodeSVtolerances(m_cvode, rel, sv_abs.n_vector()));
        return derived();
    }

    Derived &tolerance(double rel, double abs)
    {
        cv_check(CVodeSStolerances(m_cvode, rel, abs));
        return derived();
    }

    template <typename T>
    int go(double &tout, DenseBase<T> &yout, bool one_step = false)
    {
        if (!m_ready) {
            derived().prepare();
            m_ready = true;
        }

        sunvec_serial y(yout, false);
        realtype tret;
        int e = CVode(m_cvode,
                      tout,
                      y.n_vector(),
                      &tret,
                      one_step ? CV_ONE_STEP : CV_NORMAL);
        cv_check(e);
        tout = tret;
        return e;
    }

    int compute_dy(double t, Map<const VectorXd> &y, Map<VectorXd> &dy)
    {
        return m_dy(t, y, dy);
    }

  protected:
    static void err_handler(int error_code,
                            const char *module,
                            const char *function,
                            char *msg,
                            void *user_data)
    {
        std::string es;
        es.reserve(256);
        es.append("(");
        es.append(std::to_string(error_code));
        es.append(")");
        es.append(msg);
        es.append(". see[");
        es.append(module);
        es.append(":");
        es.append(function);
        es.append("]");

        throw std::runtime_error(es);
    }

    template <typename T>
    ode(multistep lmm,
        iteration iter,
        const DY_TYPE &dy,
        double t0,
        const DenseBase<T> &y0)
        : m_cvode(nullptr)
        , m_dy(dy)
        , m_y0(y0)
        , m_ready(false)
    {
        static_assert(TYPE_IS(typename T::Scalar, realtype),
                      "only support realtype scalar");

        m_cvode = CVodeCreate((int)lmm, (int)iter);
        IEXP_NOT_NULLPTR(m_cvode);

        cv_check(CVodeInit(m_cvode,
                           dy_func<ode<Derived>>::s_dy,
                           t0,
                           m_y0.n_vector()));

        init();
    }

    ~ode()
    {
        CVodeFree(&m_cvode);
    }

    void init()
    {
        cv_check(CVodeSetUserData(m_cvode, this));

        cv_check(CVodeSetErrHandlerFn(m_cvode, err_handler, nullptr));
    }

    void *m_cvode;
    const DY_TYPE m_dy;
    sunvec_serial m_y0;
    bool m_ready : 1;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_ODE__ */
