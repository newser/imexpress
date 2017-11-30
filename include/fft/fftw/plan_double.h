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

#ifndef __IEXP_FFT_FFTW_PLAN_DOUBLE__
#define __IEXP_FFT_FFTW_PLAN_DOUBLE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <fft/fftw/plan.h>

IEXP_NS_BEGIN

namespace fftw3 {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <>
class plan<double>
{
  public:
    typedef double scalar_t;
    typedef std::complex<double> complex_t;

    plan()
        : m_plan(nullptr)
    {
    }

    ~plan()
    {
        if (m_plan != nullptr) {
            fftw_destroy_plan(m_plan);
        }
    }

    // ========================================
    // c2c
    // ========================================

    void fwd(const int n, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_1d(n,
                                      (fftw_complex *)i,
                                      (fftw_complex *)o,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft(m_plan, (fftw_complex *)i, (fftw_complex *)o);
    }

    void inv(const int n, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_1d(n,
                                      (fftw_complex *)i,
                                      (fftw_complex *)o,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft(m_plan, (fftw_complex *)i, (fftw_complex *)o);
    }

    void fwd(const int n0, const int n1, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_2d(n0,
                                      n1,
                                      (fftw_complex *)i,
                                      (fftw_complex *)o,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft(m_plan, (fftw_complex *)i, (fftw_complex *)o);
    }

    void inv(const int n0, const int n1, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_2d(n0,
                                      n1,
                                      (fftw_complex *)i,
                                      (fftw_complex *)o,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft(m_plan, (fftw_complex *)i, (fftw_complex *)o);
    }

    // ========================================
    // r2c
    // ========================================

    void fwd(const int n, const scalar_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_r2c_1d(n,
                                          (scalar_t *)i,
                                          (fftw_complex *)o,
                                          FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft_r2c(m_plan, (scalar_t *)i, (fftw_complex *)o);
    }

    void inv(const int n, const complex_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_c2r_1d(n,
                                          (fftw_complex *)i,
                                          o,
                                          FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft_c2r(m_plan, (fftw_complex *)i, o);
    }

    void fwd(const int n0, const int n1, const scalar_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_r2c_2d(n0,
                                          n1,
                                          (scalar_t *)i,
                                          (fftw_complex *)o,
                                          FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft_r2c(m_plan, (scalar_t *)i, (fftw_complex *)o);
    }

    void inv(const int n0, const int n1, const complex_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_dft_c2r_2d(n0,
                                          n1,
                                          (fftw_complex *)i,
                                          o,
                                          FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_dft_c2r(m_plan, (fftw_complex *)i, o);
    }

    // ========================================
    // r2r
    // ========================================

    template <kind k>
    void fwd(const int n, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_r2r_1d(n,
                                      (scalar_t *)i,
                                      o,
                                      fwd_kind[k],
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <kind k>
    void inv(const int n, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_r2r_1d(n,
                                      (scalar_t *)i,
                                      o,
                                      inv_kind[k],
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <kind k0, kind k1>
    void fwd(const int n0, const int n1, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_r2r_2d(n0,
                                      n1,
                                      (scalar_t *)i,
                                      o,
                                      fwd_kind[k0],
                                      fwd_kind[k1],
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <kind k0, kind k1>
    void inv(const int n0, const int n1, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            m_plan = fftw_plan_r2r_2d(n0,
                                      n1,
                                      (scalar_t *)i,
                                      o,
                                      inv_kind[k0],
                                      inv_kind[k1],
                                      FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
        }

        fftw_execute_r2r(m_plan, (scalar_t *)i, o);
    }

  private:
    fftw_plan m_plan;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_FFT_FFTW_PLAN_DOUBLE__ */
