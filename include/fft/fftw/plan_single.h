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

#ifndef __IEXP_FFT_FFTW_PLAN_SINGLE__
#define __IEXP_FFT_FFTW_PLAN_SINGLE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <fft/fftw/plan.h>

IEXP_NS_BEGIN

namespace fftw3 {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define PLAN_LOCK                                                              \
    if (m_lock != nullptr) {                                                   \
        m_lock->lock();                                                        \
    }                                                                          \
    /* check again, as another thread may already create it */                 \
    if (m_plan == nullptr) {
#define PLAN_UNLOCK                                                            \
    }                                                                          \
    if (m_lock != nullptr) {                                                   \
        m_lock->unlock();                                                      \
    }

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <>
class plan<float>
{
  public:
    using scalar_t = float;
    using complex_t = std::complex<float>;

    plan(std::mutex *lock = nullptr)
        : m_plan(nullptr)
        , m_lock(lock)
    {
    }

    ~plan()
    {
        if (m_plan != nullptr) {
            fftwf_destroy_plan(m_plan);
        }
    }

    // ========================================
    // c2c
    // ========================================

    void fwd(int n, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_1d(n,
                                       (fftwf_complex *)i,
                                       (fftwf_complex *)o,
                                       FFTW_FORWARD,
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft(m_plan, (fftwf_complex *)i, (fftwf_complex *)o);
    }

    void inv(int n, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_1d(n,
                                       (fftwf_complex *)i,
                                       (fftwf_complex *)o,
                                       FFTW_BACKWARD,
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft(m_plan, (fftwf_complex *)i, (fftwf_complex *)o);
    }

    void fwd(int n0, int n1, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_2d(n0,
                                       n1,
                                       (fftwf_complex *)i,
                                       (fftwf_complex *)o,
                                       FFTW_FORWARD,
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft(m_plan, (fftwf_complex *)i, (fftwf_complex *)o);
    }

    void inv(int n0, int n1, const complex_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_2d(n0,
                                       n1,
                                       (fftwf_complex *)i,
                                       (fftwf_complex *)o,
                                       FFTW_BACKWARD,
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft(m_plan, (fftwf_complex *)i, (fftwf_complex *)o);
    }

    // ========================================
    // r2c
    // ========================================

    // num of i: n
    // num of o: (n/2 + 1)
    void fwd(int n, const scalar_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_r2c_1d(n,
                                           (scalar_t *)i,
                                           (fftwf_complex *)o,
                                           FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft_r2c(m_plan, (scalar_t *)i, (fftwf_complex *)o);
    }

    // num of i: (n/2 + 1)
    // num of o: n
    void inv(int n, const complex_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_c2r_1d(n,
                                           (fftwf_complex *)i,
                                           o,
                                           FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft_c2r(m_plan, (fftwf_complex *)i, o);
    }

    // num of i: n0 * n1
    // num of o: n0 * (n1/2 + 1)
    void fwd(int n0, int n1, const scalar_t *i, complex_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_dft_r2c_2d(n0,
                                           n1,
                                           (scalar_t *)i,
                                           (fftwf_complex *)o,
                                           FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_dft_r2c(m_plan, (scalar_t *)i, (fftwf_complex *)o);
    }

    // num of i: n0 * (n1/2 + 1)
    // num of o: n0 * n1
    void inv(int n0, int n1, const complex_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            // can not use FFTW_PRESERVE_INPUT for c2r
            fftwf_complex *__i = new fftwf_complex[n0 * n1];
            IEXP_NOT_NULLPTR(__i);

            m_plan = fftwf_plan_dft_c2r_2d(n0, n1, __i, o, FFTW_ESTIMATE);
            IEXP_NOT_NULLPTR(m_plan);

            delete[] __i;
            PLAN_UNLOCK
        }

        // Note i would be modified!!!
        fftwf_execute_dft_c2r(m_plan, (fftwf_complex *)i, o);
    }

    // ========================================
    // r2r
    // ========================================

    template <fft::kind k>
    void fwd(int n, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_r2r_1d(n,
                                       (scalar_t *)i,
                                       o,
                                       fwd_kind[k],
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <fft::kind k>
    void inv(int n, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_r2r_1d(n,
                                       (scalar_t *)i,
                                       o,
                                       inv_kind[k],
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <fft::kind k0, fft::kind k1>
    void fwd(int n0, int n1, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_r2r_2d(n0,
                                       n1,
                                       (scalar_t *)i,
                                       o,
                                       fwd_kind[k0],
                                       fwd_kind[k1],
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_r2r(m_plan, (scalar_t *)i, o);
    }

    template <fft::kind k0, fft::kind k1>
    void inv(int n0, int n1, const scalar_t *i, scalar_t *o)
    {
        if (m_plan == nullptr) {
            PLAN_LOCK
            m_plan = fftwf_plan_r2r_2d(n0,
                                       n1,
                                       (scalar_t *)i,
                                       o,
                                       inv_kind[k0],
                                       inv_kind[k1],
                                       FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            IEXP_NOT_NULLPTR(m_plan);
            PLAN_UNLOCK
        }

        fftwf_execute_r2r(m_plan, (scalar_t *)i, o);
    }

  private:
    fftwf_plan m_plan;
    std::mutex *m_lock;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#undef PLAN_LOCK
#undef PLAN_UNLOCK

#endif /* __IEXP_FFT_FFTW_PLAN_SINGLE__ */
