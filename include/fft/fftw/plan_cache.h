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

#ifndef __IEXP_FFT_FFTW_PLAN_CACHE__
#define __IEXP_FFT_FFTW_PLAN_CACHE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>
#include <fft/fftw/plan_double.h>
#include <fft/fftw/plan_float.h>

#include <map>
#include <mutex>

IEXP_NS_BEGIN

namespace fftw3 {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class plan_cache
{
  public:
    typedef plan<T> plan_t;
    typedef std::map<int64_t, plan_t> map1d_t;
    typedef std::tuple<int64_t, int64_t> key2d_t;
    typedef std::map<key2d_t, plan_t> map2d_t;

    enum how
    {
        MEASURE = FFTW_MEASURE,
        PATIENT = FFTW_PATIENT,
        EXHAUSTIVE = FFTW_EXHAUSTIVE,
    };

    // ========================================
    // 1d
    // ========================================

    void add(
        const int n, const void *i, const void *o, const bool fwd, const how h)
    {
        const int64_t k = key(n, i, o, fwd);

        std::lock_guard<std::mutex> g(m_lock1d);
        fftwl_plan p = m_plan_map1d[k];
        if (p == nullptr) {
            p = fftwl_plan_dft_1d(n,
                                  i,
                                  o,
                                  fwd ? FFTW_FORWARD : FFTW_BACKWARD,
                                  h);
            IEXP_NOT_NULLPTR(p);
            m_plan_map1d[k] = p;
        }
    }

    fftwl_plan find(const int n, const void *i, const void *o, const bool fwd)
    {
        const int64_t k = key(n, i, o, fwd);

        std::lock_guard<std::mutex> g(m_lock1d);
        return m_plan_map1d[k];
    }

    // ========================================
    // 2d
    // ========================================

    void add(const int n0,
             const int n1,
             const void *i,
             const void *o,
             const bool fwd,
             const how h)
    {
        const key2d_t k = key(n0, n1, i, o, fwd);

        std::lock_guard<std::mutex> g(m_lock2d);
        fftwl_plan p = m_plan_map2d[k];
        if (p == nullptr) {
            p = fftwl_plan_dft_2d(n0,
                                  n1,
                                  i,
                                  o,
                                  fwd ? FFTW_FORWARD : FFTW_BACKWARD,
                                  h);
            IEXP_NOT_NULLPTR(p);
            m_plan_map2d[k] = p;
        }
    }

    fftwl_plan find(const int n0,
                    const int n1,
                    const void *i,
                    const void *o,
                    const bool fwd)
    {
        const key2d_t k = key(n0, n1, i, o, fwd);

        std::lock_guard<std::mutex> g(m_lock2d);
        return m_plan_map2d[k];
    }

  private:
    int64_t key(const int n, const void *i, const void *o, const bool fwd)
    {
        const bool inplace(i == o);
        // avx512 require 512bit alignment
        const int i_align = reinterpret_cast<uintptr_t>(i) & 63;
        const int o_align = reinterpret_cast<uintptr_t>(o) & 63;

        return int64_t(fwd | (inplace << 1) | (i_align << 2) | (o_align << 8) |
                       ((int64_t)n << 32));
    }

    key2d_t key(const int n0,
                const int n1,
                const void *i,
                const void *o,
                const bool fwd)
    {
        const bool inplace(i == o);
        // avx512 require 512bit alignment
        const int i_align = reinterpret_cast<uintptr_t>(i) & 63;
        const int o_align = reinterpret_cast<uintptr_t>(o) & 63;

        return key2d_t(int64_t(fwd | (inplace << 1) | (i_align << 2) |
                               (o_align << 8) | ((int64_t)n0 << 32)),
                       int64_t(n1));
    }

    map1d_t m_plan_map1d;
    std::mutex m_lock1d;

    map2d_t m_plan_map2d;
    std::mutex m_lock2d;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_FFT_FFTW_PLAN_CACHE__ */
