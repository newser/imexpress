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
#include <fft/fftw/plan.h>

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

enum how
{
    MEASURE = FFTW_MEASURE,
    PATIENT = FFTW_PATIENT,
    EXHAUSTIVE = FFTW_EXHAUSTIVE,
};

enum scalar
{
    SINGLE,
    DOUBLE,
};

enum io
{
    C2C,
    R2C,
    C2R,
    R2R,
};

template <typename T, int dim>
struct plan_traits
{
};

template <typename T>
struct plan_traits<T, 1>
{
    using plan_t = plan<T>;
    using key_t = int64_t;
    using map_t = std::map<key_t, plan_t>;
};

template <typename T>
struct plan_traits<T, 2>
{
    using plan_t = plan<T>;
    using key_t = std::tuple<int64_t, int64_t>;
    using map_t = std::map<key_t, plan_t>;
};

template <typename I, typename O>
struct io_traits
{
};

template <>
struct io_traits<std::complex<double>, std::complex<double>>
{
    using type = double;
    static const scalar scalar = DOUBLE;
    static const io io = C2C;
};

template <>
struct io_traits<double, std::complex<double>>
{
    using type = double;
    static const scalar scalar = DOUBLE;
    static const io io = R2C;
};

template <>
struct io_traits<std::complex<double>, double>
{
    using type = double;
    static const scalar scalar = DOUBLE;
    static const io io = C2R;
};

template <>
struct io_traits<double, double>
{
    using type = double;
    static const scalar scalar = DOUBLE;
    static const io io = R2R;
};

template <>
struct io_traits<std::complex<float>, std::complex<float>>
{
    using type = float;
    static const scalar scalar = SINGLE;
    static const io io = C2C;
};

template <>
struct io_traits<float, std::complex<float>>
{
    using type = float;
    static const scalar scalar = SINGLE;
    static const io io = R2C;
};

template <>
struct io_traits<std::complex<float>, float>
{
    using type = float;
    static const scalar scalar = SINGLE;
    static const io io = C2R;
};

template <>
struct io_traits<float, float>
{
    using type = float;
    static const scalar scalar = SINGLE;
    static const io io = R2R;
};

template <typename T, int dim, kind k0, kind k1>
class plan_cache
{
  public:
    using plan_t = typename plan_traits<T, dim>::plan_t;
    using key_t = typename plan_traits<T, dim>::key_t;
    using map_t = typename plan_traits<T, dim>::map_t;

    template <typename I, typename O>
    plan_t &get(const int n,
                const I *i,
                const O *o,
                const bool fwd,
                const how h = MEASURE)
    {
        static_assert(std::is_same<T, typename io_traits<I, O>::type>::value,
                      "invalid type");

        const int64_t kval = key(n, i, o, fwd);
        {
            std::lock_guard<std::mutex> g(m_lock);

            auto r = m_plan_map.insert(
                typename map_t::value_type(kval, plan_t(&m_lock)));
            return r.first->second;
        }
    }

    template <typename I, typename O>
    plan_t &get(const int n0,
                const int n1,
                const I *i,
                const O *o,
                const bool fwd,
                const how h = MEASURE)
    {
        static_assert(std::is_same<T, typename io_traits<I, O>::type>::value,
                      "invalid type");

        const key_t kval = key(n0, n1, i, o, fwd);
        {
            std::lock_guard<std::mutex> g(m_lock);

            auto r = m_plan_map.insert(
                typename map_t::value_type(kval, plan_t(&m_lock)));
            return r.first->second;
        }
    }

  private:
    template <typename I, typename O>
    key_t key(const int n, const I *i, const O *o, const bool fwd)
    {
        const bool inplace((uintptr_t)i == (uintptr_t)o);

        return key_t(fwd | (inplace << 1) | (io_traits<I, O>::scalar << 2) |
                     (io_traits<I, O>::io << 4) | (k0 << 8) |
                     ((int64_t)n << 32));
    }

    template <typename I, typename O>
    key_t key(
        const int n0, const int n1, const I *i, const O *o, const bool fwd)
    {
        const bool inplace((uintptr_t)i == (uintptr_t)o);

        return key_t(int64_t(fwd | (inplace << 1) |
                             (io_traits<I, O>::scalar << 2) |
                             (io_traits<I, O>::io << 4) | (k0 << 8) |
                             (k1 << 12) | ((int64_t)n0 << 32)),
                     int64_t(n1));
    }

    map_t m_plan_map;
    std::mutex m_lock;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

template <kind k = KIND_NUM, typename I = void, typename O = void>
plan<typename io_traits<I, O>::type> &get_plan(
    const int n, const I *i, const O *o, const bool fwd, const how h = MEASURE)
{
    static plan_cache<typename io_traits<I, O>::type, 1, k, KIND_NUM> cache;
    return cache.get(n, i, o, fwd, h);
}

template <kind k0 = KIND_NUM,
          kind k1 = KIND_NUM,
          typename I = void,
          typename O = void>
plan<typename io_traits<I, O>::type> &get_plan(const int n0,
                                               const int n1,
                                               const I *i,
                                               const O *o,
                                               const bool fwd,
                                               const how h = MEASURE)
{
    static plan_cache<typename io_traits<I, O>::type, 2, k0, k1> cache;
    return cache.get(n0, n1, i, o, fwd, h);
}
}

IEXP_NS_END

#endif /* __IEXP_FFT_FFTW_PLAN_CACHE__ */
