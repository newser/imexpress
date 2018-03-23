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

#ifndef __IEXP_STATS_TSS__
#define __IEXP_STATS_TSS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_statistics.h>

IEXP_NS_BEGIN

namespace stats {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// TSS
// ========================================

template <typename T>
inline double tss_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double tss_impl(const double data[], size_t n)
{
    return gsl_stats_tss(data, 1, n);
}

#define DEFINE_TSS(type, name)                                                 \
    template <>                                                                \
    inline double tss_impl(const type data[], size_t n)                        \
    {                                                                          \
        return gsl_stats_##name##_sd(data, 1, n);                              \
    }
DEFINE_TSS(char, char)
DEFINE_TSS(unsigned char, uchar)
DEFINE_TSS(short, short)
DEFINE_TSS(unsigned short, ushort)
DEFINE_TSS(int, int)
DEFINE_TSS(unsigned int, uint)
DEFINE_TSS(float, float)
DEFINE_TSS(long double, long_double)
#undef DEFINE_TSS

template <typename T>
inline double tss(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return tss_impl(m_data.data(), m_data.size());
}

// ========================================
// TSS relative to specified mean
// ========================================

template <typename T>
inline double tss_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double tss_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_tss_m(data, 1, n, mean);
}

#define DEFINE_TSS_M(type, name)                                               \
    template <>                                                                \
    inline double tss_m_impl(const type data[], size_t n, double mean)         \
    {                                                                          \
        return gsl_stats_##name##_sd_m(data, 1, n, mean);                      \
    }
DEFINE_TSS_M(char, char)
DEFINE_TSS_M(unsigned char, uchar)
DEFINE_TSS_M(short, short)
DEFINE_TSS_M(unsigned short, ushort)
DEFINE_TSS_M(int, int)
DEFINE_TSS_M(unsigned int, uint)
DEFINE_TSS_M(float, float)
DEFINE_TSS_M(long double, long_double)
#undef DEFINE_TSS_M

template <typename T>
inline double tss(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return tss_m_impl(m_data.data(), m_data.size(), mean);
}

// ========================================
// weighted tss
// ========================================

template <typename T>
inline double wtss_impl(const T data[], const T w[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wtss_impl(const double data[], const double w[], size_t n)
{
    return gsl_stats_wtss(w, 1, data, 1, n);
}

template <typename T>
inline double wtss(const ArrayBase<T> &data, const ArrayBase<T> &weight)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wtss_impl(m_data.data(), m_w.data(), m_data.size());
}

// ========================================
// weighted tss relative to specified mean
// ========================================

template <typename T>
inline double wtss_impl(const T data[], const T w[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wtss_impl(const double data[],
                        const double w[],
                        size_t n,
                        double mean)
{
    return gsl_stats_wtss_m(w, 1, data, 1, n, mean);
}

template <typename T>
inline double wtss(const ArrayBase<T> &data,
                   const ArrayBase<T> &weight,
                   double mean)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wtss_impl(m_data.data(), m_w.data(), m_data.size(), mean);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_TSS__ */
