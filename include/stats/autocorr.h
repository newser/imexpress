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

#ifndef __IEXP_STATS_AUTOCORR__
#define __IEXP_STATS_AUTOCORR__

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
// autocorrelation
// ========================================

template <typename T>
inline double autocorr_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double autocorr_impl(const double data[], size_t n)
{
    return gsl_stats_lag1_autocorrelation(data, 1, n);
}

#define DEFINE_AUTOCORR(type, name)                                            \
    template <>                                                                \
    inline double autocorr_impl(const type data[], size_t n)                   \
    {                                                                          \
        return gsl_stats_##name##_lag1_autocorrelation(data, 1, n);            \
    }
DEFINE_AUTOCORR(char, char)
DEFINE_AUTOCORR(unsigned char, uchar)
DEFINE_AUTOCORR(short, short)
DEFINE_AUTOCORR(unsigned short, ushort)
DEFINE_AUTOCORR(int, int)
DEFINE_AUTOCORR(unsigned int, uint)
DEFINE_AUTOCORR(float, float)
DEFINE_AUTOCORR(long double, long_double)
#undef DEFINE_AUTOCORR

template <typename T>
inline double autocorr(const DenseBase<T> &data)
{
    typename type_eval<T>::type m_data(data.eval());
    return autocorr_impl(m_data.data(), m_data.size());
}

// ========================================
// autocorrelation with specified mean
// ========================================

template <typename T>
inline double autocorr_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double autocorr_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_lag1_autocorrelation_m(data, 1, n, mean);
}

#define DEFINE_AUTOCORR_M(type, name)                                          \
    template <>                                                                \
    inline double autocorr_m_impl(const type data[], size_t n, double mean)    \
    {                                                                          \
        return gsl_stats_##name##_lag1_autocorrelation_m(data, 1, n, mean);    \
    }
DEFINE_AUTOCORR_M(char, char)
DEFINE_AUTOCORR_M(unsigned char, uchar)
DEFINE_AUTOCORR_M(short, short)
DEFINE_AUTOCORR_M(unsigned short, ushort)
DEFINE_AUTOCORR_M(int, int)
DEFINE_AUTOCORR_M(unsigned int, uint)
DEFINE_AUTOCORR_M(float, float)
DEFINE_AUTOCORR_M(long double, long_double)
#undef DEFINE_AUTOCORR_M

template <typename T>
inline double autocorr(const DenseBase<T> &data, double mean)
{
    typename type_eval<T>::type m_data(data.eval());
    return autocorr_m_impl(m_data.data(), m_data.size(), mean);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_AUTOCORR__ */
