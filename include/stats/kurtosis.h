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

#ifndef __IEXP_STATS_KURTOSIS__
#define __IEXP_STATS_KURTOSIS__

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
// kurtosis
// ========================================

template <typename T>
inline double kurtosis_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double kurtosis_impl(const double data[], size_t n)
{
    return gsl_stats_kurtosis(data, 1, n);
}

#define DEFINE_KURTOSIS(type, name)                                            \
    template <>                                                                \
    inline double kurtosis_impl(const type data[], size_t n)                   \
    {                                                                          \
        return gsl_stats_##name##_kurtosis(data, 1, n);                        \
    }
DEFINE_KURTOSIS(char, char)
DEFINE_KURTOSIS(unsigned char, uchar)
DEFINE_KURTOSIS(short, short)
DEFINE_KURTOSIS(unsigned short, ushort)
DEFINE_KURTOSIS(int, int)
DEFINE_KURTOSIS(unsigned int, uint)
DEFINE_KURTOSIS(float, float)
DEFINE_KURTOSIS(long double, long_double)
#undef DEFINE_KURTOSIS

template <typename T>
inline double kurtosis(const DenseBase<T> &data)
{
    typename type_eval<T>::type m_data(data.eval());
    return kurtosis_impl(m_data.data(), m_data.size());
}

// ========================================
// kurtosis with specified mean and std
// ========================================

template <typename T>
inline double kurtosis_mv_impl(const T data[],
                               size_t n,
                               double mean,
                               double std)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double kurtosis_mv_impl(const double data[],
                               size_t n,
                               double mean,
                               double std)
{
    return gsl_stats_kurtosis_m_sd(data, 1, n, mean, std);
}

#define DEFINE_KURTOSIS_MV(type, name)                                         \
    template <>                                                                \
    inline double kurtosis_mv_impl(const type data[],                          \
                                   size_t n,                                   \
                                   double mean,                                \
                                   double std)                                 \
    {                                                                          \
        return gsl_stats_##name##_kurtosis_m_sd(data, 1, n, mean, std);        \
    }
DEFINE_KURTOSIS_MV(char, char)
DEFINE_KURTOSIS_MV(unsigned char, uchar)
DEFINE_KURTOSIS_MV(short, short)
DEFINE_KURTOSIS_MV(unsigned short, ushort)
DEFINE_KURTOSIS_MV(int, int)
DEFINE_KURTOSIS_MV(unsigned int, uint)
DEFINE_KURTOSIS_MV(float, float)
DEFINE_KURTOSIS_MV(long double, long_double)
#undef DEFINE_KURTOSIS_MV

template <typename T>
inline double kurtosis(const DenseBase<T> &data, double mean, double std)
{
    typename type_eval<T>::type m_data(data.eval());
    return kurtosis_mv_impl(m_data.data(), m_data.size(), mean, std);
}

////////////////////////////////////////////////////////////
// global kurtosisiants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_KURTOSIS__ */
