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

#ifndef __IEXP_STATS_SKEWNESS__
#define __IEXP_STATS_SKEWNESS__

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
// skewness
// ========================================

template <typename T>
inline double skewness_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double skewness_impl(const double data[], size_t n)
{
    return gsl_stats_skew(data, 1, n);
}

#define DEFINE_SKEWNESS(type, name)                                            \
    template <>                                                                \
    inline double skewness_impl(const type data[], size_t n)                   \
    {                                                                          \
        return gsl_stats_##name##_skew(data, 1, n);                            \
    }
DEFINE_SKEWNESS(char, char)
DEFINE_SKEWNESS(unsigned char, uchar)
DEFINE_SKEWNESS(short, short)
DEFINE_SKEWNESS(unsigned short, ushort)
DEFINE_SKEWNESS(int, int)
DEFINE_SKEWNESS(unsigned int, uint)
DEFINE_SKEWNESS(float, float)
DEFINE_SKEWNESS(long double, long_double)
#undef DEFINE_SKEWNESS

template <typename T>
inline double skewness(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return skewness_impl(m_data.data(), m_data.size());
}

// ========================================
// skewness with specified mean and std
// ========================================

template <typename T>
inline double skewness_mv_impl(const T data[],
                               size_t n,
                               double mean,
                               double std)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double skewness_mv_impl(const double data[],
                               size_t n,
                               double mean,
                               double std)
{
    return gsl_stats_skew_m_sd(data, 1, n, mean, std);
}

#define DEFINE_SKEWNESS_MV(type, name)                                         \
    template <>                                                                \
    inline double skewness_mv_impl(const type data[],                          \
                                   size_t n,                                   \
                                   double mean,                                \
                                   double std)                                 \
    {                                                                          \
        return gsl_stats_##name##_skew_m_sd(data, 1, n, mean, std);            \
    }
DEFINE_SKEWNESS_MV(char, char)
DEFINE_SKEWNESS_MV(unsigned char, uchar)
DEFINE_SKEWNESS_MV(short, short)
DEFINE_SKEWNESS_MV(unsigned short, ushort)
DEFINE_SKEWNESS_MV(int, int)
DEFINE_SKEWNESS_MV(unsigned int, uint)
DEFINE_SKEWNESS_MV(float, float)
DEFINE_SKEWNESS_MV(long double, long_double)
#undef DEFINE_SKEWNESS_MV

template <typename T>
inline double skewness(const ArrayBase<T> &data, double mean, double std)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return skewness_mv_impl(m_data.data(), m_data.size(), mean, std);
}

////////////////////////////////////////////////////////////
// global skewnessiants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_SKEWNESS__ */
