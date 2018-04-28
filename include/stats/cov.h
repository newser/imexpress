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

#ifndef __IEXP_STATS_COV__
#define __IEXP_STATS_COV__

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
// covariance
// ========================================

template <typename T>
inline double cov_impl(const T data1[], const T data2[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cov_impl(const double data1[], const double data2[], size_t n)
{
    return gsl_stats_covariance(data1, 1, data2, 1, n);
}

#define DEFINE_COV(type, name)                                                 \
    template <>                                                                \
    inline double cov_impl(const type data1[], const type data2[], size_t n)   \
    {                                                                          \
        return gsl_stats_##name##_covariance(data1, 1, data2, 1, n);           \
    }
DEFINE_COV(char, char)
DEFINE_COV(unsigned char, uchar)
DEFINE_COV(short, short)
DEFINE_COV(unsigned short, ushort)
DEFINE_COV(int, int)
DEFINE_COV(unsigned int, uint)
DEFINE_COV(float, float)
DEFINE_COV(long double, long_double)
#undef DEFINE_COV

template <typename T>
inline double cov(const DenseBase<T> &data1, const DenseBase<T> &data2)
{
    eigen_assert(data1.size() == data2.size());

    typename type_eval<T>::type m_data1(data1.eval()), m_data2(data2.eval());
    return cov_impl(m_data1.data(), m_data2.data(), m_data1.size());
}

// ========================================
// covariance
// ========================================

template <typename T>
inline double cov_impl(
    const T data1[], const T data2[], size_t n, double mean1, double mean2)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double cov_impl(const double data1[],
                       const double data2[],
                       size_t n,
                       double mean1,
                       double mean2)
{
    return gsl_stats_covariance_m(data1, 1, data2, 1, n, mean1, mean2);
}

#define DEFINE_COV_M(type, name)                                               \
    template <>                                                                \
    inline double cov_impl(const type data1[],                                 \
                           const type data2[],                                 \
                           size_t n,                                           \
                           double mean1,                                       \
                           double mean2)                                       \
    {                                                                          \
        return gsl_stats_##name##_covariance_m(data1,                          \
                                               1,                              \
                                               data2,                          \
                                               1,                              \
                                               n,                              \
                                               mean1,                          \
                                               mean2);                         \
    }
DEFINE_COV_M(char, char)
DEFINE_COV_M(unsigned char, uchar)
DEFINE_COV_M(short, short)
DEFINE_COV_M(unsigned short, ushort)
DEFINE_COV_M(int, int)
DEFINE_COV_M(unsigned int, uint)
DEFINE_COV_M(float, float)
DEFINE_COV_M(long double, long_double)
#undef DEFINE_COV_M

template <typename T>
inline double cov(const DenseBase<T> &data1,
                  const DenseBase<T> &data2,
                  double mean1,
                  double mean2)
{
    eigen_assert(data1.size() == data2.size());

    typename type_eval<T>::type m_data1(data1.eval()), m_data2(data2.eval());
    return cov_impl(m_data1.data(),
                    m_data2.data(),
                    m_data1.size(),
                    mean1,
                    mean2);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_COV__ */
