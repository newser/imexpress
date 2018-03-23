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

#ifndef __IEXP_STATS_STD__
#define __IEXP_STATS_STD__

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
// standard deviation
// ========================================

template <typename T>
inline double std_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double std_impl(const double data[], size_t n)
{
    return gsl_stats_sd(data, 1, n);
}

#define DEFINE_STD(type, name)                                                 \
    template <>                                                                \
    inline double std_impl(const type data[], size_t n)                        \
    {                                                                          \
        return gsl_stats_##name##_sd(data, 1, n);                              \
    }
DEFINE_STD(char, char)
DEFINE_STD(unsigned char, uchar)
DEFINE_STD(short, short)
DEFINE_STD(unsigned short, ushort)
DEFINE_STD(int, int)
DEFINE_STD(unsigned int, uint)
DEFINE_STD(float, float)
DEFINE_STD(long double, long_double)
#undef DEFINE_STD

template <typename T>
inline double std(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return std_impl(m_data.data(), m_data.size());
}

// ========================================
// standard deviation relative to specified mean
// ========================================

template <typename T>
inline double std_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double std_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_sd_m(data, 1, n, mean);
}

#define DEFINE_STD_M(type, name)                                               \
    template <>                                                                \
    inline double std_m_impl(const type data[], size_t n, double mean)         \
    {                                                                          \
        return gsl_stats_##name##_sd_m(data, 1, n, mean);                      \
    }
DEFINE_STD_M(char, char)
DEFINE_STD_M(unsigned char, uchar)
DEFINE_STD_M(short, short)
DEFINE_STD_M(unsigned short, ushort)
DEFINE_STD_M(int, int)
DEFINE_STD_M(unsigned int, uint)
DEFINE_STD_M(float, float)
DEFINE_STD_M(long double, long_double)
#undef DEFINE_STD_M

template <typename T>
inline double std(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return std_m_impl(m_data.data(), m_data.size(), mean);
}

// ========================================
// unbiased standard deviation
// ========================================

template <typename T>
inline double ustd_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double ustd_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_sd_with_fixed_mean(data, 1, n, mean);
}

#define DEFINE_Uustd_M(type, name)                                             \
    template <>                                                                \
    inline double ustd_m_impl(const type data[], size_t n, double mean)        \
    {                                                                          \
        return gsl_stats_##name##_sd_with_fixed_mean(data, 1, n, mean);        \
    }
DEFINE_Uustd_M(char, char) DEFINE_Uustd_M(unsigned char, uchar)
    DEFINE_Uustd_M(short, short) DEFINE_Uustd_M(unsigned short, ushort)
        DEFINE_Uustd_M(int, int) DEFINE_Uustd_M(unsigned int, uint)
            DEFINE_Uustd_M(float, float)
                DEFINE_Uustd_M(long double, long_double)
#undef DEFINE_Uustd_M

                    template <typename T>
                    inline double ustd(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return ustd_m_impl(m_data.data(), m_data.size(), mean);
}

// ========================================
// absolute standard deviation
// ========================================

template <typename T>
inline double abstd_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double abstd_impl(const double data[], size_t n)
{
    return gsl_stats_absdev(data, 1, n);
}

#define DEFINE_ABSTD(type, name)                                               \
    template <>                                                                \
    inline double abstd_impl(const type data[], size_t n)                      \
    {                                                                          \
        return gsl_stats_##name##_absdev(data, 1, n);                          \
    }
DEFINE_ABSTD(char, char)
DEFINE_ABSTD(unsigned char, uchar)
DEFINE_ABSTD(short, short)
DEFINE_ABSTD(unsigned short, ushort)
DEFINE_ABSTD(int, int)
DEFINE_ABSTD(unsigned int, uint)
DEFINE_ABSTD(float, float)
DEFINE_ABSTD(long double, long_double)
#undef DEFINE_ABSTD

template <typename T>
inline double abstd(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return abstd_impl(m_data.data(), m_data.size());
}

// ========================================
// absolute standard deviation relative to specified mean
// ========================================

template <typename T>
inline double abstd_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double abstd_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_absdev_m(data, 1, n, mean);
}

#define DEFINE_ABSTD_M(type, name)                                             \
    template <>                                                                \
    inline double abstd_m_impl(const type data[], size_t n, double mean)       \
    {                                                                          \
        return gsl_stats_##name##_absdev_m(data, 1, n, mean);                  \
    }
DEFINE_ABSTD_M(char, char)
DEFINE_ABSTD_M(unsigned char, uchar)
DEFINE_ABSTD_M(short, short)
DEFINE_ABSTD_M(unsigned short, ushort)
DEFINE_ABSTD_M(int, int)
DEFINE_ABSTD_M(unsigned int, uint)
DEFINE_ABSTD_M(float, float)
DEFINE_ABSTD_M(long double, long_double)
#undef DEFINE_ABSTD_M

template <typename T>
inline double abstd(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return abstd_m_impl(m_data.data(), m_data.size(), mean);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_STD__ */
