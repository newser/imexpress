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

#ifndef __IEXP_STATS_VAR__
#define __IEXP_STATS_VAR__

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
// variance
// ========================================

template <typename T>
inline double var_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double var_impl(const double data[], size_t n)
{
    return gsl_stats_variance(data, 1, n);
}

#define DEFINE_VAR(type, name)                                                 \
    template <>                                                                \
    inline double var_impl(const type data[], size_t n)                        \
    {                                                                          \
        return gsl_stats_##name##_variance(data, 1, n);                        \
    }
DEFINE_VAR(char, char)
DEFINE_VAR(unsigned char, uchar)
DEFINE_VAR(short, short)
DEFINE_VAR(unsigned short, ushort)
DEFINE_VAR(int, int)
DEFINE_VAR(unsigned int, uint)
DEFINE_VAR(float, float)
DEFINE_VAR(long double, long_double)
#undef DEFINE_VAR

template <typename T>
inline double var(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return var_impl(m_data.data(), m_data.size());
}

// ========================================
// variance relative to specified mean
// ========================================

template <typename T>
inline double var_m_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double var_m_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_variance_m(data, 1, n, mean);
}

#define DEFINE_VAR_M(type, name)                                               \
    template <>                                                                \
    inline double var_m_impl(const type data[], size_t n, double mean)         \
    {                                                                          \
        return gsl_stats_##name##_variance_m(data, 1, n, mean);                \
    }
DEFINE_VAR_M(char, char)
DEFINE_VAR_M(unsigned char, uchar)
DEFINE_VAR_M(short, short)
DEFINE_VAR_M(unsigned short, ushort)
DEFINE_VAR_M(int, int)
DEFINE_VAR_M(unsigned int, uint)
DEFINE_VAR_M(float, float)
DEFINE_VAR_M(long double, long_double)
#undef DEFINE_VAR_M

template <typename T>
inline double var(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return var_m_impl(m_data.data(), m_data.size(), mean);
}

// ========================================
// unbiased variance
// ========================================

template <typename T>
inline double uvar_impl(const T data[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double uvar_impl(const double data[], size_t n, double mean)
{
    return gsl_stats_variance_with_fixed_mean(data, 1, n, mean);
}

#define DEFINE_UVAR_M(type, name)                                              \
    template <>                                                                \
    inline double uvar_impl(const type data[], size_t n, double mean)          \
    {                                                                          \
        return gsl_stats_##name##_variance_with_fixed_mean(data, 1, n, mean);  \
    }
DEFINE_UVAR_M(char, char)
DEFINE_UVAR_M(unsigned char, uchar)
DEFINE_UVAR_M(short, short)
DEFINE_UVAR_M(unsigned short, ushort)
DEFINE_UVAR_M(int, int)
DEFINE_UVAR_M(unsigned int, uint)
DEFINE_UVAR_M(float, float)
DEFINE_UVAR_M(long double, long_double)
#undef DEFINE_UVAR_M

template <typename T>
inline double uvar(const ArrayBase<T> &data, double mean)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return uvar_impl(m_data.data(), m_data.size(), mean);
}

// ========================================
// weighted variance
// ========================================

template <typename T>
inline double wvar_impl(const T data[], const T w[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wvar_impl(const double data[], const double w[], size_t n)
{
    return gsl_stats_wvariance(w, 1, data, 1, n);
}

template <typename T>
inline double wvar(const ArrayBase<T> &data, const ArrayBase<T> &weight)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wvar_impl(m_data.data(), m_w.data(), m_data.size());
}

// ========================================
// weighted variance relative to specified mean
// ========================================

template <typename T>
inline double wvar_impl(const T data[], const T w[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wvar_impl(const double data[],
                        const double w[],
                        size_t n,
                        double mean)
{
    return gsl_stats_wvariance_m(w, 1, data, 1, n, mean);
}

template <typename T>
inline double wvar(const ArrayBase<T> &data,
                   const ArrayBase<T> &weight,
                   double mean)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wvar_impl(m_data.data(), m_w.data(), m_data.size(), mean);
}

// ========================================
// unbiased weighted variance
// ========================================

template <typename T>
inline double wuvar_impl(const T data[], const T w[], size_t n, double mean)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wuvar_impl(const double data[],
                         const double w[],
                         size_t n,
                         double mean)
{
    return gsl_stats_wvariance_with_fixed_mean(w, 1, data, 1, n, mean);
}

template <typename T>
inline double wuvar(const ArrayBase<T> &data,
                    const ArrayBase<T> &weight,
                    double mean)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wuvar_impl(m_data.data(), m_w.data(), m_data.size(), mean);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_VAR__ */
