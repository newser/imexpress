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

#ifndef __IEXP_STATS_CORRCOEF__
#define __IEXP_STATS_CORRCOEF__

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
// correlation coefficient
// ========================================

template <typename T>
inline double corrcoef_impl(const T data1[], const T data2[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double corrcoef_impl(const double data1[],
                            const double data2[],
                            size_t n)
{
    return gsl_stats_correlation(data1, 1, data2, 1, n);
}

#define DEFINE_CORRCOEF(type, name)                                            \
    template <>                                                                \
    inline double corrcoef_impl(const type data1[],                            \
                                const type data2[],                            \
                                size_t n)                                      \
    {                                                                          \
        return gsl_stats_##name##_correlation(data1, 1, data2, 1, n);          \
    }
DEFINE_CORRCOEF(char, char)
DEFINE_CORRCOEF(unsigned char, uchar)
DEFINE_CORRCOEF(short, short)
DEFINE_CORRCOEF(unsigned short, ushort)
DEFINE_CORRCOEF(int, int)
DEFINE_CORRCOEF(unsigned int, uint)
DEFINE_CORRCOEF(float, float)
DEFINE_CORRCOEF(long double, long_double)
#undef DEFINE_CORRCOEF

template <typename T>
inline double corrcoef(const DenseBase<T> &data1, const DenseBase<T> &data2)
{
    eigen_assert(data1.size() == data2.size());

    typename type_eval<T>::type m_data1(data1.eval()), m_data2(data2.eval());
    return corrcoef_impl(m_data1.data(), m_data2.data(), m_data1.size());
}

// ========================================
// spearman rand correlation
// ========================================

template <typename T>
inline double spearman_impl(const T data1[],
                            const T data2[],
                            size_t n,
                            double work[])
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double spearman_impl(const double data1[],
                            const double data2[],
                            size_t n,
                            double work[])
{
    return gsl_stats_spearman(data1, 1, data2, 1, n, work);
}

#define DEFINE_SPEARMAN(type, name)                                            \
    template <>                                                                \
    inline double spearman_impl(const type data1[],                            \
                                const type data2[],                            \
                                size_t n,                                      \
                                double work[])                                 \
    {                                                                          \
        return gsl_stats_##name##_spearman(data1, 1, data2, 1, n, work);       \
    }
DEFINE_SPEARMAN(char, char)
DEFINE_SPEARMAN(unsigned char, uchar)
DEFINE_SPEARMAN(short, short)
DEFINE_SPEARMAN(unsigned short, ushort)
DEFINE_SPEARMAN(int, int)
DEFINE_SPEARMAN(unsigned int, uint)
DEFINE_SPEARMAN(float, float)
DEFINE_SPEARMAN(long double, long_double)
#undef DEFINE_SPEARMAN

template <typename T>
inline double spearman(const DenseBase<T> &data1, const DenseBase<T> &data2)
{
    eigen_assert(data1.size() == data2.size());

    typename type_eval<T>::type m_data1(data1.eval()), m_data2(data2.eval());

    std::unique_ptr<double> work(new double[data1.size() * 2]);
    return spearman_impl(m_data1.data(),
                         m_data2.data(),
                         m_data1.size(),
                         work.get());
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_CORRCOEF__ */
