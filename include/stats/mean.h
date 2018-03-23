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

#ifndef __IEXP_STATS_MEAN__
#define __IEXP_STATS_MEAN__

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
// mean
// ========================================

template <typename T>
inline double mean_impl(const T data[], const size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double mean_impl(const double data[], const size_t n)
{
    return gsl_stats_mean(data, 1, n);
}

#define DEFINE_MEAN(type, name)                                                \
    template <>                                                                \
    inline double mean_impl(const type data[], const size_t n)                 \
    {                                                                          \
        return gsl_stats_##name##_mean(data, 1, n);                            \
    }
DEFINE_MEAN(char, char)
DEFINE_MEAN(unsigned char, uchar)
DEFINE_MEAN(short, short)
DEFINE_MEAN(unsigned short, ushort)
DEFINE_MEAN(int, int)
DEFINE_MEAN(unsigned int, uint)
DEFINE_MEAN(float, float)
DEFINE_MEAN(long double, long_double)
#undef DEFINE_MEAN

template <typename T>
inline double mean(const ArrayBase<T> &data)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return mean_impl(m_data.data(), m_data.size());
}

// ========================================
// weighted mean
// ========================================

template <typename T>
inline double wmean_impl(const T data[], const T w[], const size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wmean_impl(const double data[], const double w[], const size_t n)
{
    return gsl_stats_wmean(w, 1, data, 1, n);
}

template <typename T>
inline double wmean(const ArrayBase<T> &data, const ArrayBase<T> &weight)
{
    eigen_assert(IS_VEC(data) && IS_VEC(data) &&
                 (data.size() == weight.size()));

    typename type_eval<T>::type m_data(data.eval()), m_w(weight.eval());
    return wmean_impl(m_data.data(), m_w.data(), m_data.size());
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_MEAN__ */
