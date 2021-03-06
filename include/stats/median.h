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

#ifndef __IEXP_STATS_MEDIAN__
#define __IEXP_STATS_MEDIAN__

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

template <typename T>
inline double median_impl(const T data[], size_t n)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double median_impl(const double data[], size_t n)
{
    return gsl_stats_median_from_sorted_data(data, 1, n);
}

#define DEFINE_MEAN(type, name)                                                \
    template <>                                                                \
    inline double median_impl(const type data[], size_t n)                     \
    {                                                                          \
        return gsl_stats_##name##_median_from_sorted_data(data, 1, n);         \
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
inline double median(const DenseBase<T> &data)
{
    typename type_eval<T>::type m_data(data.eval());
    return median_impl(m_data.data(), m_data.size());
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_MEDIAN__ */
