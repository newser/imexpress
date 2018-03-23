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

#ifndef __IEXP_STATS_QUANTILE__
#define __IEXP_STATS_QUANTILE__

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
inline double quantile_impl(const T data[], const size_t n, double f)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double quantile_impl(const double data[], const size_t n, double f)
{
    return gsl_stats_quantile_from_sorted_data(data, 1, n, f);
}

#define DEFINE_QUANTILE(type, name)                                            \
    template <>                                                                \
    inline double quantile_impl(const type data[], const size_t n, double f)   \
    {                                                                          \
        return gsl_stats_##name##_quantile_from_sorted_data(data, 1, n, f);    \
    }
DEFINE_QUANTILE(char, char)
DEFINE_QUANTILE(unsigned char, uchar)
DEFINE_QUANTILE(short, short)
DEFINE_QUANTILE(unsigned short, ushort)
DEFINE_QUANTILE(int, int)
DEFINE_QUANTILE(unsigned int, uint)
DEFINE_QUANTILE(float, float)
DEFINE_QUANTILE(long double, long_double)
#undef DEFINE_QUANTILE

template <typename T>
inline double quantile(const ArrayBase<T> &data, double f)
{
    eigen_assert(IS_VEC(data));

    typename type_eval<T>::type m_data(data.eval());
    return quantile_impl(m_data.data(), m_data.size(), f);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_STATS_QUANTILE__ */
