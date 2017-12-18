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

#ifndef __IEXP_RAND__
#define __IEXP_RAND__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

template <typename T>
inline T rand_impl(const rng &r)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double rand_impl<double>(const rng &r)
{
    return r.uniform_double();
}

template <>
inline unsigned long rand_impl<unsigned long>(const rng &r)
{
    return r.uniform_ulong();
}

#define DEFINE_RAND_IMPL(t)                                                    \
    template <>                                                                \
    inline t rand_impl<t>(const rng &r)                                        \
    {                                                                          \
        return r.uniform_ulong(std::numeric_limits<t>::max() + 1);             \
    }
DEFINE_RAND_IMPL(long)
DEFINE_RAND_IMPL(unsigned int)
DEFINE_RAND_IMPL(int)
DEFINE_RAND_IMPL(unsigned short)
DEFINE_RAND_IMPL(short)
DEFINE_RAND_IMPL(unsigned char)
DEFINE_RAND_IMPL(char)
#undef DEFINE_RAND_IMPL

template <typename T>
inline auto rand(DenseBase<T> &x,
                 unsigned long seed = 0,
                 rng_type type = DEFAULT_RNG) -> decltype(x.derived())
{
    rng r(type, seed);

    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        data[i] = rand_impl<typename T::Scalar>(r);
    }

    return x.derived();
}

template <typename T>
inline auto rand(DenseBase<T> &x, const rng &r) -> decltype(x.derived())
{
    typename T::Scalar *data = x.derived().data();
    for (Index i = 0; i < x.size(); ++i) {
        data[i] = rand_impl<typename T::Scalar>(r);
    }

    return x.derived();
}
}

IEXP_NS_END

#endif /* __IEXP_RAND__ */
