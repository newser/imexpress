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

#ifndef __IEXP_TYPE__
#define __IEXP_TYPE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <complex>
#include <limits>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define IS_INTEGER(t) std::numeric_limits<t>::is_integer

#define IS_COMPLEX(t) is_complex<t>::value

#define TYPE_IS(t1, t2) std::is_same<t1, t2>::value

#define TYPE_CHOOSE(v, t1, t2) std::conditional<v, t1, t2>::type

#define DEFINE_TYPE_SUFFIX(def)                                                \
    def(char, char) def(unsigned char, uchar) def(short, short)                \
        def(unsigned short, ushort) def(int, int) def(unsigned int, uint)      \
            def(long, long) def(unsigned long, ulong) def(float, float)

#define UNSUPPORTED_TYPE(t) throw std::invalid_argument(typeid(t).name())

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
struct is_complex
{
    static const bool value = std::is_same<T, std::complex<float>>::value ||
                              std::is_same<T, std::complex<double>>::value ||
                              std::is_same<T, std::complex<long double>>::value;
};

template <typename T>
struct type_eval
{
    typedef typename Eigen::internal::eval<T>::type type;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_TYPE__ */
