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

#ifndef __IEXP_UTIL__
#define __IEXP_UTIL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/namespace.h>

#include <Eigen/Dense>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define IEXP_NOT_NULLPTR(p)                                                    \
    do {                                                                       \
        if ((p) == nullptr) {                                                  \
            throw std::bad_alloc();                                            \
        }                                                                      \
    } while (0)

#define TYPE_IS(t1, t2) std::is_same<t1, t2>::value

#define TYPE_CHOOSE(v, t1, t2) std::conditional<v, t1, t2>::type

#define RETURN_NAN_OR_THROW(t) throw(t)

#define RETURN_OR_THROW(t) throw(t)

#define IS_VEC(v) (((v).rows() == 1) || ((v).cols() == 1))

#define VEC_SAME_SIZE(v1, v2)                                                  \
    (IS_VEC(v1) && IS_VEC(v2) && ((v1).size() == (v2).size()))

#define MATRIX_SAME_SIZE(m1, m2)                                               \
    (((m1).rows() == (m2).rows()) && ((m1).cols() == (m2).cols()))

#define UNSUPPORTED_TYPE(t) throw std::invalid_argument(typeid(t).name())

#define IS_COMPLEX(t) is_complex<t>::value

#define SCALAR(t) std::conditional<is_complex<t>::value, t::value_type, t>::type

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

/*
template <typename T, typename U>
struct type_is
{
    static const bool value = false;
};

template <typename T>
struct type_is<T, T>
{
    static const bool value = true;
};
*/

template <typename T>
struct type_eval
{
    typedef typename Eigen::internal::eval<T>::type type;
};

/*
template <typename T1, typename T2, bool V>
struct type_choose
{
    typedef void type;
};

template <typename T1, typename T2>
struct type_choose<T1, T2, true>
{
    typedef T1 type;
};

template <typename T1, typename T2>
struct type_choose<T1, T2, false>
{
    typedef T2 type;
};
*/

template <typename T>
struct is_complex
{
    static const bool value = std::is_same<T, std::complex<float>>::value ||
                              std::is_same<T, std::complex<double>>::value ||
                              std::is_same<T, std::complex<long double>>::value;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_UTIL__ */
