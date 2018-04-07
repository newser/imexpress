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

#define IS_LREF(t) std::is_lvalue_reference<t>::value

#define IS_RREF(t) std::is_rvalue_reference<t>::value

#define TYPE_IS(t1, t2) std::is_same<t1, t2>::value

#define TYPE_CHOOSE(v, t1, t2) std::conditional<v, t1, t2>::type

#define DEFINE_TYPE_SUFFIX(def)                                                \
    def(char, char) def(unsigned char, uchar) def(short, short)                \
        def(unsigned short, ushort) def(int, int) def(unsigned int, uint)      \
            def(long, long) def(unsigned long, ulong) def(float, float)

#define UNSUPPORTED_TYPE(t) throw std::invalid_argument(typeid(t).name())

#define SCALAR(t) std::conditional<is_complex<t>::value, t::value_type, t>::type

#define IS_MATRIX(t) is_matrix<t>::value

#define IS_ARRAY(t) is_array<t>::value

#define TP1_OF(t) typename Eigen::internal::traits<t>::Scalar

#define TP2_OF(t) Eigen::internal::traits<t>::RowsAtCompileTime

#define TP3_OF(t) Eigen::internal::traits<t>::ColsAtCompileTime

//#define TP4_OF(t) Eigen::internal::traits<t>::Flags
#define TP4_OF(t)                                                              \
    (Eigen::internal::traits<t>::Flags & RowMajorBit) ? RowMajor : ColMajor

#define TP5_OF(t) Eigen::internal::traits<t>::MaxRowsAtCompileTime

#define TP6_OF(t) Eigen::internal::traits<t>::MaxColsAtCompileTime

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
    using type = typename Eigen::internal::eval<T>::type;
};

template <typename T>
struct is_matrix
{
    static const bool value =
        std::is_same<typename Eigen::internal::traits<T>::XprKind,
                     MatrixXpr>::value;
};

template <typename T>
struct is_array
{
    static const bool value =
        std::is_same<typename Eigen::internal::traits<T>::XprKind,
                     ArrayXpr>::value;
};

template <typename T>
struct dense_derive
{
    using matrix = Matrix<TP1_OF(T),
                          TP2_OF(T),
                          TP3_OF(T),
                          TP4_OF(T),
                          TP5_OF(T),
                          TP6_OF(T)>;
    using array =
        Array<TP1_OF(T), TP2_OF(T), TP3_OF(T), TP4_OF(T), TP5_OF(T), TP6_OF(T)>;
    using type = typename std::conditional<IS_MATRIX(T), matrix, array>::type;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_TYPE__ */
