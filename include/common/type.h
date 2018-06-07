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

#define UNSUPPORTED_TYPE(t)                                                    \
    throw std::invalid_argument(std::string("type \"")                         \
                                    .append(typeid(t).name())                  \
                                    .append("\" is not implemented yet"))

#define SCALAR(t) std::conditional<is_complex<t>::value, t::value_type, t>::type

#define IS_MATRIX(t) is_matrix<t>::value

#define IS_ARRAY(t) is_array<t>::value

#define IS_DYNAMIC(t) is_dynamic<t>::value

#define TP1(t) typename Eigen::internal::traits<t>::Scalar

#define TP2(t) Eigen::internal::traits<t>::RowsAtCompileTime

#define TP3(t) Eigen::internal::traits<t>::ColsAtCompileTime

//#define TP4(t) Eigen::internal::traits<t>::Flags
#define TP4(t)                                                                 \
    (Eigen::internal::traits<t>::Flags & RowMajorBit) ? RowMajor : ColMajor

#define TP5(t) Eigen::internal::traits<t>::MaxRowsAtCompileTime

#define TP6(t) Eigen::internal::traits<t>::MaxColsAtCompileTime

#define TYPE_BOOL(expr)                                                        \
    typename TYPE_CHOOSE(expr, std::true_type, std::false_type)

#define M2V_NUM(type, x) ((TP4(T) == RowMajor) ? (x).rows() : (x).cols())

#define M2V_DIM(type, x) ((TP4(T) == RowMajor) ? (x).cols() : (x).rows())

#define DEFINE_DERIVED                                                         \
    Derived &derived()                                                         \
    {                                                                          \
        return *static_cast<Derived *>(this);                                  \
    }                                                                          \
                                                                               \
    const Derived &derived() const                                             \
    {                                                                          \
        return *static_cast<const Derived *>(this);                            \
    }

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
struct is_dynamic
{
    static const bool value =
        ((Eigen::internal::traits<T>::MaxRowsAtCompileTime == Dynamic) ||
         (Eigen::internal::traits<T>::MaxColsAtCompileTime == Dynamic));
};

template <typename T,
          typename _Scalar = TP1(T),
          int _Rows = TP2(T),
          int _Cols = TP3(T),
          int _Options = TP4(T),
          int _MaxRows = TP5(T),
          int _MaxCols = TP6(T)>
struct dense_derive
{
    using matrix = Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using array = Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
    using type = typename std::conditional<IS_MATRIX(T), matrix, array>::type;
};

using RowMatrixXd = Matrix<double, Dynamic, Dynamic, RowMajor>;

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_TYPE__ */
