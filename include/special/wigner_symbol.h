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

#ifndef __IEXP_WIGNER_SYMBOL__
#define __IEXP_WIGNER_SYMBOL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coupling.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// wigner 3j symbol
// ========================================

template <typename T>
inline double wigner3j_impl(
    const T ja, const T jb, const T jc, const T ma, const T mb, const T mc)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner3j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int ma,
                            const int mb,
                            const int mc)
{
    return gsl_sf_coupling_3j(ja, jb, jc, ma, mb, mc);
}

template <typename T>
double wigner3j(const ArrayBase<T> &c)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 2) && (c.cols() == 3));

    return wigner3j_impl(c(0, 0), c(0, 1), c(0, 2), c(1, 0), c(1, 1), c(1, 2));
}

template <typename T, typename U>
inline double wigner3j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T ma,
                              const T mb,
                              const T mc,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner3j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int ma,
                              const int mb,
                              const int mc,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_3j_e(ja, jb, jc, ma, mb, mc, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("wigner3j"));
}

template <typename T>
double wigner3j(const ArrayBase<T> &c, double &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 2) && (c.cols() == 3));

    return wigner3j_e_impl(c(0, 0),
                           c(0, 1),
                           c(0, 2),
                           c(1, 0),
                           c(1, 1),
                           c(1, 2),
                           e);
}

// ========================================
// wigner 6j syjeol
// ========================================

template <typename T>
inline double wigner6j_impl(
    const T ja, const T jb, const T jc, const T jd, const T je, const T jf)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner6j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int jd,
                            const int je,
                            const int jf)
{
    return gsl_sf_coupling_6j(ja, jb, jc, jd, je, jf);
}

template <typename T>
double wigner6j(const ArrayBase<T> &c)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 2) && (c.cols() == 3));

    return wigner6j_impl(c(0, 0), c(0, 1), c(0, 2), c(1, 0), c(1, 1), c(1, 2));
}

template <typename T, typename U>
inline double wigner6j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T jd,
                              const T je,
                              const T jf,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner6j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int jd,
                              const int je,
                              const int jf,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_6j_e(ja, jb, jc, jd, je, jf, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("wigner6j"));
}

template <typename T>
double wigner6j(const ArrayBase<T> &c, double &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 2) && (c.cols() == 3));

    return wigner6j_e_impl(c(0, 0),
                           c(0, 1),
                           c(0, 2),
                           c(1, 0),
                           c(1, 1),
                           c(1, 2),
                           e);
}

// ========================================
// wigner 9j syjeol
// ========================================

template <typename T>
inline double wigner9j_impl(const T ja,
                            const T jb,
                            const T jc,
                            const T jd,
                            const T je,
                            const T jf,
                            const T jg,
                            const T jh,
                            const T ji)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner9j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int jd,
                            const int je,
                            const int jf,
                            const int jg,
                            const int jh,
                            const int ji)
{
    return gsl_sf_coupling_9j(ja, jb, jc, jd, je, jf, jg, jh, ji);
}

template <typename T>
double wigner9j(const ArrayBase<T> &c)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 3) && (c.cols() == 3));

    return wigner9j_impl(c(0, 0),
                         c(0, 1),
                         c(0, 2),
                         c(1, 0),
                         c(1, 1),
                         c(1, 2),
                         c(2, 0),
                         c(2, 1),
                         c(2, 2));
}

template <typename T, typename U>
inline double wigner9j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T jd,
                              const T je,
                              const T jf,
                              const T jg,
                              const T jh,
                              const T ji,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner9j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int jd,
                              const int je,
                              const int jf,
                              const int jg,
                              const int jh,
                              const int ji,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_9j_e(ja, jb, jc, jd, je, jf, jg, jh, ji, &r) ==
        GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("wigner9j"));
}

template <typename T>
double wigner9j(const ArrayBase<T> &c, double &e)
{
    static_assert(TYPE_IS(typename T::Scalar, int) ||
                      TYPE_IS(typename T::Scalar, unsigned int),
                  "c can only be int or uint array");
    eigen_assert((c.rows() == 3) && (c.cols() == 3));

    return wigner9j_e_impl(c(0, 0),
                           c(0, 1),
                           c(0, 2),
                           c(1, 0),
                           c(1, 1),
                           c(1, 2),
                           c(2, 0),
                           c(2, 1),
                           c(2, 2),
                           e);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_WIGNER_SYMBOL__ */
