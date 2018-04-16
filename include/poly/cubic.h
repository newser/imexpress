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

#ifndef __IEXP_POLY_CUBIC__
#define __IEXP_POLY_CUBIC__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <math/inf_nan.h>

#include <gsl/gsl_poly.h>

IEXP_NS_BEGIN

namespace poly {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// solve in real field
// ========================================

inline int solve_cubic(
    double a, double b, double c, double d, double &x0, double &x1, double &x2)
{
    eigen_assert(a != 0);
    b /= a;
    c /= a;
    d /= a;

    int r = gsl_poly_solve_cubic(b, c, d, &x0, &x1, &x2);
    if (r == 0) {
        x0 = NAN;
        x1 = NAN;
        x2 = NAN;
    } else if (r == 1) {
        x1 = NAN;
        x2 = NAN;
    } else if (r == 2) {
        x2 = NAN;
    }
    return r;
}

inline std::tuple<double, double, double> solve_cubic(double a,
                                                      double b,
                                                      double c,
                                                      double d)
{
    double x0, x1, x2;
    solve_cubic(a, b, c, d, x0, x1, x2);
    return std::make_tuple(x0, x1, x2);
}

// ========================================
// solve in complex field
// ========================================

inline int complex_solve_cubic(double a,
                               double b,
                               double c,
                               double d,
                               std::complex<double> &x0,
                               std::complex<double> &x1,
                               std::complex<double> &x2)
{
    eigen_assert(a != 0);
    b /= a;
    c /= a;
    d /= a;

    int r = gsl_poly_complex_solve_cubic(b,
                                         c,
                                         d,
                                         (gsl_complex *)&x0,
                                         (gsl_complex *)&x1,
                                         (gsl_complex *)&x2);
    if (r == 0) {
        x0 = std::complex<double>(NAN, NAN);
        x1 = std::complex<double>(NAN, NAN);
        x2 = std::complex<double>(NAN, NAN);
    } else if (r == 1) {
        x1 = std::complex<double>(NAN, NAN);
        x2 = std::complex<double>(NAN, NAN);
    } else if (r == 2) {
        x2 = std::complex<double>(NAN, NAN);
    }
    return r;
}

inline std::
    tuple<std::complex<double>, std::complex<double>, std::complex<double>>
    complex_solve_cubic(double a, double b, double c, double d)
{
    std::complex<double> x0, x1, x2;
    complex_solve_cubic(a, b, c, d, x0, x1, x2);
    return std::make_tuple(x0, x1, x2);
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_POLY_CUBIC__ */
