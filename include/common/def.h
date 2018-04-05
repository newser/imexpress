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

#ifndef __IEXP_DEF__
#define __IEXP_DEF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <gsl/config.h>
#include <gsl/gsl_mode.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define DEFINE_TYPE_SUFFIX(def)                                                \
    def(char, char) def(unsigned char, uchar) def(short, short)                \
        def(unsigned short, ushort) def(int, int) def(unsigned int, uint)      \
            def(long, long) def(unsigned long, ulong) def(float, float)

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

enum class precision
{
    // Double-precision, a relative accuracy of approximately 2 * 10^{-16}.
    DOUBLE = GSL_PREC_DOUBLE,
    // Single-precision, a relative accuracy of approximately 10^{-7}.
    SINGLE = GSL_PREC_SINGLE,
    // Approximate values, a relative accuracy of approximately 5 * 10^{-4}.
    APPROX = GSL_PREC_APPROX,
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_DEF__ */
