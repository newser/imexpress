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

#ifndef __IEXP_DAE_ERROR__
#define __IEXP_DAE_ERROR__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define cv_check(expr)                                                         \
    do {                                                                       \
        int __e = expr;                                                        \
        if (__e < 0) {                                                         \
            iexp::dae::cv_throw(__e);                                          \
        }                                                                      \
    } while (0)

#define ls_check(expr)                                                         \
    do {                                                                       \
        int __e = expr;                                                        \
        if (__e < 0) {                                                         \
            iexp::dae::ls_throw(__e);                                          \
        }                                                                      \
    } while (0)

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

extern void cv_throw(int e);

extern void ls_throw(int e);
}

IEXP_NS_END

#endif /* __IEXP_DAE_ERROR__ */
