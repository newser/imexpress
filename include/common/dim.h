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

#ifndef __IEXP_DIM__
#define __IEXP_DIM__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define IS_VEC(v) (((v).rows() == 1) || ((v).cols() == 1))

#define VEC_SAME_SIZE(v1, v2)                                                  \
    (IS_VEC(v1) && IS_VEC(v2) && ((v1).size() == (v2).size()))

#define MATRIX_SAME_SIZE(m1, m2)                                               \
    (((m1).rows() == (m2).rows()) && ((m1).cols() == (m2).cols()))

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_DIM__ */
