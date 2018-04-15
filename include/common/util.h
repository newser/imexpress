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

#define RETURN_NAN_OR_THROW(t) throw(t)

#define RETURN_OR_THROW(t) throw(t)

#define IEXP_NOT_COPYABLE(name)                                                \
    name(const name &) = delete;                                               \
    name &operator=(const name &) = delete;

#define ASSIGN_IJ_ROW(dst, src, src_expr)                                      \
    do {                                                                       \
        int __idx = 0;                                                         \
        for (int i = 0; i < (src).rows(); ++i) {                               \
            for (int j = 0; j < (src).cols(); ++j) {                           \
                (dst)[__idx++] = (src_expr);                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

#define ASSIGN_IJ_COL(dst, src, src_expr)                                      \
    do {                                                                       \
        int __idx = 0;                                                         \
        for (int j = 0; j < (src).cols(); ++j) {                               \
            for (int i = 0; i < (src).rows(); ++i) {                           \
                (dst)[__idx++] = (src_expr);                                   \
            }                                                                  \
        }                                                                      \
    } while (0)

#define ASSIGN_IJ(row_major, dst, src, src_expr)                               \
    do {                                                                       \
        if ((row_major)) {                                                     \
            ASSIGN_IJ_ROW(dst, src, src_expr);                                 \
        } else {                                                               \
            ASSIGN_IJ_COL(dst, src, src_expr);                                 \
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

IEXP_NS_END

#endif /* __IEXP_UTIL__ */
