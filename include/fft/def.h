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

#ifndef __IEXP_FFT_DEF__
#define __IEXP_FFT_DEF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

IEXP_NS_BEGIN

namespace fft {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define IS_DCT(k) (((k) >= DCT_I) && ((k) <= DCT_IV))

#define IS_DST(k) (((k) >= DST_I) && ((k) <= DST_IV))

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

enum kind
{
    DCT_I,
    DCT_II,
    DCT_III,
    DCT_IV,
    DST_I,
    DST_II,
    DST_III,
    DST_IV,

    KIND_NUM
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_FFT_DEF__ */
