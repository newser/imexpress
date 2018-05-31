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

#ifndef __IEXP_DAE_DEF__
#define __IEXP_DAE_DEF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <cvode/cvode.h>
#include <sundials/sundials_iterative.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

enum class precondition
{
    NONE = PREC_NONE,
    LEFT = PREC_LEFT,
    RIGHT = PREC_RIGHT,
};

enum class gram_schmidt
{
    MODIFIED = MODIFIED_GS,
    CLASSICAL = CLASSICAL_GS,
};

enum class multistep
{
    ADAMS = CV_ADAMS,
    BDF = CV_BDF,
};

enum class iteration
{
    FUNCTIONAL = CV_FUNCTIONAL,
    NEWTON = CV_NEWTON,
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_DEF__ */
