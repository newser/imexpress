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

#ifndef __IEXP_INTEGRAL_TYPE__
#define __IEXP_INTEGRAL_TYPE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

<<<<<<< HEAD
#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
=======
#include <common/namespace.h>

#include <gsl/gsl_errno.h>
>>>>>>> c0b684af6635f4febd01a9e71af32e55409b2db2

IEXP_NS_BEGIN

namespace integral
{
    ////////////////////////////////////////////////////////////
    // macro definition
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // type definition
    ////////////////////////////////////////////////////////////

    enum key
    {
        GAUSS15 = GSL_INTEG_GAUSS15,
        GAUSS21 = GSL_INTEG_GAUSS21,
        GAUSS31 = GSL_INTEG_GAUSS31,
        GAUSS41 = GSL_INTEG_GAUSS41,
        GAUSS51 = GSL_INTEG_GAUSS51,
        GAUSS61 = GSL_INTEG_GAUSS61,
    };

    ////////////////////////////////////////////////////////////
    // global variants
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // interface declaration
    ////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_INTEGRAL_TYPE__ */
