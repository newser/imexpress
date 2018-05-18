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

#ifndef __IEXP_HYDROGENIC_R__
#define __IEXP_HYDROGENIC_R__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_m2vnum_2d.h>
#include <common/template_m2vnum_4d.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// hydrogenicR_1
// ========================================

DEFINE_TEMPLATE_M2VNUM_2D(
    hydroR1, double, double, double, double, gsl_sf_hydrogenicR_1)

DEFINE_TEMPLATE_M2VNUM_2D_E(
    hydroR1, double, double, double, double, gsl_sf_hydrogenicR_1_e)

// ========================================
// hydrogenicR
// ========================================

DEFINE_TEMPLATE_M2VNUM_4D(
    hydroR, double, double, int, int, double, double, gsl_sf_hydrogenicR)

DEFINE_TEMPLATE_M2VNUM_4D_E(
    hydroR, double, double, int, int, double, double, gsl_sf_hydrogenicR_e)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_HYDROGENIC_R__ */
