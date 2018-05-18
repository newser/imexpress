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

#ifndef __IEXP_TRIGONOMETRY__
#define __IEXP_TRIGONOMETRY__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>
#include <common/template_m2vnum_2d.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// sin
DEFINE_TEMPLATE_FOREACH(sin, double, double, gsl_sf_sin)
DEFINE_TEMPLATE_FOREACH_E(sin, double, double, gsl_sf_sin_e)

// cos
DEFINE_TEMPLATE_FOREACH(cos, double, double, gsl_sf_cos)
DEFINE_TEMPLATE_FOREACH_E(cos, double, double, gsl_sf_cos_e)

// hypot
DEFINE_TEMPLATE_M2VNUM_2D(hypot, double, double, double, gsl_sf_hypot)
DEFINE_TEMPLATE_M2VNUM_2D_E(hypot, double, double, double, gsl_sf_hypot_e)

// sinc
DEFINE_TEMPLATE_FOREACH(sinc, double, double, gsl_sf_sinc)
DEFINE_TEMPLATE_FOREACH_E(sinc, double, double, gsl_sf_sinc_e)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_TRIGONOMETRY__ */
