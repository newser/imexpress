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

#ifndef __IEXP_SF_TRANSPORT__
#define __IEXP_SF_TRANSPORT__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_transport.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// transport2
DEFINE_TEMPLATE_FOREACH(transport2, double, double, gsl_sf_transport_2)

DEFINE_TEMPLATE_FOREACH_E(transport2, double, double, gsl_sf_transport_2_e)

// transport3
DEFINE_TEMPLATE_FOREACH(transport3, double, double, gsl_sf_transport_3)

DEFINE_TEMPLATE_FOREACH_E(transport3, double, double, gsl_sf_transport_3_e)

// transport_4
DEFINE_TEMPLATE_FOREACH(transport4, double, double, gsl_sf_transport_4)

DEFINE_TEMPLATE_FOREACH_E(transport4, double, double, gsl_sf_transport_4_e)

// transport_5
DEFINE_TEMPLATE_FOREACH(transport5, double, double, gsl_sf_transport_5)

DEFINE_TEMPLATE_FOREACH_E(transport5, double, double, gsl_sf_transport_5_e)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_SF_TRANSPORT__ */
