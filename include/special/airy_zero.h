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

#ifndef __IEXP_SPECIAL_AIRY_ZERO__
#define __IEXP_SPECIAL_AIRY_ZERO__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/template_foreach.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// location of n-th 0 of airy Ai
// ========================================

DEFINE_TEMPLATE_FOREACH(airy_n0_Ai, double, unsigned int, gsl_sf_airy_zero_Ai)

DEFINE_TEMPLATE_FOREACH_E(airy_n0_Ai,
                          double,
                          unsigned int,
                          gsl_sf_airy_zero_Ai_e)

// ========================================
// location of n-th 0 of airy Ai derivative
// ========================================

DEFINE_TEMPLATE_FOREACH(airy_n0_Ai_deriv,
                        double,
                        unsigned int,
                        gsl_sf_airy_zero_Ai_deriv)

DEFINE_TEMPLATE_FOREACH_E(airy_n0_Ai_deriv,
                          double,
                          unsigned int,
                          gsl_sf_airy_zero_Ai_deriv_e)

// ========================================
// location of n-th 0 of airy Bi
// ========================================

DEFINE_TEMPLATE_FOREACH(airy_n0_Bi, double, unsigned int, gsl_sf_airy_zero_Bi)

DEFINE_TEMPLATE_FOREACH_E(airy_n0_Bi,
                          double,
                          unsigned int,
                          gsl_sf_airy_zero_Bi_e)

// ========================================
// location of n-th 0 of airy Bi derivative
// ========================================

DEFINE_TEMPLATE_FOREACH(airy_n0_Bi_deriv,
                        double,
                        unsigned int,
                        gsl_sf_airy_zero_Bi_deriv)

DEFINE_TEMPLATE_FOREACH_E(airy_n0_Bi_deriv,
                          double,
                          unsigned int,
                          gsl_sf_airy_zero_Bi_deriv_e)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_SPECIAL_AIRY_ZERO__ */
