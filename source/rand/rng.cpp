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

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <rand/rng.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// internal macro
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// internal type
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// extern declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variant
////////////////////////////////////////////////////////////

static const gsl_rng_type *s_type_map[] = {
    gsl_rng_borosh13,
    gsl_rng_coveyou,
    gsl_rng_cmrg,
    gsl_rng_fishman18,
    gsl_rng_fishman20,
    gsl_rng_fishman2x,
    gsl_rng_gfsr4,
    gsl_rng_knuthran,
    gsl_rng_knuthran2,
    gsl_rng_knuthran2002,
    gsl_rng_lecuyer21,
    gsl_rng_minstd,
    gsl_rng_mrg,
    gsl_rng_mt19937,
    gsl_rng_mt19937_1999,
    gsl_rng_mt19937_1998,
    gsl_rng_r250,
    gsl_rng_ran0,
    gsl_rng_ran1,
    gsl_rng_ran2,
    gsl_rng_ran3,
    gsl_rng_rand,
    gsl_rng_rand48,
    gsl_rng_random128_bsd,
    gsl_rng_random128_glibc2,
    gsl_rng_random128_libc5,
    gsl_rng_random256_bsd,
    gsl_rng_random256_glibc2,
    gsl_rng_random256_libc5,
    gsl_rng_random32_bsd,
    gsl_rng_random32_glibc2,
    gsl_rng_random32_libc5,
    gsl_rng_random64_bsd,
    gsl_rng_random64_glibc2,
    gsl_rng_random64_libc5,
    gsl_rng_random8_bsd,
    gsl_rng_random8_glibc2,
    gsl_rng_random8_libc5,
    gsl_rng_random_bsd,
    gsl_rng_random_glibc2,
    gsl_rng_random_libc5,
    gsl_rng_randu,
    gsl_rng_ranf,
    gsl_rng_ranlux,
    gsl_rng_ranlux389,
    gsl_rng_ranlxd1,
    gsl_rng_ranlxd2,
    gsl_rng_ranlxs0,
    gsl_rng_ranlxs1,
    gsl_rng_ranlxs2,
    gsl_rng_ranmar,
    gsl_rng_slatec,
    gsl_rng_taus,
    gsl_rng_taus2,
    gsl_rng_taus113,
    gsl_rng_transputer,
    gsl_rng_tt800,
    gsl_rng_uni,
    gsl_rng_uni32,
    gsl_rng_vax,
    gsl_rng_waterman14,
    gsl_rng_zuf,
};

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface implementation
////////////////////////////////////////////////////////////

rng::rng(rng_type type, unsigned long seed)
{
    eigen_assert(type < RNG_TYPE_NUM);
    m_rng = gsl_rng_alloc(s_type_map[type]);
    IEXP_NOT_NULLPTR(m_rng);

    if (seed != 0) {
        gsl_rng_set(m_rng, seed);
    }
}
}

IEXP_NS_END
