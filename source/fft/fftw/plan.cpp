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

#include <fft/fftw/plan.h>

IEXP_NS_BEGIN

namespace fftw3 {

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

fftw_r2r_kind fwd_kind[KIND_NUM] = {
    FFTW_REDFT00,
    FFTW_REDFT10,
    FFTW_REDFT01,
    FFTW_REDFT11,
    FFTW_RODFT00,
    FFTW_RODFT10,
    FFTW_RODFT01,
    FFTW_RODFT11,
};

fftw_r2r_kind inv_kind[KIND_NUM] = {
    FFTW_REDFT00,
    FFTW_REDFT01,
    FFTW_REDFT10,
    FFTW_REDFT11,
    FFTW_RODFT00,
    FFTW_RODFT01,
    FFTW_RODFT10,
    FFTW_RODFT11,
};

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface implementation
////////////////////////////////////////////////////////////
}

IEXP_NS_END
