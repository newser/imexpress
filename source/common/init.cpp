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

#include <common/init.h>

#include <gsl/gsl_errno.h>

#include <string>

IEXP_NS_BEGIN

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

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

static void iexp_gsl_error_handler(const char *reason,
                                   const char *file,
                                   int line,
                                   int gsl_errno);

////////////////////////////////////////////////////////////
// interface implementation
////////////////////////////////////////////////////////////

bool init()
{
    gsl_set_error_handler(iexp_gsl_error_handler);

    return true;
}

void exit()
{
}

void iexp_gsl_error_handler(const char *reason,
                            const char *file,
                            int line,
                            int gsl_errno)
{
    std::string es;
    es.reserve(256);
    es.append(gsl_strerror(gsl_errno));
    es.append(":");
    es.append(reason);
    es.append(". see[");
    es.append(file);
    es.append(":");
    es.append(std::to_string(line));
    es.append("]");

    throw std::runtime_error(es);
}

IEXP_NS_END
