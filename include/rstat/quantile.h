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

#ifndef __IEXP_QUANTILE__
#define __IEXP_QUANTILE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>
#include <common/util.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rstat.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class quantile
{
  public:
    quantile(double f)
    {
        m_qw = gsl_rstat_quantile_alloc(f);
        IEXP_NOT_NULLPTR(m_qw);
    }

    ~quantile()
    {
        gsl_rstat_quantile_free(m_qw);
    }

    void reset()
    {
        if (gsl_rstat_quantile_reset(m_qw) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("quantile reset"));
        }
    }

    void add(double x) const
    {
        if (gsl_rstat_quantile_add(x, m_qw) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("quantile add"));
        }
    }

    double get()
    {
        return gsl_rstat_quantile_get(m_qw);
    }

  private:
    quantile(const quantile &rs) = delete;
    quantile(quantile &&rs) = delete;

    quantile &operator=(const quantile &rs) = delete;
    quantile &operator=(const quantile &&c) = delete;

    gsl_rstat_quantile_workspace *m_qw;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_QUANTILE__ */
