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

#ifndef __IEXP_RSTAT__
#define __IEXP_RSTAT__

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

class rstat
{
  public:
    rstat()
    {
        m_rw = gsl_rstat_alloc();
        IEXP_NOT_NULLPTR(m_rw);
    }

    ~rstat()
    {
        gsl_rstat_free(m_rw);
    }

    void reset()
    {
        if (gsl_rstat_reset(m_rw) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("rstat reset"));
        }
    }

    void add(double x) const
    {
        if (gsl_rstat_add(x, m_rw) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("rstat add"));
        }
    }

    size_t size()
    {
        return gsl_rstat_n(m_rw);
    }

    double min()
    {
        return gsl_rstat_min(m_rw);
    }

    double max()
    {
        return gsl_rstat_max(m_rw);
    }

    double mean()
    {
        return gsl_rstat_mean(m_rw);
    }

    double var()
    {
        return gsl_rstat_variance(m_rw);
    }

    double std()
    {
        return gsl_rstat_sd(m_rw);
    }

    double std_mean()
    {
        return gsl_rstat_sd_mean(m_rw);
    }

    double rms()
    {
        return gsl_rstat_rms(m_rw);
    }

    double skewness()
    {
        return gsl_rstat_skew(m_rw);
    }

    double kurtosis()
    {
        return gsl_rstat_kurtosis(m_rw);
    }

    double median()
    {
        return gsl_rstat_median(m_rw);
    }

  private:
    rstat(const rstat &rs) = delete;
    rstat(rstat &&rs) = delete;

    rstat &operator=(const rstat &rs) = delete;
    rstat &operator=(const rstat &&c) = delete;

    gsl_rstat_workspace *m_rw;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_RSTAT__ */
