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
        gsl_rstat_reset(m_rw);
    }

    void add(double x)
    {
        gsl_rstat_add(x, m_rw);
    }

    size_t size() const
    {
        return gsl_rstat_n(m_rw);
    }

    double min() const
    {
        return gsl_rstat_min(m_rw);
    }

    double max() const
    {
        return gsl_rstat_max(m_rw);
    }

    double mean() const
    {
        return gsl_rstat_mean(m_rw);
    }

    double var() const
    {
        return gsl_rstat_variance(m_rw);
    }

    double std() const
    {
        return gsl_rstat_sd(m_rw);
    }

    double std_mean() const
    {
        return gsl_rstat_sd_mean(m_rw);
    }

    double rms() const
    {
        return gsl_rstat_rms(m_rw);
    }

    double skewness() const
    {
        return gsl_rstat_skew(m_rw);
    }

    double kurtosis() const
    {
        return gsl_rstat_kurtosis(m_rw);
    }

    double median() const
    {
        return gsl_rstat_median(m_rw);
    }

  private:
    rstat(const rstat &) = delete;
    rstat &operator=(const rstat &) = delete;

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
