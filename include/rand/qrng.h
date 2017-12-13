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

#ifndef __IEXP_QRNG__
#define __IEXP_QRNG__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_qrng.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

enum qrng_type
{
    NIEDERREITER_2,
    SOBOL,
    HALTON,
    REVERSEHALTON,

    QRNG_TYPE_NUM
};

class qrng
{
  public:
    qrng(qrng_type type, unsigned int dim);

    ~qrng()
    {
        if (m_qrng != nullptr) {
            gsl_qrng_free(m_qrng);
        }
    }

    qrng(const qrng &other)
    {
        m_qrng = gsl_qrng_clone(other.m_qrng);
        eigen_assert(m_qrng != nullptr);
    }

    qrng(qrng &&other)
    {
        m_qrng = other.m_qrng;
        other.m_qrng = nullptr;
    }

    qrng &operator=(const qrng &other)
    {
        gsl_qrng_memcpy(m_qrng, other.m_qrng);
        return *this;
    }

    qrng &operator=(qrng &&other)
    {
        if (m_qrng != nullptr) {
            gsl_qrng_free(m_qrng);
        }
        m_qrng = other.m_qrng;
        other.m_qrng = nullptr;
        return *this;
    }

    void reset()
    {
        gsl_qrng_init(m_qrng);
    }

    void next(double x[]) const
    {
        gsl_qrng_get(m_qrng, x);
    }

    const char *name() const
    {
        return gsl_qrng_name(m_qrng);
    }

  private:
    gsl_qrng *m_qrng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_QRNG__ */
