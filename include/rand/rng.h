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

#ifndef __IEXP_RNG__
#define __IEXP_RNG__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_rng.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define DEFAULT_RNG_TYPE iexp::rand::rng::type::MT19937

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class rng
{
  public:
    enum class type
    {
        BOROSH13,
        COVEYOU,
        CMRG,
        FISHMAN18,
        FISHMAN20,
        FISHMAN2X,
        GFSR4,
        KNUTHRAN,
        KNUTHRAN2,
        KNUTHRAN2002,
        LECUYER21,
        MINSTD,
        MRG,
        MT19937,
        MT19937_1999,
        MT19937_1998,
        R250,
        RAN0,
        RAN1,
        RAN2,
        RAN3,
        RAND,
        RAND48,
        RANDOM128_BSD,
        RANDOM128_GLIBC2,
        RANDOM128_LIBC5,
        RANDOM256_BSD,
        RANDOM256_GLIBC2,
        RANDOM256_LIBC5,
        RANDOM32_BSD,
        RANDOM32_GLIBC2,
        RANDOM32_LIBC5,
        RANDOM64_BSD,
        RANDOM64_GLIBC2,
        RANDOM64_LIBC5,
        RANDOM8_BSD,
        RANDOM8_GLIBC2,
        RANDOM8_LIBC5,
        RANDOM_BSD,
        RANDOM_GLIBC2,
        RANDOM_LIBC5,
        RANDU,
        RANF,
        RANLUX,
        RANLUX389,
        RANLXD1,
        RANLXD2,
        RANLXS0,
        RANLXS1,
        RANLXS2,
        RANMAR,
        SLATEC,
        TAUS,
        TAUS2,
        TAUS113,
        TRANSPUTER,
        TT800,
        UNI,
        UNI32,
        VAX,
        WATERMAN14,
        ZUF,
    };

  public:
    rng(type type = DEFAULT_RNG_TYPE, unsigned long seed = 0);

    ~rng()
    {
        if (m_rng != nullptr) {
            gsl_rng_free(m_rng);
        }
    }

    rng(const rng &rhs)
    {
        m_rng = gsl_rng_clone(rhs.m_rng);
        eigen_assert(m_rng != nullptr);
    }

    rng(rng &&rhs)
    {
        m_rng = rhs.m_rng;
        rhs.m_rng = nullptr;
    }

    rng &operator=(const rng &rhs)
    {
        gsl_rng_memcpy(m_rng, rhs.m_rng);
        return *this;
    }

    rng &operator=(rng &&rhs)
    {
        if (m_rng != nullptr) {
            gsl_rng_free(m_rng);
        }
        m_rng = rhs.m_rng;
        rhs.m_rng = nullptr;
        return *this;
    }

    rng &seed(unsigned long seed)
    {
        gsl_rng_set(m_rng, seed);
        return *this;
    }

    unsigned long uniform_ulong() const
    {
        return gsl_rng_get(m_rng);
    }

    // [0, n - 1]
    unsigned long uniform_ulong(unsigned long n) const
    {
        n = std::min(n, max());
        return gsl_rng_uniform_int(m_rng, n);
    }

    // [0, 1)
    double uniform_double() const
    {
        return gsl_rng_uniform(m_rng);
    }

    // (0, 1)
    double uniform_pos_double() const
    {
        return gsl_rng_uniform_pos(m_rng);
    }

    const char *name() const
    {
        return gsl_rng_name(m_rng);
    }

    unsigned long min() const
    {
        return gsl_rng_min(m_rng);
    }

    unsigned long max() const
    {
        return gsl_rng_max(m_rng);
    }

    gsl_rng *gsl() const
    {
        return m_rng;
    }

  private:
    gsl_rng *m_rng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RNG__ */
