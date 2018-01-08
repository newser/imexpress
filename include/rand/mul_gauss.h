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

#ifndef __IEXP_RAND_MUL_GAUSS__
#define __IEXP_RAND_MUL_GAUSS__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class mgauss_rng
{
  public:
    mgauss_rng(size_t k,
               double mu[],
               double cov[],
               rng_type type = DEFAULT_RNG,
               unsigned long seed = 0)
        : m_mu_block{k, mu}
        , m_mu{k,
               1,
               mu,
               &m_mu_block,
               0}
        , m_rng(type, seed)
    {
        m_L = gsl_matrix_alloc(k, k);
        IEXP_NOT_NULLPTR(m_L);
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                gsl_matrix_set(m_L, i, j, cov[i * k + j]);
            }
        }

        gsl_linalg_cholesky_decomp1(m_L);
    }

    ~mgauss_rng()
    {
        gsl_matrix_free(m_L);
    }

    void seed(unsigned long seed)
    {
        m_rng.seed(seed);
    }

    void next(double x[])
    {
        gsl_block b{m_mu_block.size, x};
        gsl_vector result{m_mu_block.size,
                          1,
                          x,
                          &b,
                          0};
        gsl_ran_multivariate_gaussian(m_rng.gsl(), &m_mu, m_L, &result);
    }

  private:
    mgauss_rng(const mgauss_rng &) = delete;
    mgauss_rng &operator=(const mgauss_rng &other) = delete;

    gsl_matrix *m_L;
    gsl_block m_mu_block;
    gsl_vector m_mu;
    rng m_rng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_MUL_GAUSS__ */
