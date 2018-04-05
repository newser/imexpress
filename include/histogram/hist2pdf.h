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

#ifndef __IEXP_HIST2PDF__
#define __IEXP_HIST2PDF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <histogram/hist2.h>
#include <rand/rng.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class hist2pdf
{
  public:
    hist2pdf(const hist2 &h,
             rand::rng_type type = rand::DEFAULT_RNG,
             unsigned long seed = 0)
        : m_ghp(nullptr)
        , m_xrng(type, seed)
        , m_yrng(type, seed)
    {
        m_ghp = gsl_histogram2d_pdf_alloc(h.xsize(), h.ysize());
        IEXP_NOT_NULLPTR(m_ghp);

        if (gsl_histogram2d_pdf_init(m_ghp, h.m_gh2) != GSL_SUCCESS) {
            RETURN_OR_THROW(std::runtime_error("hist2pdf"));
        }
    }

    ~hist2pdf()
    {
        gsl_histogram2d_pdf_free(m_ghp);
    }

    void next(double &x, double &y)
    {
        gsl_histogram2d_pdf_sample(m_ghp,
                                   m_xrng.uniform_double(),
                                   m_yrng.uniform_double(),
                                   &x,
                                   &y);
    }

    template <typename T, typename U>
    void next(DenseBase<T> &x, DenseBase<U> &y)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "T must be double type");
        static_assert(TYPE_IS(typename U::Scalar, double),
                      "U must be double type");

        eigen_assert(MATRIX_SAME_SIZE(x, y));
        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                next(x(i, j), y(i, j));
            }
        }
    }

  private:
    hist2pdf(const hist2pdf &h) = delete;
    hist2pdf(hist2pdf &&h) = delete;

    hist2pdf &operator=(const hist2pdf &rs) = delete;
    hist2pdf &operator=(const hist2pdf &&c) = delete;

    gsl_histogram2d_pdf *m_ghp;
    rand::rng m_xrng, m_yrng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_HIST2PDF__ */
