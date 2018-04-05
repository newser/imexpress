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

#ifndef __IEXP_HISTPDF__
#define __IEXP_HISTPDF__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <histogram/hist.h>
#include <rand/rng.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class histpdf
{
  public:
    histpdf(const hist &h,
            rand::rng_type type = rand::DEFAULT_RNG,
            unsigned long seed = 0)
        : m_ghp(nullptr)
        , m_rng(type, seed)
    {
        m_ghp = gsl_histogram_pdf_alloc(h.size());
        IEXP_NOT_NULLPTR(m_ghp);

        gsl_histogram_pdf_init(m_ghp, h.m_gh);
    }

    ~histpdf()
    {
        gsl_histogram_pdf_free(m_ghp);
    }

    double next()
    {
        return gsl_histogram_pdf_sample(m_ghp, m_rng.uniform_double());
    }

    template <typename T>
    void next(DenseBase<T> &x)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "must be double type");
        // static_assert(IS_LREF(decltype(x(0, 0))), "not lvalue");

        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                x(i, j) = next();
            }
        }
    }

  private:
    histpdf(const histpdf &h) = delete;
    histpdf(histpdf &&h) = delete;

    histpdf &operator=(const histpdf &rs) = delete;
    histpdf &operator=(const histpdf &&c) = delete;

    gsl_histogram_pdf *m_ghp;
    rand::rng m_rng;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_HISTPDF__ */
