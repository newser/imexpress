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

#ifndef __IEXP_HIST__
#define __IEXP_HIST__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_histogram.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class histpdf;

class hist
{
    friend class histpdf;

  public:
    hist(const double range[], size_t n)
    {
        eigen_assert(n > 1);
        m_gh = gsl_histogram_alloc(n - 1);
        IEXP_NOT_NULLPTR(m_gh);

        gsl_histogram_set_ranges(m_gh, range, n);
    }

    hist(const std::initializer_list<double> &range)
    {
        eigen_assert(range.size() > 1);
        m_gh = gsl_histogram_alloc(range.size() - 1);
        IEXP_NOT_NULLPTR(m_gh);

        gsl_histogram_set_ranges(m_gh, range.begin(), range.size());
    }

    template <typename T>
    hist(const DenseBase<T> &range)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "must be double type");

        eigen_assert(range.size() > 1);
        m_gh = gsl_histogram_alloc(range.size() - 1);
        IEXP_NOT_NULLPTR(m_gh);

        typename type_eval<T>::type e_r(range.eval());
        gsl_histogram_set_ranges(m_gh, e_r.data(), e_r.size());
    }

    hist(size_t n, double min, double max)
    {
        m_gh = gsl_histogram_alloc(n);
        IEXP_NOT_NULLPTR(m_gh);

        gsl_histogram_set_ranges_uniform(m_gh, min, max);
    }

    hist(const hist &h)
    {
        m_gh = gsl_histogram_clone(h.m_gh);
        IEXP_NOT_NULLPTR(m_gh);
    }

    hist(hist &&h)
    {
        m_gh = h.m_gh;
        h.m_gh = nullptr;
    }

    ~hist()
    {
        if (m_gh != nullptr) {
            gsl_histogram_free(m_gh);
        }
    }

    hist &operator=(const hist &h)
    {
        eigen_assert(gsl_histogram_bins(m_gh) == gsl_histogram_bins(h.m_gh));
        gsl_histogram_memcpy(m_gh, h.m_gh);
        return *this;
    }

    hist &operator=(hist &&h)
    {
        if (m_gh != nullptr) {
            gsl_histogram_free(m_gh);
        }
        m_gh = h.m_gh;
        h.m_gh = nullptr;
        return *this;
    }

    double operator[](size_t i)
    {
        eigen_assert(i < gsl_histogram_bins(m_gh));
        return gsl_histogram_get(m_gh, i);
    }

    bool operator==(const hist &h) const
    {
        return gsl_histogram_equal_bins_p(m_gh, h.m_gh) == 1;
    }

    hist &operator<<(double x)
    {
        gsl_histogram_increment(m_gh, x);
        return *this;
    }

    template <typename T>
    hist &operator<<(const DenseBase<T> &x)
    {
        add(x);
        return *this;
    }

    hist &operator+=(const hist &h)
    {
        eigen_assert(gsl_histogram_bins(m_gh) == gsl_histogram_bins(h.m_gh));
        gsl_histogram_add(m_gh, h.m_gh);
        return *this;
    }

    hist &operator+=(double offset)
    {
        gsl_histogram_shift(m_gh, offset);
        return *this;
    }

    hist &operator-=(const hist &h)
    {
        eigen_assert(gsl_histogram_bins(m_gh) == gsl_histogram_bins(h.m_gh));
        gsl_histogram_sub(m_gh, h.m_gh);
        return *this;
    }

    hist &operator-=(double offset)
    {
        gsl_histogram_shift(m_gh, -offset);
        return *this;
    }

    hist &operator*=(const hist &h)
    {
        eigen_assert(gsl_histogram_bins(m_gh) == gsl_histogram_bins(h.m_gh));
        gsl_histogram_mul(m_gh, h.m_gh);
        return *this;
    }

    hist &operator*=(double scale)
    {
        gsl_histogram_scale(m_gh, scale);
        return *this;
    }

    hist &operator/=(const hist &h)
    {
        eigen_assert(gsl_histogram_bins(m_gh) == gsl_histogram_bins(h.m_gh));
        gsl_histogram_div(m_gh, h.m_gh);
        return *this;
    }

    hist &operator/=(double scale)
    {
        gsl_histogram_scale(m_gh, 1 / scale);
        return *this;
    }

    void reset()
    {
        gsl_histogram_reset(m_gh);
    }

    void add(double x, double weight = 1.0)
    {
        gsl_histogram_accumulate(m_gh, x, weight);
    }

    template <typename T>
    void add(const DenseBase<T> &x)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "must be double type");

        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                gsl_histogram_accumulate(m_gh, x(i, j), 1.0);
            }
        }
    }

    template <typename T, typename U>
    void add(const DenseBase<T> &x, const DenseBase<U> &weight)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "T must be double type");
        static_assert(TYPE_IS(typename U::Scalar, double),
                      "U must be double type");

        eigen_assert((x.rows() == weight.rows()) &&
                     (x.cols() == weight.cols()));
        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                gsl_histogram_accumulate(m_gh, x(i, j), weight(i, j));
            }
        }
    }

    void range(size_t i, double &lower, double &upper) const
    {
        eigen_assert(i < gsl_histogram_bins(m_gh));
        gsl_histogram_get_range(m_gh, i, &lower, &upper);
    }

    double max() const
    {
        return gsl_histogram_max(m_gh);
    }

    double min() const
    {
        return gsl_histogram_min(m_gh);
    }

    size_t size() const
    {
        return gsl_histogram_bins(m_gh);
    }

    size_t find(double x) const
    {
        try {
            size_t i;
            gsl_histogram_find(m_gh, x, &i);
            return i;
        } catch (...) {
            return (size_t)-1;
        }
    }

    double max_val() const
    {
        return gsl_histogram_max_val(m_gh);
    }

    size_t max_idx() const
    {
        return gsl_histogram_max_bin(m_gh);
    }

    double min_val() const
    {
        return gsl_histogram_min_val(m_gh);
    }

    size_t min_idx() const
    {
        return gsl_histogram_min_bin(m_gh);
    }

    double mean() const
    {
        return gsl_histogram_mean(m_gh);
    }

    double std() const
    {
        return gsl_histogram_sigma(m_gh);
    }

    double sum() const
    {
        return gsl_histogram_sum(m_gh);
    }

  private:
    gsl_histogram *m_gh;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_HIST__ */
