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

#ifndef __IEXP_HIST2__
#define __IEXP_HIST2__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_histogram2d.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class hist2pdf;

class hist2
{
    friend class hist2pdf;

  public:
    hist2(const double xrange[], size_t nx, const double yrange[], size_t ny)
    {
        eigen_assert((nx > 1) && (ny > 1));
        m_gh2 = gsl_histogram2d_alloc(nx - 1, ny - 1);
        IEXP_NOT_NULLPTR(m_gh2);

        gsl_histogram2d_set_ranges(m_gh2, xrange, nx, yrange, ny);
    }

    hist2(const std::initializer_list<double> &xrange,
          const std::initializer_list<double> &yrange)
    {
        eigen_assert((xrange.size() > 1) && (yrange.size() > 1));
        m_gh2 = gsl_histogram2d_alloc(xrange.size() - 1, yrange.size() - 1);
        IEXP_NOT_NULLPTR(m_gh2);

        gsl_histogram2d_set_ranges(m_gh2,
                                   xrange.begin(),
                                   xrange.size(),
                                   yrange.begin(),
                                   yrange.size());
    }

    template <typename T, typename U>
    hist2(const DenseBase<T> &xrange, const DenseBase<U> &yrange)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "T must be double type");
        static_assert(TYPE_IS(typename U::Scalar, double),
                      "U must be double type");

        eigen_assert((xrange.size() > 1) && (yrange.size() > 1));
        m_gh2 = gsl_histogram2d_alloc(xrange.size() - 1, yrange.size() - 1);
        IEXP_NOT_NULLPTR(m_gh2);

        typename type_eval<T>::type e_xr(xrange.eval());
        typename type_eval<U>::type e_yr(yrange.eval());
        gsl_histogram2d_set_ranges(m_gh2,
                                   e_xr.data(),
                                   e_xr.size(),
                                   e_yr.data(),
                                   e_yr.size());
    }

    hist2(size_t nx,
          double xmin,
          double xmax,
          size_t ny,
          double ymin,
          double ymax)
    {
        m_gh2 = gsl_histogram2d_alloc(nx, ny);
        IEXP_NOT_NULLPTR(m_gh2);

        gsl_histogram2d_set_ranges_uniform(m_gh2, xmin, xmax, ymin, ymax);
    }

    hist2(const hist2 &h)
    {
        m_gh2 = gsl_histogram2d_clone(h.m_gh2);
        IEXP_NOT_NULLPTR(m_gh2);
    }

    hist2(hist2 &&h)
    {
        m_gh2 = h.m_gh2;
        h.m_gh2 = nullptr;
    }

    ~hist2()
    {
        if (m_gh2 != nullptr) {
            gsl_histogram2d_free(m_gh2);
        }
    }

    hist2 &operator=(const hist2 &h)
    {
        eigen_assert(gsl_histogram2d_nx(m_gh2) == gsl_histogram2d_nx(h.m_gh2));
        eigen_assert(gsl_histogram2d_ny(m_gh2) == gsl_histogram2d_ny(h.m_gh2));
        gsl_histogram2d_memcpy(m_gh2, h.m_gh2);
        return *this;
    }

    hist2 &operator=(hist2 &&h)
    {
        if (m_gh2 != nullptr) {
            gsl_histogram2d_free(m_gh2);
        }
        m_gh2 = h.m_gh2;
        h.m_gh2 = nullptr;
        return *this;
    }

    bool operator==(const hist2 &h) const
    {
        return gsl_histogram2d_equal_bins_p(m_gh2, h.m_gh2) == 1;
    }

    hist2 &operator+=(const hist2 &h)
    {
        eigen_assert(gsl_histogram2d_nx(m_gh2) == gsl_histogram2d_nx(h.m_gh2));
        eigen_assert(gsl_histogram2d_ny(m_gh2) == gsl_histogram2d_ny(h.m_gh2));
        gsl_histogram2d_add(m_gh2, h.m_gh2);
        return *this;
    }

    hist2 &operator+=(double offset)
    {
        gsl_histogram2d_shift(m_gh2, offset);
        return *this;
    }

    hist2 &operator-=(const hist2 &h)
    {
        eigen_assert(gsl_histogram2d_nx(m_gh2) == gsl_histogram2d_nx(h.m_gh2));
        eigen_assert(gsl_histogram2d_ny(m_gh2) == gsl_histogram2d_ny(h.m_gh2));
        gsl_histogram2d_sub(m_gh2, h.m_gh2);
        return *this;
    }

    hist2 &operator-=(double offset)
    {
        gsl_histogram2d_shift(m_gh2, -offset);
        return *this;
    }

    hist2 &operator*=(const hist2 &h)
    {
        eigen_assert(gsl_histogram2d_nx(m_gh2) == gsl_histogram2d_nx(h.m_gh2));
        eigen_assert(gsl_histogram2d_ny(m_gh2) == gsl_histogram2d_ny(h.m_gh2));
        gsl_histogram2d_mul(m_gh2, h.m_gh2);
        return *this;
    }

    hist2 &operator*=(double scale)
    {
        gsl_histogram2d_scale(m_gh2, scale);
        return *this;
    }

    hist2 &operator/=(const hist2 &h)
    {
        eigen_assert(gsl_histogram2d_nx(m_gh2) == gsl_histogram2d_nx(h.m_gh2));
        eigen_assert(gsl_histogram2d_ny(m_gh2) == gsl_histogram2d_ny(h.m_gh2));
        gsl_histogram2d_div(m_gh2, h.m_gh2);
        return *this;
    }

    hist2 &operator/=(double scale)
    {
        gsl_histogram2d_scale(m_gh2, 1 / scale);
        return *this;
    }

    hist2 &reset()
    {
        gsl_histogram2d_reset(m_gh2);
        return *this;
    }

    double get(size_t i, size_t j)
    {
        eigen_assert(i < gsl_histogram2d_nx(m_gh2));
        eigen_assert(j < gsl_histogram2d_ny(m_gh2));
        return gsl_histogram2d_get(m_gh2, i, j);
    }

    hist2 &add(double x, double y, double weight = 1.0)
    {
        gsl_histogram2d_accumulate(m_gh2, x, y, weight);
        return *this;
    }

    template <typename T, typename U>
    hist2 &add(const DenseBase<T> &x, const DenseBase<U> &y)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "T must be double type");
        static_assert(TYPE_IS(typename U::Scalar, double),
                      "U must be double type");

        eigen_assert((x.rows() == y.rows()) && (x.cols() == y.cols()));
        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                gsl_histogram2d_accumulate(m_gh2, x(i, j), y(i, j), 1.0);
            }
        }
        return *this;
    }

    template <typename T, typename U, typename V>
    hist2 &add(const DenseBase<T> &x,
               const DenseBase<U> &y,
               const DenseBase<V> &weight)
    {
        static_assert(TYPE_IS(typename T::Scalar, double),
                      "T must be double type");
        static_assert(TYPE_IS(typename U::Scalar, double),
                      "U must be double type");
        static_assert(TYPE_IS(typename V::Scalar, double),
                      "V must be double type");

        eigen_assert(MATRIX_SAME_SIZE(x, y));
        eigen_assert(MATRIX_SAME_SIZE(x, weight));
        for (Index i = 0; i < x.rows(); ++i) {
            for (Index j = 0; j < x.cols(); ++j) {
                gsl_histogram2d_accumulate(m_gh2,
                                           x(i, j),
                                           y(i, j),
                                           weight(i, j));
            }
        }
        return *this;
    }

    void xrange(size_t i, double &lower, double &upper) const
    {
        eigen_assert(i < gsl_histogram2d_nx(m_gh2));
        gsl_histogram2d_get_xrange(m_gh2, i, &lower, &upper);
    }

    void yrange(size_t j, double &lower, double &upper) const
    {
        eigen_assert(j < gsl_histogram2d_ny(m_gh2));
        gsl_histogram2d_get_yrange(m_gh2, j, &lower, &upper);
    }

    double xmax() const
    {
        return gsl_histogram2d_xmax(m_gh2);
    }

    double ymax() const
    {
        return gsl_histogram2d_ymax(m_gh2);
    }

    double xmin() const
    {
        return gsl_histogram2d_xmin(m_gh2);
    }

    double ymin() const
    {
        return gsl_histogram2d_ymin(m_gh2);
    }

    size_t xsize() const
    {
        return gsl_histogram2d_nx(m_gh2);
    }

    size_t ysize() const
    {
        return gsl_histogram2d_ny(m_gh2);
    }

    bool find(double x, double y, size_t &i, size_t &j) const
    {
        try {
            gsl_histogram2d_find(m_gh2, x, y, &i, &j);
            return true;
        } catch (...) {
            return false;
        }
    }

    double max_val() const
    {
        return gsl_histogram2d_max_val(m_gh2);
    }

    void max_idx(size_t &i, size_t &j) const
    {
        gsl_histogram2d_max_bin(m_gh2, &i, &j);
    }

    double min_val() const
    {
        return gsl_histogram2d_min_val(m_gh2);
    }

    void min_idx(size_t &i, size_t &j) const
    {
        gsl_histogram2d_min_bin(m_gh2, &i, &j);
    }

    double xmean() const
    {
        return gsl_histogram2d_xmean(m_gh2);
    }

    double ymean() const
    {
        return gsl_histogram2d_ymean(m_gh2);
    }

    double xstd() const
    {
        return gsl_histogram2d_xsigma(m_gh2);
    }

    double ystd() const
    {
        return gsl_histogram2d_ysigma(m_gh2);
    }

    double sum() const
    {
        return gsl_histogram2d_sum(m_gh2);
    }

    double cov() const
    {
        return gsl_histogram2d_cov(m_gh2);
    }

  private:
    gsl_histogram2d *m_gh2;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_HIST2__ */
