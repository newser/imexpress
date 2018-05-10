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

#ifndef __IEXP_SF_DILOG__
#define __IEXP_SF_DILOG__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_dilog.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class dilog_functor
    : public functor_foreach<dilog_functor<T>, T, typename T::Scalar>
{
  public:
    dilog_functor(const T &x)
        : functor_foreach<dilog_functor<T>, T, typename T::Scalar>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_dilog(x);
    }

    std::complex<double> foreach_impl(std::complex<double> x) const
    {
        // compute norm and arg here may introduce error and impact perf. can
        // implement api with norm and arg as its parameters
        double norm = std::norm(x);
        double arg = std::arg(x);
        gsl_sf_result re, im;
        gsl_sf_complex_dilog_e(norm, arg, &re, &im);
        return std::complex<double>(re.val, im.val);
    }
};

template <typename T>
inline CwiseNullaryOp<dilog_functor<T>, typename dilog_functor<T>::ResultType>
dilog(const DenseBase<T> &x)
{
    using ResultType = typename dilog_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dilog_functor<T>(x.derived()));
}

template <typename T, typename U>
class dilog_e_functor
    : public functor_foreach_e<dilog_e_functor<T, U>, T, U, double>
{
  public:
    dilog_e_functor(const T &x, U &e)
        : functor_foreach_e<dilog_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_dilog_e(x, &r);
        e = r.err;
        return r.val;
    }

    std::complex<double> foreach_e_impl(std::complex<double> x,
                                        std::complex<double> &e) const
    {
        double norm = std::norm(x);
        double arg = std::arg(x);
        gsl_sf_result re, im;
        gsl_sf_complex_dilog_e(norm, arg, &re, &im);
        e.real(re.err);
        e.imag(im.err);
        return std::complex<double>(re.val, im.val);
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<dilog_e_functor<T, U>,
                      typename dilog_e_functor<T, U>::ResultType>
dilog(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename dilog_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   dilog_e_functor<T, U>(x.derived(),
                                                         e.derived()));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_SF_DILOG__ */
