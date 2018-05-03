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

#ifndef __IEXP_SPECIAL_AIRY__
#define __IEXP_SPECIAL_AIRY__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <common/def.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_airy.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// airy Ai
// ========================================

template <bool scaled, precision p, typename T>
class airy_Ai_functor
    : public functor_foreach<airy_Ai_functor<scaled, p, T>, T, double>
{
  public:
    airy_Ai_functor(const T &x)
        : functor_foreach<airy_Ai_functor<scaled, p, T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return foreach_impl(x, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_impl(double x, std::true_type) const
    {
        return gsl_sf_airy_Ai_scaled(x, (int)p);
    }

    double foreach_impl(double x, std::false_type) const
    {
        return gsl_sf_airy_Ai(x, (int)p);
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void>
inline CwiseNullaryOp<airy_Ai_functor<scaled, p, T>,
                      typename airy_Ai_functor<scaled, p, T>::ResultType>
airy_Ai(const DenseBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Ai_functor<scaled, p, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_Ai_functor<scaled, p, T>(x.derived()));
}

template <bool scaled, precision p, typename T, typename U>
class airy_Ai_e_functor
    : public functor_foreach_e<airy_Ai_e_functor<scaled, p, T, U>, T, U, double>
{
  public:
    airy_Ai_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_Ai_e_functor<scaled, p, T, U>, T, U, double>(x,
                                                                              e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return foreach_e_impl(x, e, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_e_impl(double x, double &e, std::true_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_scaled_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(double x, double &e, std::false_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<airy_Ai_e_functor<scaled, p, T, U>,
                      typename airy_Ai_e_functor<scaled, p, T, U>::ResultType>
airy_Ai(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Ai_e_functor<scaled, p, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Ai_e_functor<scaled, p, T, U>(x.derived(),
                                                       e.derived()));
}

// ========================================
// airy Ai derivative
// ========================================

template <bool scaled, precision p, typename T>
class airy_Ai_deriv_functor
    : public functor_foreach<airy_Ai_deriv_functor<scaled, p, T>, T, double>
{
  public:
    airy_Ai_deriv_functor(const T &x)
        : functor_foreach<airy_Ai_deriv_functor<scaled, p, T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return foreach_impl(x, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_impl(double x, std::true_type) const
    {
        return gsl_sf_airy_Ai_deriv_scaled(x, (int)p);
    }

    double foreach_impl(double x, std::false_type) const
    {
        return gsl_sf_airy_Ai_deriv(x, (int)p);
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void>
inline CwiseNullaryOp<airy_Ai_deriv_functor<scaled, p, T>,
                      typename airy_Ai_deriv_functor<scaled, p, T>::ResultType>
airy_Ai_deriv(const DenseBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Ai_deriv_functor<scaled, p, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_Ai_deriv_functor<scaled, p, T>(
                                       x.derived()));
}

template <bool scaled, precision p, typename T, typename U>
class airy_Ai_deriv_e_functor
    : public functor_foreach_e<airy_Ai_deriv_e_functor<scaled, p, T, U>,
                               T,
                               U,
                               double>
{
  public:
    airy_Ai_deriv_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_Ai_deriv_e_functor<scaled, p, T, U>,
                            T,
                            U,
                            double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return foreach_e_impl(x, e, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_e_impl(double x, double &e, std::true_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_deriv_scaled_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(double x, double &e, std::false_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Ai_deriv_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    airy_Ai_deriv_e_functor<scaled, p, T, U>,
    typename airy_Ai_deriv_e_functor<scaled, p, T, U>::ResultType>
airy_Ai_deriv(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType =
        typename airy_Ai_deriv_e_functor<scaled, p, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Ai_deriv_e_functor<scaled, p, T, U>(x.derived(),
                                                             e.derived()));
}

// ========================================
// airy Bi
// ========================================

template <bool scaled, precision p, typename T>
class airy_Bi_functor
    : public functor_foreach<airy_Bi_functor<scaled, p, T>, T, double>
{
  public:
    airy_Bi_functor(const T &x)
        : functor_foreach<airy_Bi_functor<scaled, p, T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return foreach_impl(x, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_impl(double x, std::true_type) const
    {
        return gsl_sf_airy_Bi_scaled(x, (int)p);
    }

    double foreach_impl(double x, std::false_type) const
    {
        return gsl_sf_airy_Bi(x, (int)p);
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void>
inline CwiseNullaryOp<airy_Bi_functor<scaled, p, T>,
                      typename airy_Bi_functor<scaled, p, T>::ResultType>
airy_Bi(const DenseBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Bi_functor<scaled, p, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_Bi_functor<scaled, p, T>(x.derived()));
}

template <bool scaled, precision p, typename T, typename U>
class airy_Bi_e_functor
    : public functor_foreach_e<airy_Bi_e_functor<scaled, p, T, U>, T, U, double>
{
  public:
    airy_Bi_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_Bi_e_functor<scaled, p, T, U>, T, U, double>(x,
                                                                              e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return foreach_e_impl(x, e, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_e_impl(double x, double &e, std::true_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_scaled_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(double x, double &e, std::false_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<airy_Bi_e_functor<scaled, p, T, U>,
                      typename airy_Bi_e_functor<scaled, p, T, U>::ResultType>
airy_Bi(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Bi_e_functor<scaled, p, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Bi_e_functor<scaled, p, T, U>(x.derived(),
                                                       e.derived()));
}

// ========================================
// airy Bi derivative
// ========================================

template <bool scaled, precision p, typename T>
class airy_Bi_deriv_functor
    : public functor_foreach<airy_Bi_deriv_functor<scaled, p, T>, T, double>
{
  public:
    airy_Bi_deriv_functor(const T &x)
        : functor_foreach<airy_Bi_deriv_functor<scaled, p, T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return foreach_impl(x, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_impl(double x, std::true_type) const
    {
        return gsl_sf_airy_Bi_deriv_scaled(x, (int)p);
    }

    double foreach_impl(double x, std::false_type) const
    {
        return gsl_sf_airy_Bi_deriv(x, (int)p);
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void>
inline CwiseNullaryOp<airy_Bi_deriv_functor<scaled, p, T>,
                      typename airy_Bi_deriv_functor<scaled, p, T>::ResultType>
airy_Bi_deriv(const DenseBase<T> &x)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename airy_Bi_deriv_functor<scaled, p, T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_Bi_deriv_functor<scaled, p, T>(
                                       x.derived()));
}

template <bool scaled, precision p, typename T, typename U>
class airy_Bi_deriv_e_functor
    : public functor_foreach_e<airy_Bi_deriv_e_functor<scaled, p, T, U>,
                               T,
                               U,
                               double>
{
  public:
    airy_Bi_deriv_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_Bi_deriv_e_functor<scaled, p, T, U>,
                            T,
                            U,
                            double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        return foreach_e_impl(x, e, TYPE_BOOL(scaled)());
    }

  private:
    double foreach_e_impl(double x, double &e, std::true_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_deriv_scaled_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }

    double foreach_e_impl(double x, double &e, std::false_type) const
    {
        gsl_sf_result r;
        gsl_sf_airy_Bi_deriv_e(x, (int)p, &r);
        e = r.err;
        return r.val;
    }
};

template <bool scaled = false,
          precision p = precision::DOUBLE,
          typename T = void,
          typename U = void>
inline CwiseNullaryOp<
    airy_Bi_deriv_e_functor<scaled, p, T, U>,
    typename airy_Bi_deriv_e_functor<scaled, p, T, U>::ResultType>
airy_Bi_deriv(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType =
        typename airy_Bi_deriv_e_functor<scaled, p, T, U>::ResultType;
    return ResultType::
        NullaryExpr(x.rows(),
                    x.cols(),
                    airy_Bi_deriv_e_functor<scaled, p, T, U>(x.derived(),
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

#endif /* __IEXP_SPECIAL_AIRY__ */
