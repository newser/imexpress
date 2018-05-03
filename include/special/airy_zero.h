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

#ifndef __IEXP_SPECIAL_AIRY_ZERO__
#define __IEXP_SPECIAL_AIRY_ZERO__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

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
// location of n-th 0 of airy Ai
// ========================================

template <typename T>
class airy_n0_Ai_functor
    : public functor_foreach<airy_n0_Ai_functor<T>, T, double>
{
  public:
    airy_n0_Ai_functor(const T &x)
        : functor_foreach<airy_n0_Ai_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(unsigned int x) const
    {
        return gsl_sf_airy_zero_Ai(x);
    }
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Ai_functor<T>,
                      typename airy_n0_Ai_functor<T>::ResultType>
airy_n0_Ai(const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Ai_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Ai_functor<T>(x.derived()));
}

template <typename T, typename U>
class airy_n0_Ai_e_functor
    : public functor_foreach_e<airy_n0_Ai_e_functor<T, U>, T, U, double>
{
  public:
    airy_n0_Ai_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_n0_Ai_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(unsigned int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Ai_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Ai_e_functor<T, U>,
                      typename airy_n0_Ai_e_functor<T, U>::ResultType>
airy_n0_Ai(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Ai_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Ai_e_functor<T, U>(x.derived(),
                                                              e.derived()));
}

// ========================================
// location of n-th 0 of airy Ai derivative
// ========================================

template <typename T>
class airy_n0_Ai_deriv_functor
    : public functor_foreach<airy_n0_Ai_deriv_functor<T>, T, double>
{
  public:
    airy_n0_Ai_deriv_functor(const T &x)
        : functor_foreach<airy_n0_Ai_deriv_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(unsigned int x) const
    {
        return gsl_sf_airy_zero_Ai_deriv(x);
    }
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Ai_deriv_functor<T>,
                      typename airy_n0_Ai_deriv_functor<T>::ResultType>
airy_n0_Ai_deriv(const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Ai_deriv_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Ai_deriv_functor<T>(x.derived()));
}

template <typename T, typename U>
class airy_n0_Ai_deriv_e_functor
    : public functor_foreach_e<airy_n0_Ai_deriv_e_functor<T, U>, T, U, double>
{
  public:
    airy_n0_Ai_deriv_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_n0_Ai_deriv_e_functor<T, U>, T, U, double>(x,
                                                                            e)
    {
    }

    double foreach_e_impl(unsigned int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Ai_deriv_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Ai_deriv_e_functor<T, U>,
                      typename airy_n0_Ai_deriv_e_functor<T, U>::ResultType>
airy_n0_Ai_deriv(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Ai_deriv_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Ai_deriv_e_functor<T,
                                                              U>(x.derived(),
                                                                 e.derived()));
}

// ========================================
// location of n-th 0 of airy Bi
// ========================================

template <typename T>
class airy_n0_Bi_functor
    : public functor_foreach<airy_n0_Bi_functor<T>, T, double>
{
  public:
    airy_n0_Bi_functor(const T &x)
        : functor_foreach<airy_n0_Bi_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(unsigned int x) const
    {
        return gsl_sf_airy_zero_Bi(x);
    }
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Bi_functor<T>,
                      typename airy_n0_Bi_functor<T>::ResultType>
airy_n0_Bi(const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Bi_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Bi_functor<T>(x.derived()));
}

template <typename T, typename U>
class airy_n0_Bi_e_functor
    : public functor_foreach_e<airy_n0_Bi_e_functor<T, U>, T, U, double>
{
  public:
    airy_n0_Bi_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_n0_Bi_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(unsigned int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Bi_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Bi_e_functor<T, U>,
                      typename airy_n0_Bi_e_functor<T, U>::ResultType>
airy_n0_Bi(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Bi_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Bi_e_functor<T, U>(x.derived(),
                                                              e.derived()));
}

// ========================================
// location of n-th 0 of airy Bi derivative
// ========================================

template <typename T>
class airy_n0_Bi_deriv_functor
    : public functor_foreach<airy_n0_Bi_deriv_functor<T>, T, double>
{
  public:
    airy_n0_Bi_deriv_functor(const T &x)
        : functor_foreach<airy_n0_Bi_deriv_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(unsigned int x) const
    {
        return gsl_sf_airy_zero_Bi_deriv(x);
    }
};

template <typename T>
inline CwiseNullaryOp<airy_n0_Bi_deriv_functor<T>,
                      typename airy_n0_Bi_deriv_functor<T>::ResultType>
airy_n0_Bi_deriv(const DenseBase<T> &x)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Bi_deriv_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Bi_deriv_functor<T>(x.derived()));
}

template <typename T, typename U>
class airy_n0_Bi_deriv_e_functor
    : public functor_foreach_e<airy_n0_Bi_deriv_e_functor<T, U>, T, U, double>
{
  public:
    airy_n0_Bi_deriv_e_functor(const T &x, U &e)
        : functor_foreach_e<airy_n0_Bi_deriv_e_functor<T, U>, T, U, double>(x,
                                                                            e)
    {
    }

    double foreach_e_impl(unsigned int x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_airy_zero_Bi_deriv_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<airy_n0_Bi_deriv_e_functor<T, U>,
                      typename airy_n0_Bi_deriv_e_functor<T, U>::ResultType>
airy_n0_Bi_deriv(const DenseBase<T> &x, DenseBase<U> &e)
{
    static_assert(IS_INTEGER(typename T::Scalar),
                  "only support integer scalar");

    using ResultType = typename airy_n0_Bi_deriv_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   airy_n0_Bi_deriv_e_functor<T,
                                                              U>(x.derived(),
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

#endif /* __IEXP_SPECIAL_AIRY_ZERO__ */
