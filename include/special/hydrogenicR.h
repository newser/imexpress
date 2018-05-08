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

#ifndef __IEXP_HYDROGENIC_R__
#define __IEXP_HYDROGENIC_R__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coulomb.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// hydrogenicR_1
// ========================================

template <typename T>
class hydroR1_functor : public functor_m2vnum_2d<hydroR1_functor<T>, T, double>
{
  public:
    hydroR1_functor(const T &z_r)
        : functor_m2vnum_2d<hydroR1_functor<T>, T, double>(z_r)
    {
    }

    double m2vnum_impl(double z, double r) const
    {
        return gsl_sf_hydrogenicR_1(z, r);
    }
};

template <typename T>
inline CwiseNullaryOp<hydroR1_functor<T>,
                      typename hydroR1_functor<T>::ResultType>
hydroR1(const DenseBase<T> &z_r)
{
    using ResultType = typename hydroR1_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, z_r),
                                   M2VNUM_COL(T, z_r),
                                   hydroR1_functor<T>(z_r.derived()));
}

template <typename T, typename U>
class hydroR1_e_functor
    : public functor_m2vnum_2d_e<hydroR1_e_functor<T, U>, T, U, double>
{
  public:
    hydroR1_e_functor(const T &z_r, U &e)
        : functor_m2vnum_2d_e<hydroR1_e_functor<T, U>, T, U, double>(z_r, e)
    {
    }

    double m2vnum_e_impl(double z, double r, double &e) const
    {
        gsl_sf_result ret;
        gsl_sf_hydrogenicR_1_e(z, r, &ret);
        e = ret.err;
        return ret.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<hydroR1_e_functor<T, U>,
                      typename hydroR1_e_functor<T, U>::ResultType>
hydroR1(const DenseBase<T> &z_r, DenseBase<U> &e)
{
    using ResultType = typename hydroR1_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, z_r),
                                   M2VNUM_COL(T, z_r),
                                   hydroR1_e_functor<T, U>(z_r.derived(),
                                                           e.derived()));
}

// ========================================
// hydrogenicR
// ========================================

template <typename T>
class hydroR_functor : public functor_m2vnum_4d<hydroR_functor<T>, T, double>
{
  public:
    hydroR_functor(const T &n_l_z_r)
        : functor_m2vnum_4d<hydroR_functor<T>, T, double>(n_l_z_r)
    {
    }

    double m2vnum_impl(double n, double l, double z, double r) const
    {
        // if converting between double and int becomes perf bottleneck, we can
        // implement another api which can specify n and l
        return gsl_sf_hydrogenicR((int)n, (int)l, z, r);
    }
};

template <typename T>
inline CwiseNullaryOp<hydroR_functor<T>, typename hydroR_functor<T>::ResultType>
hydroR(const DenseBase<T> &n_l_z_r)
{
    using ResultType = typename hydroR_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, n_l_z_r),
                                   M2VNUM_COL(T, n_l_z_r),
                                   hydroR_functor<T>(n_l_z_r.derived()));
}

template <typename T, typename U>
class hydroR_e_functor
    : public functor_m2vnum_4d_e<hydroR_e_functor<T, U>, T, U, double>
{
  public:
    hydroR_e_functor(const T &n_l_z_r, U &e)
        : functor_m2vnum_4d_e<hydroR_e_functor<T, U>, T, U, double>(n_l_z_r, e)
    {
    }

    double m2vnum_e_impl(
        double n, double l, double z, double r, double &e) const
    {
        gsl_sf_result ret;
        gsl_sf_hydrogenicR_e(n, l, z, r, &ret);
        e = ret.err;
        return ret.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<hydroR_e_functor<T, U>,
                      typename hydroR_e_functor<T, U>::ResultType>
hydroR(const DenseBase<T> &n_l_z_r, DenseBase<U> &e)
{
    using ResultType = typename hydroR_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, n_l_z_r),
                                   M2VNUM_COL(T, n_l_z_r),
                                   hydroR_e_functor<T, U>(n_l_z_r.derived(),
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

#endif /* __IEXP_HYDROGENIC_R__ */
