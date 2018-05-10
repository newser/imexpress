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

#ifndef __IEXP_WIGNER_SYMBOL__
#define __IEXP_WIGNER_SYMBOL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coupling.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// wigner 3j symbol
// ========================================

template <typename T>
class wigner3j_functor
    : public functor_m2vnum_6d<wigner3j_functor<T>, T, double>
{
  public:
    wigner3j_functor(const T &jm)
        : functor_m2vnum_6d<wigner3j_functor<T>, T, double>(jm)
    {
    }

    double m2vnum_impl(int ja, int jb, int jc, int ma, int mb, int mc) const
    {
        return gsl_sf_coupling_3j(ja, jb, jc, ma, mb, mc);
    }
};

template <typename T>
inline CwiseNullaryOp<wigner3j_functor<T>,
                      typename wigner3j_functor<T>::ResultType>
wigner3j(const DenseBase<T> &jm)
{
    using ResultType = typename wigner3j_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, jm),
                                   M2VNUM_COL(T, jm),
                                   wigner3j_functor<T>(jm.derived()));
}

template <typename T, typename U>
class wigner3j_e_functor
    : public functor_m2vnum_6d_e<wigner3j_e_functor<T, U>, T, U, double>
{
  public:
    wigner3j_e_functor(const T &jm, U &e)
        : functor_m2vnum_6d_e<wigner3j_e_functor<T, U>, T, U, double>(jm, e)
    {
    }

    double m2vnum_e_impl(
        int ja, int jb, int jc, int ma, int mb, int mc, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_coupling_3j_e(ja, jb, jc, ma, mb, mc, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner3j_e_functor<T, U>,
                      typename wigner3j_e_functor<T, U>::ResultType>
wigner3j(const DenseBase<T> &jm, DenseBase<U> &e)
{
    using ResultType = typename wigner3j_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, jm),
                                   M2VNUM_COL(T, jm),
                                   wigner3j_e_functor<T, U>(jm.derived(),
                                                            e.derived()));
}

// ========================================
// wigner 6j syjeol
// ========================================

template <typename T>
class wigner6j_functor
    : public functor_m2vnum_6d<wigner6j_functor<T>, T, double>
{
  public:
    wigner6j_functor(const T &j)
        : functor_m2vnum_6d<wigner6j_functor<T>, T, double>(j)
    {
    }

    double m2vnum_impl(int ja, int jb, int jc, int jd, int je, int jf) const
    {
        return gsl_sf_coupling_6j(ja, jb, jc, jd, je, jf);
    }
};

template <typename T>
inline CwiseNullaryOp<wigner6j_functor<T>,
                      typename wigner6j_functor<T>::ResultType>
wigner6j(const DenseBase<T> &j)
{
    using ResultType = typename wigner6j_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, j),
                                   M2VNUM_COL(T, j),
                                   wigner6j_functor<T>(j.derived()));
}

template <typename T, typename U>
class wigner6j_e_functor
    : public functor_m2vnum_6d_e<wigner6j_e_functor<T, U>, T, U, double>
{
  public:
    wigner6j_e_functor(const T &j, U &e)
        : functor_m2vnum_6d_e<wigner6j_e_functor<T, U>, T, U, double>(j, e)
    {
    }

    double m2vnum_e_impl(
        int ja, int jb, int jc, int jd, int je, int jf, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_coupling_6j_e(ja, jb, jc, jd, je, jf, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner6j_e_functor<T, U>,
                      typename wigner6j_e_functor<T, U>::ResultType>
wigner6j(const DenseBase<T> &j, DenseBase<U> &e)
{
    using ResultType = typename wigner6j_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, j),
                                   M2VNUM_COL(T, j),
                                   wigner6j_e_functor<T, U>(j.derived(),
                                                            e.derived()));
}

// ========================================
// wigner 9j syjeol
// ========================================

template <typename T>
class wigner9j_functor
    : public functor_m2vnum_9d<wigner9j_functor<T>, T, double>
{
  public:
    wigner9j_functor(const T &j)
        : functor_m2vnum_9d<wigner9j_functor<T>, T, double>(j)
    {
    }

    double m2vnum_impl(
        int ja, int jb, int jc, int jd, int je, int jf, int jg, int jh, int ji)
        const
    {
        return gsl_sf_coupling_9j(ja, jb, jc, jd, je, jf, jg, jh, ji);
    }
};

template <typename T>
inline CwiseNullaryOp<wigner9j_functor<T>,
                      typename wigner9j_functor<T>::ResultType>
wigner9j(const DenseBase<T> &j)
{
    using ResultType = typename wigner9j_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, j),
                                   M2VNUM_COL(T, j),
                                   wigner9j_functor<T>(j.derived()));
}

template <typename T, typename U>
class wigner9j_e_functor
    : public functor_m2vnum_9d_e<wigner9j_e_functor<T, U>, T, U, double>
{
  public:
    wigner9j_e_functor(const T &j, U &e)
        : functor_m2vnum_9d_e<wigner9j_e_functor<T, U>, T, U, double>(j, e)
    {
    }

    double m2vnum_e_impl(int ja,
                         int jb,
                         int jc,
                         int jd,
                         int je,
                         int jf,
                         int jg,
                         int jh,
                         int ji,
                         double &e) const
    {
        gsl_sf_result r;
        gsl_sf_coupling_9j_e(ja, jb, jc, jd, je, jf, jg, jh, ji, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner9j_e_functor<T, U>,
                      typename wigner9j_e_functor<T, U>::ResultType>
wigner9j(const DenseBase<T> &j, DenseBase<U> &e)
{
    using ResultType = typename wigner9j_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, j),
                                   M2VNUM_COL(T, j),
                                   wigner9j_e_functor<T, U>(j.derived(),
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

#endif /* __IEXP_WIGNER_SYMBOL__ */
