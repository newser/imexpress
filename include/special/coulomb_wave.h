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

#ifndef __IEXP_COULOMB_WAVE__
#define __IEXP_COULOMB_WAVE__

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
// coulomb wave function, F
// ========================================

template <typename T>
class coulw_f_functor : public functor_m2vnum_4d<coulw_f_functor<T>, T, double>
{
  public:
    coulw_f_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_f_functor<T>, T, double>(eta_x_l_k)
    {
    }

    double m2vnum_impl(double eta, double x, double l, double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return F.val;
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_f_functor<T>,
                      typename coulw_f_functor<T>::ResultType>
coulw_f(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_f_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_f_functor<T>(eta_x_l_k.derived()));
}

template <typename T, typename U>
class coulw_f_e_functor
    : public functor_m2vnum_4d_e<coulw_f_e_functor<T, U>, T, U, double>
{
  public:
    coulw_f_e_functor(const T &eta_x_l_k, U &e)
        : functor_m2vnum_4d_e<coulw_f_e_functor<T, U>, T, U, double>(eta_x_l_k,
                                                                     e)
    {
    }

    double m2vnum_e_impl(
        double eta, double x, double l, double k, double &e) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        e = F.err;
        return F.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_f_e_functor<T, U>,
                      typename coulw_f_e_functor<T, U>::ResultType>
coulw_f(const DenseBase<T> &eta_x_l_k, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_f_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_f_e_functor<T, U>(eta_x_l_k.derived(),
                                                           e.derived()));
}

// ========================================
// coulomb wave function, F derivative
// ========================================

template <typename T>
class coulw_df_functor
    : public functor_m2vnum_4d<coulw_df_functor<T>, T, double>
{
  public:
    coulw_df_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_df_functor<T>, T, double>(eta_x_l_k)
    {
    }

    double m2vnum_impl(double eta, double x, double l, double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return Fp.val;
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_df_functor<T>,
                      typename coulw_df_functor<T>::ResultType>
coulw_df(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_df_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_df_functor<T>(eta_x_l_k.derived()));
}

template <typename T, typename U>
class coulw_df_e_functor
    : public functor_m2vnum_4d_e<coulw_df_e_functor<T, U>, T, U, double>
{
  public:
    using ResultType = typename dense_derive_m2vnum<T, double>::type;

    coulw_df_e_functor(const T &eta_x_l_k, U &e)
        : functor_m2vnum_4d_e<coulw_df_e_functor<T, U>, T, U, double>(eta_x_l_k,
                                                                      e)
    {
    }

    double m2vnum_e_impl(
        double eta, double x, double l, double k, double &e) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        e = Fp.err;
        return Fp.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_df_e_functor<T, U>,
                      typename coulw_df_e_functor<T, U>::ResultType>
coulw_df(const DenseBase<T> &eta_x_l_k, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_df_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_df_e_functor<T, U>(eta_x_l_k.derived(),
                                                            e.derived()));
}

// ========================================
// coulomb wave function, G
// ========================================

template <typename T>
class coulw_g_functor : public functor_m2vnum_4d<coulw_g_functor<T>, T, double>
{
  public:
    coulw_g_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_g_functor<T>, T, double>(eta_x_l_k)
    {
    }

    double m2vnum_impl(double eta, double x, double l, double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return G.val;
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_g_functor<T>,
                      typename coulw_g_functor<T>::ResultType>
coulw_g(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_g_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_g_functor<T>(eta_x_l_k.derived()));
}

template <typename T, typename U>
class coulw_g_e_functor
    : public functor_m2vnum_4d_e<coulw_g_e_functor<T, U>, T, U, double>
{
  public:
    coulw_g_e_functor(const T &eta_x_l_k, U &e)
        : functor_m2vnum_4d_e<coulw_g_e_functor<T, U>, T, U, double>(eta_x_l_k,
                                                                     e)
    {
    }

    double m2vnum_e_impl(
        double eta, double x, double l, double k, double &e) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        e = G.err;
        return G.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_g_e_functor<T, U>,
                      typename coulw_g_e_functor<T, U>::ResultType>
coulw_g(const DenseBase<T> &eta_x_l_k, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_g_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_g_e_functor<T, U>(eta_x_l_k.derived(),
                                                           e.derived()));
}

// ========================================
// coulomb wave function, G derivative
// ========================================

template <typename T>
class coulw_dg_functor
    : public functor_m2vnum_4d<coulw_dg_functor<T>, T, double>
{
  public:
    coulw_dg_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_dg_functor<T>, T, double>(eta_x_l_k)
    {
    }

    double m2vnum_impl(double eta, double x, double l, double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return Gp.val;
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_dg_functor<T>,
                      typename coulw_dg_functor<T>::ResultType>
coulw_dg(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_dg_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_dg_functor<T>(eta_x_l_k.derived()));
}

template <typename T, typename U>
class coulw_dg_e_functor
    : public functor_m2vnum_4d_e<coulw_dg_e_functor<T, U>, T, U, double>
{
  public:
    using ResultType = typename dense_derive_m2vnum<T, double>::type;

    coulw_dg_e_functor(const T &eta_x_l_k, U &e)
        : functor_m2vnum_4d_e<coulw_dg_e_functor<T, U>, T, U, double>(eta_x_l_k,
                                                                      e)
    {
    }

    double m2vnum_e_impl(
        double eta, double x, double l, double k, double &e) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        e = Gp.err;
        return Gp.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_dg_e_functor<T, U>,
                      typename coulw_dg_e_functor<T, U>::ResultType>
coulw_dg(const DenseBase<T> &eta_x_l_k, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_dg_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_dg_e_functor<T, U>(eta_x_l_k.derived(),
                                                            e.derived()));
}

// ========================================
// coulomb wave function, F and G
// ========================================

template <typename T>
class coulw_fg_functor : public functor_m2vnum_4d<coulw_fg_functor<T>,
                                                  T,
                                                  std::tuple<double, double>>
{
  public:
    coulw_fg_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_fg_functor<T>, T, std::tuple<double, double>>(
              eta_x_l_k)
    {
    }

    std::tuple<double, double> m2vnum_impl(double eta,
                                           double x,
                                           double l,
                                           double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return std::make_tuple(F.val, G.val);
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_fg_functor<T>,
                      typename coulw_fg_functor<T>::ResultType>
coulw_fg(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_fg_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_fg_functor<T>(eta_x_l_k.derived()));
}

// ========================================
// coulomb wave function, F and G and derivative
// ========================================

template <typename T>
class coulw_dfg_functor
    : public functor_m2vnum_4d<coulw_dfg_functor<T>,
                               T,
                               std::tuple<double, double, double, double>>
{
  public:
    coulw_dfg_functor(const T &eta_x_l_k)
        : functor_m2vnum_4d<coulw_dfg_functor<T>,
                            T,
                            std::tuple<double, double, double, double>>(
              eta_x_l_k)
    {
    }

    std::tuple<double, double, double, double> m2vnum_impl(double eta,
                                                           double x,
                                                           double l,
                                                           double k) const
    {
        gsl_sf_result F, Fp, G, Gp;
        double eF, eG;
        gsl_sf_coulomb_wave_FG_e(eta, x, l, k, &F, &Fp, &G, &Gp, &eF, &eG);
        return std::make_tuple(F.val, G.val, Fp.val, Gp.val);
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_dfg_functor<T>,
                      typename coulw_dfg_functor<T>::ResultType>
coulw_dfg(const DenseBase<T> &eta_x_l_k)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_dfg_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, eta_x_l_k),
                                   M2VNUM_COL(T, eta_x_l_k),
                                   coulw_dfg_functor<T>(eta_x_l_k.derived()));
}

// ========================================
// coulomb wave function, normalization constant
// ========================================

template <typename T>
class coulw_cl_functor
    : public functor_m2vnum_2d<coulw_cl_functor<T>, T, double>
{
  public:
    coulw_cl_functor(const T &l_eta)
        : functor_m2vnum_2d<coulw_cl_functor<T>, T, double>(l_eta)
    {
    }

    double m2vnum_impl(double l, double eta) const
    {
        gsl_sf_result r;
        gsl_sf_coulomb_CL_e(l, eta, &r);
        return r.val;
    }
};

template <typename T>
inline CwiseNullaryOp<coulw_cl_functor<T>,
                      typename coulw_cl_functor<T>::ResultType>
coulw_cl(const DenseBase<T> &l_eta)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_cl_functor<T>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, l_eta),
                                   M2VNUM_COL(T, l_eta),
                                   coulw_cl_functor<T>(l_eta.derived()));
}

template <typename T, typename U>
class coulw_cl_e_functor
    : public functor_m2vnum_2d_e<coulw_cl_e_functor<T, U>, T, U, double>
{
  public:
    coulw_cl_e_functor(const T &l_eta, U &e)
        : functor_m2vnum_2d_e<coulw_cl_e_functor<T, U>, T, U, double>(l_eta, e)
    {
    }

    double m2vnum_e_impl(double l, double eta, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_coulomb_CL_e(l, eta, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_cl_e_functor<T, U>,
                      typename coulw_cl_e_functor<T, U>::ResultType>
coulw_cl(const DenseBase<T> &l_eta, DenseBase<U> &e)
{
    static_assert(TYPE_IS(typename T::Scalar, double),
                  "only support double scalar");

    using ResultType = typename coulw_cl_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, l_eta),
                                   M2VNUM_COL(T, l_eta),
                                   coulw_cl_e_functor<T, U>(l_eta.derived(),
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

#endif /* __IEXP_COULOMB_WAVE__ */
