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

#ifndef __IEXP_ELLIPTIC_INTEGRAL__
#define __IEXP_ELLIPTIC_INTEGRAL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_ellint.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// legendre first kind, complete
// ========================================

template <typename T, precision p>
class ellint_k_functor
    : public functor_foreach<ellint_k_functor<T, p>, T, double>
{
  public:
    ellint_k_functor(const T &k)
        : functor_foreach<ellint_k_functor<T, p>, T, double>(k)
    {
    }

    double foreach_impl(double k) const
    {
        return gsl_sf_ellint_Kcomp(k, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_k_functor<T, p>,
                      typename ellint_k_functor<T, p>::ResultType>
ellint_k(const DenseBase<T> &k)
{
    using ResultType = typename ellint_k_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(k.rows(),
                                   k.cols(),
                                   ellint_k_functor<T, p>(k.derived()));
}

template <typename T, typename U, precision p>
class ellint_k_e_functor
    : public functor_foreach_e<ellint_k_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_k_e_functor(const T &k, U &e)
        : functor_foreach_e<ellint_k_e_functor<T, U, p>, T, U, double>(k, e)
    {
    }

    double foreach_e_impl(double k, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_Kcomp_e(k, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_k_e_functor<T, U, p>,
                      typename ellint_k_e_functor<T, U, p>::ResultType>
ellint_k(const DenseBase<T> &k, DenseBase<U> &e)
{
    using ResultType = typename ellint_k_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(k.rows(),
                                   k.cols(),
                                   ellint_k_e_functor<T, U, p>(k.derived(),
                                                               e.derived()));
}

// ========================================
// legendre first kind, incomplete
// ========================================

template <typename T, precision p>
class ellint_ki_functor
    : public functor_m2vnum_2d<ellint_ki_functor<T, p>, T, double>
{
  public:
    ellint_ki_functor(const T &phi_k)
        : functor_m2vnum_2d<ellint_ki_functor<T, p>, T, double>(phi_k)
    {
    }

    double m2vnum_impl(double phi, double k) const
    {
        return gsl_sf_ellint_F(phi, k, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_ki_functor<T, p>,
                      typename ellint_ki_functor<T, p>::ResultType>
ellint_ki(const DenseBase<T> &phi_k)
{
    using ResultType = typename ellint_ki_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k),
                                   M2VNUM_COL(T, phi_k),
                                   ellint_ki_functor<T, p>(phi_k.derived()));
}

template <typename T, typename U, precision p>
class ellint_ki_e_functor
    : public functor_m2vnum_2d_e<ellint_ki_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_ki_e_functor(const T &phi_k, U &e)
        : functor_m2vnum_2d_e<ellint_ki_e_functor<T, U, p>, T, U, double>(phi_k,
                                                                          e)
    {
    }

    double m2vnum_e_impl(double phi, double k, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_F_e(phi, k, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_ki_e_functor<T, U, p>,
                      typename ellint_ki_e_functor<T, U, p>::ResultType>
ellint_ki(const DenseBase<T> &phi_k, DenseBase<U> &e)
{
    using ResultType = typename ellint_ki_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k),
                                   M2VNUM_COL(T, phi_k),
                                   ellint_ki_e_functor<T, U, p>(phi_k.derived(),
                                                                e.derived()));
}

// ========================================
// legendre second kind, complete
// ========================================

template <typename T, precision p>
class ellint_e_functor
    : public functor_foreach<ellint_e_functor<T, p>, T, double>
{
  public:
    ellint_e_functor(const T &k)
        : functor_foreach<ellint_e_functor<T, p>, T, double>(k)
    {
    }

    double foreach_impl(double k) const
    {
        return gsl_sf_ellint_Ecomp(k, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_e_functor<T, p>,
                      typename ellint_e_functor<T, p>::ResultType>
ellint_e(const DenseBase<T> &k)
{
    using ResultType = typename ellint_e_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(k.rows(),
                                   k.cols(),
                                   ellint_e_functor<T, p>(k.derived()));
}

template <typename T, typename U, precision p>
class ellint_e_e_functor
    : public functor_foreach_e<ellint_e_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_e_e_functor(const T &k, U &e)
        : functor_foreach_e<ellint_e_e_functor<T, U, p>, T, U, double>(k, e)
    {
    }

    double foreach_e_impl(double k, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_Ecomp_e(k, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_e_e_functor<T, U, p>,
                      typename ellint_e_e_functor<T, U, p>::ResultType>
ellint_e(const DenseBase<T> &k, DenseBase<U> &e)
{
    using ResultType = typename ellint_e_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(k.rows(),
                                   k.cols(),
                                   ellint_e_e_functor<T, U, p>(k.derived(),
                                                               e.derived()));
}

// ========================================
// legendre second kind, incomplete
// ========================================

template <typename T, precision p>
class ellint_ei_functor
    : public functor_m2vnum_2d<ellint_ei_functor<T, p>, T, double>
{
  public:
    ellint_ei_functor(const T &phi_k)
        : functor_m2vnum_2d<ellint_ei_functor<T, p>, T, double>(phi_k)
    {
    }

    double m2vnum_impl(double phi, double k) const
    {
        return gsl_sf_ellint_E(phi, k, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_ei_functor<T, p>,
                      typename ellint_ei_functor<T, p>::ResultType>
ellint_ei(const DenseBase<T> &phi_k)
{
    using ResultType = typename ellint_ei_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k),
                                   M2VNUM_COL(T, phi_k),
                                   ellint_ei_functor<T, p>(phi_k.derived()));
}

template <typename T, typename U, precision p>
class ellint_ei_e_functor
    : public functor_m2vnum_2d_e<ellint_ei_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_ei_e_functor(const T &phi_k, U &e)
        : functor_m2vnum_2d_e<ellint_ei_e_functor<T, U, p>, T, U, double>(phi_k,
                                                                          e)
    {
    }

    double m2vnum_e_impl(double phi, double k, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_E_e(phi, k, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_ei_e_functor<T, U, p>,
                      typename ellint_ei_e_functor<T, U, p>::ResultType>
ellint_ei(const DenseBase<T> &phi_k, DenseBase<U> &e)
{
    using ResultType = typename ellint_ei_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k),
                                   M2VNUM_COL(T, phi_k),
                                   ellint_ei_e_functor<T, U, p>(phi_k.derived(),
                                                                e.derived()));
}

// ========================================
// legendre third kind, complete
// ========================================

template <typename T, precision p>
class ellint_p_functor
    : public functor_m2vnum_2d<ellint_p_functor<T, p>, T, double>
{
  public:
    ellint_p_functor(const T &k_n)
        : functor_m2vnum_2d<ellint_p_functor<T, p>, T, double>(k_n)
    {
    }

    double m2vnum_impl(double k, double n) const
    {
        return gsl_sf_ellint_Pcomp(k, n, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_p_functor<T, p>,
                      typename ellint_p_functor<T, p>::ResultType>
ellint_p(const DenseBase<T> &k_n)
{
    using ResultType = typename ellint_p_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, k_n),
                                   M2VNUM_COL(T, k_n),
                                   ellint_p_functor<T, p>(k_n.derived()));
}

template <typename T, typename U, precision p>
class ellint_p_e_functor
    : public functor_m2vnum_2d_e<ellint_p_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_p_e_functor(const T &k_n, U &e)
        : functor_m2vnum_2d_e<ellint_p_e_functor<T, U, p>, T, U, double>(k_n, e)
    {
    }

    double m2vnum_e_impl(double k, double n, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_Pcomp_e(k, n, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_p_e_functor<T, U, p>,
                      typename ellint_p_e_functor<T, U, p>::ResultType>
ellint_p(const DenseBase<T> &k_n, DenseBase<U> &e)
{
    using ResultType = typename ellint_p_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, k_n),
                                   M2VNUM_COL(T, k_n),
                                   ellint_p_e_functor<T, U, p>(k_n.derived(),
                                                               e.derived()));
}

// ========================================
// legendre third kind, incomplete
// ========================================

template <typename T, precision p>
class ellint_pi_functor
    : public functor_m2vnum_3d<ellint_pi_functor<T, p>, T, double>
{
  public:
    ellint_pi_functor(const T &phi_k_n)
        : functor_m2vnum_3d<ellint_pi_functor<T, p>, T, double>(phi_k_n)
    {
    }

    double m2vnum_impl(double phi, double k, double n) const
    {
        return gsl_sf_ellint_P(phi, k, n, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_pi_functor<T, p>,
                      typename ellint_pi_functor<T, p>::ResultType>
ellint_pi(const DenseBase<T> &phi_k_n)
{
    using ResultType = typename ellint_pi_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k_n),
                                   M2VNUM_COL(T, phi_k_n),
                                   ellint_pi_functor<T, p>(phi_k_n.derived()));
}

template <typename T, typename U, precision p>
class ellint_pi_e_functor
    : public functor_m2vnum_3d_e<ellint_pi_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_pi_e_functor(const T &phi_k_n, U &e)
        : functor_m2vnum_3d_e<ellint_pi_e_functor<T, U, p>,
                              T,
                              U,
                              double>(phi_k_n, e)
    {
    }

    double m2vnum_e_impl(double phi, double k, double n, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_P_e(phi, k, n, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_pi_e_functor<T, U, p>,
                      typename ellint_pi_e_functor<T, U, p>::ResultType>
ellint_pi(const DenseBase<T> &phi_k_n, DenseBase<U> &e)
{
    using ResultType = typename ellint_pi_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, phi_k_n),
                                   M2VNUM_COL(T, phi_k_n),
                                   ellint_pi_e_functor<T, U, p>(phi_k_n
                                                                    .derived(),
                                                                e.derived()));
}

// ========================================
// carlson RC
// ========================================

template <typename T, precision p>
class ellint_rc_functor
    : public functor_m2vnum_2d<ellint_rc_functor<T, p>, T, double>
{
  public:
    ellint_rc_functor(const T &x_y)
        : functor_m2vnum_2d<ellint_rc_functor<T, p>, T, double>(x_y)
    {
    }

    double m2vnum_impl(double x, double y) const
    {
        return gsl_sf_ellint_RC(x, y, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_rc_functor<T, p>,
                      typename ellint_rc_functor<T, p>::ResultType>
ellint_rc(const DenseBase<T> &x_y)
{
    using ResultType = typename ellint_rc_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y),
                                   M2VNUM_COL(T, x_y),
                                   ellint_rc_functor<T, p>(x_y.derived()));
}

template <typename T, typename U, precision p>
class ellint_rc_e_functor
    : public functor_m2vnum_2d_e<ellint_rc_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_rc_e_functor(const T &x_y, U &e)
        : functor_m2vnum_2d_e<ellint_rc_e_functor<T, U, p>, T, U, double>(x_y,
                                                                          e)
    {
    }

    double m2vnum_e_impl(double x, double y, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_RC_e(x, y, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_rc_e_functor<T, U, p>,
                      typename ellint_rc_e_functor<T, U, p>::ResultType>
ellint_rc(const DenseBase<T> &x_y, DenseBase<U> &e)
{
    using ResultType = typename ellint_rc_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y),
                                   M2VNUM_COL(T, x_y),
                                   ellint_rc_e_functor<T, U, p>(x_y.derived(),
                                                                e.derived()));
}

// ========================================
// carlson RD
// ========================================

template <typename T, precision p>
class ellint_rd_functor
    : public functor_m2vnum_3d<ellint_rd_functor<T, p>, T, double>
{
  public:
    ellint_rd_functor(const T &x_y_z)
        : functor_m2vnum_3d<ellint_rd_functor<T, p>, T, double>(x_y_z)
    {
    }

    double m2vnum_impl(double x, double y, double z) const
    {
        return gsl_sf_ellint_RD(x, y, z, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_rd_functor<T, p>,
                      typename ellint_rd_functor<T, p>::ResultType>
ellint_rd(const DenseBase<T> &x_y_z)
{
    using ResultType = typename ellint_rd_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z),
                                   M2VNUM_COL(T, x_y_z),
                                   ellint_rd_functor<T, p>(x_y_z.derived()));
}

template <typename T, typename U, precision p>
class ellint_rd_e_functor
    : public functor_m2vnum_3d_e<ellint_rd_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_rd_e_functor(const T &x_y_z, U &e)
        : functor_m2vnum_3d_e<ellint_rd_e_functor<T, U, p>, T, U, double>(x_y_z,
                                                                          e)
    {
    }

    double m2vnum_e_impl(double x, double y, double z, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_RD_e(x, y, z, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_rd_e_functor<T, U, p>,
                      typename ellint_rd_e_functor<T, U, p>::ResultType>
ellint_rd(const DenseBase<T> &x_y_z, DenseBase<U> &e)
{
    using ResultType = typename ellint_rd_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z),
                                   M2VNUM_COL(T, x_y_z),
                                   ellint_rd_e_functor<T, U, p>(x_y_z.derived(),
                                                                e.derived()));
}

// ========================================
// carlson RF
// ========================================

template <typename T, precision p>
class ellint_rf_functor
    : public functor_m2vnum_3d<ellint_rf_functor<T, p>, T, double>
{
  public:
    ellint_rf_functor(const T &x_y_z)
        : functor_m2vnum_3d<ellint_rf_functor<T, p>, T, double>(x_y_z)
    {
    }

    double m2vnum_impl(double x, double y, double z) const
    {
        return gsl_sf_ellint_RF(x, y, z, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_rf_functor<T, p>,
                      typename ellint_rf_functor<T, p>::ResultType>
ellint_rf(const DenseBase<T> &x_y_z)
{
    using ResultType = typename ellint_rf_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z),
                                   M2VNUM_COL(T, x_y_z),
                                   ellint_rf_functor<T, p>(x_y_z.derived()));
}

template <typename T, typename U, precision p>
class ellint_rf_e_functor
    : public functor_m2vnum_3d_e<ellint_rf_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_rf_e_functor(const T &x_y_z, U &e)
        : functor_m2vnum_3d_e<ellint_rf_e_functor<T, U, p>, T, U, double>(x_y_z,
                                                                          e)
    {
    }

    double m2vnum_e_impl(double x, double y, double z, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_RF_e(x, y, z, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_rf_e_functor<T, U, p>,
                      typename ellint_rf_e_functor<T, U, p>::ResultType>
ellint_rf(const DenseBase<T> &x_y_z, DenseBase<U> &e)
{
    using ResultType = typename ellint_rf_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z),
                                   M2VNUM_COL(T, x_y_z),
                                   ellint_rf_e_functor<T, U, p>(x_y_z.derived(),
                                                                e.derived()));
}

// ========================================
// carlson RJ
// ========================================

template <typename T, precision p>
class ellint_rj_functor
    : public functor_m2vnum_4d<ellint_rj_functor<T, p>, T, double>
{
  public:
    ellint_rj_functor(const T &x_y_z_p)
        : functor_m2vnum_4d<ellint_rj_functor<T, p>, T, double>(x_y_z_p)
    {
    }

    double m2vnum_impl(double x, double y, double z, double u) const
    {
        return gsl_sf_ellint_RJ(x, y, z, u, (gsl_mode_t)p);
    }
};

template <precision p = precision::DOUBLE, typename T = void>
inline CwiseNullaryOp<ellint_rj_functor<T, p>,
                      typename ellint_rj_functor<T, p>::ResultType>
ellint_rj(const DenseBase<T> &x_y_z_p)
{
    using ResultType = typename ellint_rj_functor<T, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z_p),
                                   M2VNUM_COL(T, x_y_z_p),
                                   ellint_rj_functor<T, p>(x_y_z_p.derived()));
}

template <typename T, typename U, precision p>
class ellint_rj_e_functor
    : public functor_m2vnum_4d_e<ellint_rj_e_functor<T, U, p>, T, U, double>
{
  public:
    ellint_rj_e_functor(const T &x_y_z_p, U &e)
        : functor_m2vnum_4d_e<ellint_rj_e_functor<T, U, p>,
                              T,
                              U,
                              double>(x_y_z_p, e)
    {
    }

    double m2vnum_e_impl(
        double x, double y, double z, double u, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_ellint_RJ_e(x, y, z, u, (gsl_mode_t)p, &r);
        e = r.err;
        return r.val;
    }
};

template <precision p = precision::DOUBLE, typename T = void, typename U = void>
inline CwiseNullaryOp<ellint_rj_e_functor<T, U, p>,
                      typename ellint_rj_e_functor<T, U, p>::ResultType>
ellint_rj(const DenseBase<T> &x_y_z_p, DenseBase<U> &e)
{
    using ResultType = typename ellint_rj_e_functor<T, U, p>::ResultType;
    return ResultType::NullaryExpr(M2VNUM_ROW(T, x_y_z_p),
                                   M2VNUM_COL(T, x_y_z_p),
                                   ellint_rj_e_functor<T, U, p>(x_y_z_p
                                                                    .derived(),
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

#endif /* __IEXP_ELLIPTIC_INTEGRAL__ */
