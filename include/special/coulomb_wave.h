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
inline T coulw_f_impl(const T l, const T eta, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_f_impl(const double l, const double eta, const double x)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        return F.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_f"));
}

template <typename T>
class coulw_f_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_f_functor(const T &l, const T &eta, const T &x)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_f_impl(m_l(i, j), m_eta(i, j), m_x(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<coulw_f_functor<T>,
                      typename coulw_f_functor<T>::ArrayType>
coulw_f(const ArrayBase<T> &l, const ArrayBase<T> &eta, const ArrayBase<T> &x)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    typedef typename coulw_f_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_f_functor<T>(l.derived(),
                                                     eta.derived(),
                                                     x.derived()));
}

template <typename T>
inline T coulw_f_e_impl(const T l, const T eta, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_f_e_impl(const double l,
                             const double eta,
                             const double x,
                             double &e)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        e = F.err;
        return F.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_f"));
}

template <typename T, typename U>
class coulw_f_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_f_e_functor(const T &l, const T &eta, const T &x, U &e)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_f_e_impl(m_l(i, j), m_eta(i, j), m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_f_e_functor<T, U>,
                      typename coulw_f_e_functor<T, U>::ArrayType>
coulw_f(const ArrayBase<T> &l,
        const ArrayBase<T> &eta,
        const ArrayBase<T> &x,
        ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));
    eigen_assert(MATRIX_SAME_SIZE(l, e));

    typedef typename coulw_f_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_f_e_functor<T, U>(l.derived(),
                                                          eta.derived(),
                                                          x.derived(),
                                                          e.derived()));
}

// ========================================
// coulomb wave function, F derivative
// ========================================

template <typename T>
inline T coulw_df_impl(const T l, const T eta, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_df_impl(const double l, const double eta, const double x)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        return Fp.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_df"));
}

template <typename T>
class coulw_df_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_df_functor(const T &l, const T &eta, const T &x)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_df_impl(m_l(i, j), m_eta(i, j), m_x(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<coulw_df_functor<T>,
                      typename coulw_df_functor<T>::ArrayType>
coulw_df(const ArrayBase<T> &l, const ArrayBase<T> &eta, const ArrayBase<T> &x)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    typedef typename coulw_df_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_df_functor<T>(l.derived(),
                                                      eta.derived(),
                                                      x.derived()));
}

template <typename T>
inline T coulw_df_e_impl(const T l, const T eta, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_df_e_impl(const double l,
                              const double eta,
                              const double x,
                              double &e)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        e = Fp.err;
        return Fp.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_df"));
}

template <typename T, typename U>
class coulw_df_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_df_e_functor(const T &l, const T &eta, const T &x, U &e)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_df_e_impl(m_l(i, j), m_eta(i, j), m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_df_e_functor<T, U>,
                      typename coulw_df_e_functor<T, U>::ArrayType>
coulw_df(const ArrayBase<T> &l,
         const ArrayBase<T> &eta,
         const ArrayBase<T> &x,
         ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));
    eigen_assert(MATRIX_SAME_SIZE(l, e));

    typedef typename coulw_df_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_df_e_functor<T, U>(l.derived(),
                                                           eta.derived(),
                                                           x.derived(),
                                                           e.derived()));
}

// ========================================
// coulomb wave function, G
// ========================================

template <typename T>
inline T coulw_g_impl(const T l, const T eta, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_g_impl(const double l, const double eta, const double x)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        return G.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_g"));
}

template <typename T>
class coulw_g_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_g_functor(const T &l, const T &eta, const T &x)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_g_impl(m_l(i, j), m_eta(i, j), m_x(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<coulw_g_functor<T>,
                      typename coulw_g_functor<T>::ArrayType>
coulw_g(const ArrayBase<T> &l, const ArrayBase<T> &eta, const ArrayBase<T> &x)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    typedef typename coulw_g_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_g_functor<T>(l.derived(),
                                                     eta.derived(),
                                                     x.derived()));
}

template <typename T>
inline T coulw_g_e_impl(const T l, const T eta, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_g_e_impl(const double l,
                             const double eta,
                             const double x,
                             double &e)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        e = G.err;
        return G.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_g"));
}

template <typename T, typename U>
class coulw_g_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_g_e_functor(const T &l, const T &eta, const T &x, U &e)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_g_e_impl(m_l(i, j), m_eta(i, j), m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_g_e_functor<T, U>,
                      typename coulw_g_e_functor<T, U>::ArrayType>
coulw_g(const ArrayBase<T> &l,
        const ArrayBase<T> &eta,
        const ArrayBase<T> &x,
        ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));
    eigen_assert(MATRIX_SAME_SIZE(l, e));

    typedef typename coulw_g_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_g_e_functor<T, U>(l.derived(),
                                                          eta.derived(),
                                                          x.derived(),
                                                          e.derived()));
}

// ========================================
// coulomb wave function, G derivative
// ========================================

template <typename T>
inline T coulw_dg_impl(const T l, const T eta, const T x)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_dg_impl(const double l, const double eta, const double x)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        return Gp.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_dg"));
}

template <typename T>
class coulw_dg_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_dg_functor(const T &l, const T &eta, const T &x)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_dg_impl(m_l(i, j), m_eta(i, j), m_x(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
};

template <typename T>
inline CwiseNullaryOp<coulw_dg_functor<T>,
                      typename coulw_dg_functor<T>::ArrayType>
coulw_dg(const ArrayBase<T> &l, const ArrayBase<T> &eta, const ArrayBase<T> &x)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    typedef typename coulw_dg_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_dg_functor<T>(l.derived(),
                                                      eta.derived(),
                                                      x.derived()));
}

template <typename T>
inline T coulw_dg_e_impl(const T l, const T eta, const T x, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_dg_e_impl(const double l,
                              const double eta,
                              const double x,
                              double &e)
{
    gsl_sf_result F, Fp, G, Gp;
    double eF, eG;
    int s = gsl_sf_coulomb_wave_FG_e(eta, x, l, 0, &F, &Fp, &G, &Gp, &eF, &eG);
    if (s == GSL_SUCCESS) {
        e = Gp.err;
        return Gp.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_dg"));
}

template <typename T, typename U>
class coulw_dg_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_dg_e_functor(const T &l, const T &eta, const T &x, U &e)
        : m_l(l)
        , m_eta(eta)
        , m_x(x)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_dg_e_impl(m_l(i, j), m_eta(i, j), m_x(i, j), m_e(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    const T &m_x;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_dg_e_functor<T, U>,
                      typename coulw_dg_e_functor<T, U>::ArrayType>
coulw_dg(const ArrayBase<T> &l,
         const ArrayBase<T> &eta,
         const ArrayBase<T> &x,
         ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));
    eigen_assert(MATRIX_SAME_SIZE(l, e));

    typedef typename coulw_dg_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_dg_e_functor<T, U>(l.derived(),
                                                           eta.derived(),
                                                           x.derived(),
                                                           e.derived()));
}

// ========================================
// coulomb wave function, F and G
// ========================================

template <typename T>
inline void coulw_fg_impl(
    const T l, const int k, const T eta, const T x, T &F, T &G)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void coulw_fg_impl(const double l,
                          const int k,
                          const double eta,
                          const double x,
                          double &F,
                          double &G)
{
    gsl_sf_result rF, rFp, rG, rGp;
    double eF, eG;
    if (gsl_sf_coulomb_wave_FG_e(eta,
                                 x,
                                 l,
                                 k,
                                 &rF,
                                 &rFp,
                                 &rG,
                                 &rGp,
                                 &eF,
                                 &eG) == GSL_SUCCESS) {
        F = rF.val;
        G = rG.val;
    }
}

template <typename T, typename U>
inline void coulw_fg(const ArrayBase<T> &l,
                     const ArrayBase<U> &k,
                     const ArrayBase<T> &eta,
                     const ArrayBase<T> &x,
                     ArrayBase<T> &F,
                     ArrayBase<T> &G)
{
    static_assert(TYPE_IS(typename U::Scalar, int) ||
                      TYPE_IS(typename U::Scalar, unsigned int),
                  "k can only be int or uint array");

    eigen_assert(MATRIX_SAME_SIZE(l, k));
    eigen_assert(MATRIX_SAME_SIZE(l, k));
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    // F.derived().resize(l.rows(), l.cols());
    // G.derived().resize(l.rows(), l.cols());

    // loop unroll?
    // care row/col major?
    for (Index i = 0; i < l.cols(); ++i) {
        for (Index j = 0; j < l.rows(); ++j) {
            coulw_fg_impl(l(j, i),
                          k(j, i),
                          eta(j, i),
                          x(j, i),
                          F.derived()(j, i),
                          G.derived()(j, i));
        }
    }
}

// ========================================
// coulomb wave function, F and G and derivative
// ========================================

template <typename T>
inline void coulw_dfg_impl(
    const T l, const int k, const T eta, const T x, T &F, T &G, T &dF, T &dG)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void coulw_dfg_impl(const double l,
                           const int k,
                           const double eta,
                           const double x,
                           double &F,
                           double &G,
                           double &dF,
                           double &dG)
{
    gsl_sf_result rF, rFp, rG, rGp;
    double eF, eG;
    if (gsl_sf_coulomb_wave_FG_e(eta,
                                 x,
                                 l,
                                 k,
                                 &rF,
                                 &rFp,
                                 &rG,
                                 &rGp,
                                 &eF,
                                 &eG) == GSL_SUCCESS) {
        F = rF.val;
        G = rG.val;
        dF = rFp.val;
        dG = rGp.val;
    }
}

template <typename T, typename U>
inline void coulw_dfg(const ArrayBase<T> &l,
                      const ArrayBase<U> &k,
                      const ArrayBase<T> &eta,
                      const ArrayBase<T> &x,
                      ArrayBase<T> &F,
                      ArrayBase<T> &G,
                      ArrayBase<T> &dF,
                      ArrayBase<T> &dG)
{
    static_assert(TYPE_IS(typename U::Scalar, int) ||
                      TYPE_IS(typename U::Scalar, unsigned int),
                  "k can only be int or uint array");

    eigen_assert(MATRIX_SAME_SIZE(l, k));
    eigen_assert(MATRIX_SAME_SIZE(l, k));
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, x));

    // F.derived().resize(l.rows(), l.cols());
    // G.derived().resize(l.rows(), l.cols());
    // dF.derived().resize(l.rows(), l.cols());
    // dG.derived().resize(l.rows(), l.cols());

    // loop unroll?
    // care row/col major?
    for (Index i = 0; i < l.cols(); ++i) {
        for (Index j = 0; j < l.rows(); ++j) {
            coulw_dfg_impl(l(j, i),
                           k(j, i),
                           eta(j, i),
                           x(j, i),
                           F.derived()(j, i),
                           G.derived()(j, i),
                           dF.derived()(j, i),
                           dG.derived()(j, i));
        }
    }
}

// ========================================
// coulomb wave function, normalization constant
// ========================================

template <typename T>
inline T coulw_cl_impl(const T l, const T eta)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_cl_impl(const double l, const double eta)
{
    gsl_sf_result r;
    if (gsl_sf_coulomb_CL_e(l, eta, &r) == GSL_SUCCESS) {
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_cl"));
}

template <typename T>
class coulw_cl_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_cl_functor(const T &l, const T &eta)
        : m_l(l)
        , m_eta(eta)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_cl_impl(m_l(i, j), m_eta(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
};

template <typename T>
inline CwiseNullaryOp<coulw_cl_functor<T>,
                      typename coulw_cl_functor<T>::ArrayType>
coulw_cl(const ArrayBase<T> &l, const ArrayBase<T> &eta)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));

    typedef typename coulw_cl_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_cl_functor<T>(l.derived(),
                                                      eta.derived()));
}

template <typename T>
inline T coulw_cl_e_impl(const T l, const T eta, T &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double coulw_cl_e_impl(const double l, const double eta, double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coulomb_CL_e(l, eta, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    RETURN_NAN_OR_THROW(std::runtime_error("coulw_cl"));
}

template <typename T, typename U>
class coulw_cl_e_functor
{
  public:
    typedef Array<typename T::Scalar,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    coulw_cl_e_functor(const T &l, const T &eta, U &e)
        : m_l(l)
        , m_eta(eta)
        , m_e(e)
    {
    }

    const typename T::Scalar operator()(Index i, Index j) const
    {
        return coulw_cl_e_impl(m_l(i, j), m_eta(i, j), m_e(i, j));
    }

  private:
    const T &m_l;
    const T &m_eta;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<coulw_cl_e_functor<T, U>,
                      typename coulw_cl_e_functor<T, U>::ArrayType>
coulw_cl(const ArrayBase<T> &l, const ArrayBase<T> &eta, ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(l, eta));
    eigen_assert(MATRIX_SAME_SIZE(l, e));

    typedef typename coulw_cl_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(l.rows(),
                                  l.cols(),
                                  coulw_cl_e_functor<T, U>(l.derived(),
                                                           eta.derived(),
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
