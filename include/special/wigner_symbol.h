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
inline double wigner3j_impl(
    const T ja, const T jb, const T jc, const T ma, const T mb, const T mc)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner3j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int ma,
                            const int mb,
                            const int mc)
{
    return gsl_sf_coupling_3j(ja, jb, jc, ma, mb, mc);
}

template <typename T>
class wigner3j_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner3j_functor(const T &ja,
                     const T &jb,
                     const T &jc,
                     const T &ma,
                     const T &mb,
                     const T &mc)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_ma(ma)
        , m_mb(mb)
        , m_mc(mc)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner3j_impl(m_ja(i, j),
                             m_jb(i, j),
                             m_jc(i, j),
                             m_ma(i, j),
                             m_mb(i, j),
                             m_mc(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_ma;
    const T &m_mb;
    const T &m_mc;
};

template <typename T>
inline CwiseNullaryOp<wigner3j_functor<T>,
                      typename wigner3j_functor<T>::ArrayType>
wigner3j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &ma,
         const ArrayBase<T> &mb,
         const ArrayBase<T> &mc)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, ma));
    eigen_assert(MATRIX_SAME_SIZE(ja, mb));
    eigen_assert(MATRIX_SAME_SIZE(ja, mc));

    typedef typename wigner3j_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner3j_functor<T>(ja.derived(),
                                                      jb.derived(),
                                                      jc.derived(),
                                                      ma.derived(),
                                                      mb.derived(),
                                                      mc.derived()));
}

template <typename T, typename U>
inline double wigner3j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T ma,
                              const T mb,
                              const T mc,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner3j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int ma,
                              const int mb,
                              const int mc,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_3j_e(ja, jb, jc, ma, mb, mc, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("wigner3j"));
}

template <typename T, typename U>
class wigner3j_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner3j_e_functor(const T &ja,
                       const T &jb,
                       const T &jc,
                       const T &ma,
                       const T &mb,
                       const T &mc,
                       U &e)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_ma(ma)
        , m_mb(mb)
        , m_mc(mc)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner3j_e_impl(m_ja(i, j),
                               m_jb(i, j),
                               m_jc(i, j),
                               m_ma(i, j),
                               m_mb(i, j),
                               m_mc(i, j),
                               m_e(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_ma;
    const T &m_mb;
    const T &m_mc;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner3j_e_functor<T, U>,
                      typename wigner3j_e_functor<T, U>::ArrayType>
wigner3j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &ma,
         const ArrayBase<T> &mb,
         const ArrayBase<T> &mc,
         ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, ma));
    eigen_assert(MATRIX_SAME_SIZE(ja, mb));
    eigen_assert(MATRIX_SAME_SIZE(ja, mc));

    typedef typename wigner3j_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner3j_e_functor<T, U>(ja.derived(),
                                                           jb.derived(),
                                                           jc.derived(),
                                                           ma.derived(),
                                                           mb.derived(),
                                                           mc.derived(),
                                                           e.derived()));
}

// ========================================
// wigner 6j syjeol
// ========================================

template <typename T>
inline double wigner6j_impl(
    const T ja, const T jb, const T jc, const T jd, const T je, const T jf)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner6j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int jd,
                            const int je,
                            const int jf)
{
    return gsl_sf_coupling_6j(ja, jb, jc, jd, je, jf);
}

template <typename T>
class wigner6j_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner6j_functor(const T &ja,
                     const T &jb,
                     const T &jc,
                     const T &jd,
                     const T &je,
                     const T &jf)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_jd(jd)
        , m_je(je)
        , m_jf(jf)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner6j_impl(m_ja(i, j),
                             m_jb(i, j),
                             m_jc(i, j),
                             m_jd(i, j),
                             m_je(i, j),
                             m_jf(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_jd;
    const T &m_je;
    const T &m_jf;
};

template <typename T>
inline CwiseNullaryOp<wigner6j_functor<T>,
                      typename wigner6j_functor<T>::ArrayType>
wigner6j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &jd,
         const ArrayBase<T> &je,
         const ArrayBase<T> &jf)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, jd));
    eigen_assert(MATRIX_SAME_SIZE(ja, je));
    eigen_assert(MATRIX_SAME_SIZE(ja, jf));

    typedef typename wigner6j_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner6j_functor<T>(ja.derived(),
                                                      jb.derived(),
                                                      jc.derived(),
                                                      jd.derived(),
                                                      je.derived(),
                                                      jf.derived()));
}

template <typename T, typename U>
inline double wigner6j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T jd,
                              const T je,
                              const T jf,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner6j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int jd,
                              const int je,
                              const int jf,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_6j_e(ja, jb, jc, jd, je, jf, &r) == GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("wigner6j"));
}

template <typename T, typename U>
class wigner6j_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner6j_e_functor(const T &ja,
                       const T &jb,
                       const T &jc,
                       const T &jd,
                       const T &je,
                       const T &jf,
                       U &e)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_jd(jd)
        , m_je(je)
        , m_jf(jf)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner6j_e_impl(m_ja(i, j),
                               m_jb(i, j),
                               m_jc(i, j),
                               m_jd(i, j),
                               m_je(i, j),
                               m_jf(i, j),
                               m_e(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_jd;
    const T &m_je;
    const T &m_jf;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner6j_e_functor<T, U>,
                      typename wigner6j_e_functor<T, U>::ArrayType>
wigner6j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &jd,
         const ArrayBase<T> &je,
         const ArrayBase<T> &jf,
         ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, jd));
    eigen_assert(MATRIX_SAME_SIZE(ja, je));
    eigen_assert(MATRIX_SAME_SIZE(ja, jf));

    typedef typename wigner6j_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner6j_e_functor<T, U>(ja.derived(),
                                                           jb.derived(),
                                                           jc.derived(),
                                                           jd.derived(),
                                                           je.derived(),
                                                           jf.derived(),
                                                           e.derived()));
}

// ========================================
// wigner 9j syjeol
// ========================================

template <typename T>
inline double wigner9j_impl(const T ja,
                            const T jb,
                            const T jc,
                            const T jd,
                            const T je,
                            const T jf,
                            const T jg,
                            const T jh,
                            const T ji)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner9j_impl(const int ja,
                            const int jb,
                            const int jc,
                            const int jd,
                            const int je,
                            const int jf,
                            const int jg,
                            const int jh,
                            const int ji)
{
    return gsl_sf_coupling_9j(ja, jb, jc, jd, je, jf, jg, jh, ji);
}

template <typename T>
class wigner9j_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner9j_functor(const T &ja,
                     const T &jb,
                     const T &jc,
                     const T &jd,
                     const T &je,
                     const T &jf,
                     const T &jg,
                     const T &jh,
                     const T &ji)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_jd(jd)
        , m_je(je)
        , m_jf(jf)
        , m_jg(jg)
        , m_jh(jh)
        , m_ji(ji)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner9j_impl(m_ja(i, j),
                             m_jb(i, j),
                             m_jc(i, j),
                             m_jd(i, j),
                             m_je(i, j),
                             m_jf(i, j),
                             m_jg(i, j),
                             m_jh(i, j),
                             m_ji(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_jd;
    const T &m_je;
    const T &m_jf;
    const T &m_jg;
    const T &m_jh;
    const T &m_ji;
};

template <typename T>
inline CwiseNullaryOp<wigner9j_functor<T>,
                      typename wigner9j_functor<T>::ArrayType>
wigner9j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &jd,
         const ArrayBase<T> &je,
         const ArrayBase<T> &jf,
         const ArrayBase<T> &jg,
         const ArrayBase<T> &jh,
         const ArrayBase<T> &ji)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, jd));
    eigen_assert(MATRIX_SAME_SIZE(ja, je));
    eigen_assert(MATRIX_SAME_SIZE(ja, jf));
    eigen_assert(MATRIX_SAME_SIZE(ja, jg));
    eigen_assert(MATRIX_SAME_SIZE(ja, jh));
    eigen_assert(MATRIX_SAME_SIZE(ja, ji));

    typedef typename wigner9j_functor<T>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner9j_functor<T>(ja.derived(),
                                                      jb.derived(),
                                                      jc.derived(),
                                                      jd.derived(),
                                                      je.derived(),
                                                      jf.derived(),
                                                      jg.derived(),
                                                      jh.derived(),
                                                      ji.derived()));
}

template <typename T, typename U>
inline double wigner9j_e_impl(const T ja,
                              const T jb,
                              const T jc,
                              const T jd,
                              const T je,
                              const T jf,
                              const T jg,
                              const T jh,
                              const T ji,
                              U &e)
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline double wigner9j_e_impl(const int ja,
                              const int jb,
                              const int jc,
                              const int jd,
                              const int je,
                              const int jf,
                              const int jg,
                              const int jh,
                              const int ji,
                              double &e)
{
    gsl_sf_result r;
    if (gsl_sf_coupling_9j_e(ja, jb, jc, jd, je, jf, jg, jh, ji, &r) ==
        GSL_SUCCESS) {
        e = r.err;
        return r.val;
    }
    THROW_OR_RETURN_NAN(std::runtime_error("wigner9j"));
}

template <typename T, typename U>
class wigner9j_e_functor
{
  public:
    typedef Array<double,
                  T::RowsAtCompileTime,
                  T::ColsAtCompileTime,
                  T::Flags & RowMajorBit ? RowMajor : ColMajor,
                  T::MaxRowsAtCompileTime,
                  T::MaxColsAtCompileTime>
        ArrayType;

    wigner9j_e_functor(const T &ja,
                       const T &jb,
                       const T &jc,
                       const T &jd,
                       const T &je,
                       const T &jf,
                       const T &jg,
                       const T &jh,
                       const T &ji,
                       U &e)
        : m_ja(ja)
        , m_jb(jb)
        , m_jc(jc)
        , m_jd(jd)
        , m_je(je)
        , m_jf(jf)
        , m_jg(jg)
        , m_jh(jh)
        , m_ji(ji)
        , m_e(e)
    {
    }

    const double operator()(Index i, Index j) const
    {
        return wigner9j_e_impl(m_ja(i, j),
                               m_jb(i, j),
                               m_jc(i, j),
                               m_jd(i, j),
                               m_je(i, j),
                               m_jf(i, j),
                               m_jg(i, j),
                               m_jh(i, j),
                               m_ji(i, j),
                               m_e(i, j));
    }

  private:
    const T &m_ja;
    const T &m_jb;
    const T &m_jc;
    const T &m_jd;
    const T &m_je;
    const T &m_jf;
    const T &m_jg;
    const T &m_jh;
    const T &m_ji;
    U &m_e;
};

template <typename T, typename U>
inline CwiseNullaryOp<wigner9j_e_functor<T, U>,
                      typename wigner9j_e_functor<T, U>::ArrayType>
wigner9j(const ArrayBase<T> &ja,
         const ArrayBase<T> &jb,
         const ArrayBase<T> &jc,
         const ArrayBase<T> &jd,
         const ArrayBase<T> &je,
         const ArrayBase<T> &jf,
         const ArrayBase<T> &jg,
         const ArrayBase<T> &jh,
         const ArrayBase<T> &ji,
         ArrayBase<U> &e)
{
    eigen_assert(MATRIX_SAME_SIZE(ja, jb));
    eigen_assert(MATRIX_SAME_SIZE(ja, jc));
    eigen_assert(MATRIX_SAME_SIZE(ja, jd));
    eigen_assert(MATRIX_SAME_SIZE(ja, je));
    eigen_assert(MATRIX_SAME_SIZE(ja, jf));
    eigen_assert(MATRIX_SAME_SIZE(ja, jg));
    eigen_assert(MATRIX_SAME_SIZE(ja, jh));
    eigen_assert(MATRIX_SAME_SIZE(ja, ji));

    typedef typename wigner9j_e_functor<T, U>::ArrayType ArrayType;
    return ArrayType::NullaryExpr(ja.rows(),
                                  ja.cols(),
                                  wigner9j_e_functor<T, U>(ja.derived(),
                                                           jb.derived(),
                                                           jc.derived(),
                                                           jd.derived(),
                                                           je.derived(),
                                                           jf.derived(),
                                                           jg.derived(),
                                                           jh.derived(),
                                                           ji.derived(),
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
