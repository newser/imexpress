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

#ifndef __IEXP_EXP_INTEGRAL__
#define __IEXP_EXP_INTEGRAL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>

IEXP_NS_BEGIN

namespace sf {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

// ========================================
// expint_1
// ========================================

template <typename T>
class expint_1_functor : public functor_foreach<expint_1_functor<T>, T, double>
{
  public:
    expint_1_functor(const T &x)
        : functor_foreach<expint_1_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_E1(x);
    }
};

template <typename T>
inline CwiseNullaryOp<expint_1_functor<T>,
                      typename expint_1_functor<T>::ResultType>
expint_1(const DenseBase<T> &x)
{
    using ResultType = typename expint_1_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_1_functor<T>(x.derived()));
}

template <typename T, typename U>
class expint_1_e_functor
    : public functor_foreach_e<expint_1_e_functor<T, U>, T, U, double>
{
  public:
    expint_1_e_functor(const T &x, U &e)
        : functor_foreach_e<expint_1_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_E1_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_1_e_functor<T, U>,
                      typename expint_1_e_functor<T, U>::ResultType>
expint_1(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_1_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_1_e_functor<T, U>(x.derived(),
                                                            e.derived()));
}

// ========================================
// expint_2
// ========================================

template <typename T>
class expint_2_functor : public functor_foreach<expint_2_functor<T>, T, double>
{
  public:
    expint_2_functor(const T &x)
        : functor_foreach<expint_2_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_E2(x);
    }
};

template <typename T>
inline CwiseNullaryOp<expint_2_functor<T>,
                      typename expint_2_functor<T>::ResultType>
expint_2(const DenseBase<T> &x)
{
    using ResultType = typename expint_2_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_2_functor<T>(x.derived()));
}

template <typename T, typename U>
class expint_2_e_functor
    : public functor_foreach_e<expint_2_e_functor<T, U>, T, U, double>
{
  public:
    expint_2_e_functor(const T &x, U &e)
        : functor_foreach_e<expint_2_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_E2_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_2_e_functor<T, U>,
                      typename expint_2_e_functor<T, U>::ResultType>
expint_2(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_2_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_2_e_functor<T, U>(x.derived(),
                                                            e.derived()));
}

// ========================================
// expint
// ========================================

template <typename T>
class expint_functor : public functor_foreach<expint_functor<T>, T, double>
{
  public:
    expint_functor(int n, const T &x)
        : functor_foreach<expint_functor<T>, T, double>(x)
        , m_n(n)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_En(m_n, x);
    }

  private:
    int m_n;
};

template <typename T>
inline CwiseNullaryOp<expint_functor<T>, typename expint_functor<T>::ResultType>
expint(int n, const DenseBase<T> &x)
{
    using ResultType = typename expint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_functor<T>(n, x.derived()));
}

template <typename T, typename U>
class expint_e_functor
    : public functor_foreach_e<expint_e_functor<T, U>, T, U, double>
{
  public:
    expint_e_functor(int n, const T &x, U &e)
        : functor_foreach_e<expint_e_functor<T, U>, T, U, double>(x, e)
        , m_n(n)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_En_e(m_n, x, &r);
        e = r.err;
        return r.val;
    }

  private:
    int m_n;
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_e_functor<T, U>,
                      typename expint_e_functor<T, U>::ResultType>
expint(int n, const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_e_functor<T, U>(n,
                                                          x.derived(),
                                                          e.derived()));
}

// ========================================
// expint_ei
// ========================================

template <typename T>
class expint_ei_functor
    : public functor_foreach<expint_ei_functor<T>, T, double>
{
  public:
    expint_ei_functor(const T &x)
        : functor_foreach<expint_ei_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_Ei(x);
    }
};

template <typename T>
inline CwiseNullaryOp<expint_ei_functor<T>,
                      typename expint_ei_functor<T>::ResultType>
expint_ei(const DenseBase<T> &x)
{
    using ResultType = typename expint_ei_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_ei_functor<T>(x.derived()));
}

template <typename T, typename U>
class expint_ei_e_functor
    : public functor_foreach_e<expint_ei_e_functor<T, U>, T, U, double>
{
  public:
    expint_ei_e_functor(const T &x, U &e)
        : functor_foreach_e<expint_ei_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_Ei_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_ei_e_functor<T, U>,
                      typename expint_ei_e_functor<T, U>::ResultType>
expint_ei(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_ei_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_ei_e_functor<T, U>(x.derived(),
                                                             e.derived()));
}

// ========================================
// expint_ei3
// ========================================

template <typename T>
class expint_ei3_functor
    : public functor_foreach<expint_ei3_functor<T>, T, double>
{
  public:
    expint_ei3_functor(const T &x)
        : functor_foreach<expint_ei3_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_expint_3(x);
    }
};

template <typename T>
inline CwiseNullaryOp<expint_ei3_functor<T>,
                      typename expint_ei3_functor<T>::ResultType>
expint_ei3(const DenseBase<T> &x)
{
    using ResultType = typename expint_ei3_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_ei3_functor<T>(x.derived()));
}

template <typename T, typename U>
class expint_ei3_e_functor
    : public functor_foreach_e<expint_ei3_e_functor<T, U>, T, U, double>
{
  public:
    expint_ei3_e_functor(const T &x, U &e)
        : functor_foreach_e<expint_ei3_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_expint_3_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<expint_ei3_e_functor<T, U>,
                      typename expint_ei3_e_functor<T, U>::ResultType>
expint_ei3(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename expint_ei3_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   expint_ei3_e_functor<T, U>(x.derived(),
                                                              e.derived()));
}

// ========================================
// sinh integral
// ========================================

template <typename T>
class sinhint_functor : public functor_foreach<sinhint_functor<T>, T, double>
{
  public:
    sinhint_functor(const T &x)
        : functor_foreach<sinhint_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_Shi(x);
    }
};

template <typename T>
inline CwiseNullaryOp<sinhint_functor<T>,
                      typename sinhint_functor<T>::ResultType>
sinhint(const DenseBase<T> &x)
{
    using ResultType = typename sinhint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sinhint_functor<T>(x.derived()));
}

template <typename T, typename U>
class sinhint_e_functor
    : public functor_foreach_e<sinhint_e_functor<T, U>, T, U, double>
{
  public:
    sinhint_e_functor(const T &x, U &e)
        : functor_foreach_e<sinhint_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_Shi_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<sinhint_e_functor<T, U>,
                      typename sinhint_e_functor<T, U>::ResultType>
sinhint(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sinhint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sinhint_e_functor<T, U>(x.derived(),
                                                           e.derived()));
}

// ========================================
// cosh integral
// ========================================

template <typename T>
class coshint_functor : public functor_foreach<coshint_functor<T>, T, double>
{
  public:
    coshint_functor(const T &x)
        : functor_foreach<coshint_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_Chi(x);
    }
};

template <typename T>
inline CwiseNullaryOp<coshint_functor<T>,
                      typename coshint_functor<T>::ResultType>
coshint(const DenseBase<T> &x)
{
    using ResultType = typename coshint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   coshint_functor<T>(x.derived()));
}

template <typename T, typename U>
class coshint_e_functor
    : public functor_foreach_e<coshint_e_functor<T, U>, T, U, double>
{
  public:
    coshint_e_functor(const T &x, U &e)
        : functor_foreach_e<coshint_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_Chi_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<coshint_e_functor<T, U>,
                      typename coshint_e_functor<T, U>::ResultType>
coshint(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename coshint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   coshint_e_functor<T, U>(x.derived(),
                                                           e.derived()));
}

// ========================================
// sin integral
// ========================================

template <typename T>
class sinint_functor : public functor_foreach<sinint_functor<T>, T, double>
{
  public:
    sinint_functor(const T &x)
        : functor_foreach<sinint_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_Si(x);
    }
};

template <typename T>
inline CwiseNullaryOp<sinint_functor<T>, typename sinint_functor<T>::ResultType>
sinint(const DenseBase<T> &x)
{
    using ResultType = typename sinint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sinint_functor<T>(x.derived()));
}

template <typename T, typename U>
class sinint_e_functor
    : public functor_foreach_e<sinint_e_functor<T, U>, T, U, double>
{
  public:
    sinint_e_functor(const T &x, U &e)
        : functor_foreach_e<sinint_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_Si_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<sinint_e_functor<T, U>,
                      typename sinint_e_functor<T, U>::ResultType>
sinint(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename sinint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   sinint_e_functor<T, U>(x.derived(),
                                                          e.derived()));
}

// ========================================
// cos integral
// ========================================

template <typename T>
class cosint_functor : public functor_foreach<cosint_functor<T>, T, double>
{
  public:
    cosint_functor(const T &x)
        : functor_foreach<cosint_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_Ci(x);
    }
};

template <typename T>
inline CwiseNullaryOp<cosint_functor<T>, typename cosint_functor<T>::ResultType>
cosint(const DenseBase<T> &x)
{
    using ResultType = typename cosint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   cosint_functor<T>(x.derived()));
}

template <typename T, typename U>
class cosint_e_functor
    : public functor_foreach_e<cosint_e_functor<T, U>, T, U, double>
{
  public:
    cosint_e_functor(const T &x, U &e)
        : functor_foreach_e<cosint_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_Ci_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<cosint_e_functor<T, U>,
                      typename cosint_e_functor<T, U>::ResultType>
cosint(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename cosint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   cosint_e_functor<T, U>(x.derived(),
                                                          e.derived()));
}

// ========================================
// cos integral
// ========================================

template <typename T>
class atanint_functor : public functor_foreach<atanint_functor<T>, T, double>
{
  public:
    atanint_functor(const T &x)
        : functor_foreach<atanint_functor<T>, T, double>(x)
    {
    }

    double foreach_impl(double x) const
    {
        return gsl_sf_atanint(x);
    }
};

template <typename T>
inline CwiseNullaryOp<atanint_functor<T>,
                      typename atanint_functor<T>::ResultType>
atanint(const DenseBase<T> &x)
{
    using ResultType = typename atanint_functor<T>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   atanint_functor<T>(x.derived()));
}

template <typename T, typename U>
class atanint_e_functor
    : public functor_foreach_e<atanint_e_functor<T, U>, T, U, double>
{
  public:
    atanint_e_functor(const T &x, U &e)
        : functor_foreach_e<atanint_e_functor<T, U>, T, U, double>(x, e)
    {
    }

    double foreach_e_impl(double x, double &e) const
    {
        gsl_sf_result r;
        gsl_sf_atanint_e(x, &r);
        e = r.err;
        return r.val;
    }
};

template <typename T, typename U>
inline CwiseNullaryOp<atanint_e_functor<T, U>,
                      typename atanint_e_functor<T, U>::ResultType>
atanint(const DenseBase<T> &x, DenseBase<U> &e)
{
    using ResultType = typename atanint_e_functor<T, U>::ResultType;
    return ResultType::NullaryExpr(x.rows(),
                                   x.cols(),
                                   atanint_e_functor<T, U>(x.derived(),
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

#endif /* __IEXP_EXP_INTEGRAL__ */
