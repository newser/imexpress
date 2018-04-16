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

#ifndef __IEXP_RAND_CHOOSE__
#define __IEXP_RAND_CHOOSE__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <bool row_form, typename T>
class choose_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    choose_functor(const T &x, size_t k, rng &r)
        : m_result(new Scalar[k])
    {
        typename type_eval<T>::type m_x(x.eval());
        gsl_ran_choose(r.gsl(),
                       m_result.get(),
                       k,
                       const_cast<Scalar *>(m_x.data()),
                       m_x.size(),
                       sizeof(typename T::Scalar));
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<choose_functor<row_form, T>,
                      typename choose_functor<row_form, T>::ResultType>
choose(const DenseBase<T> &x,
       size_t k,
       unsigned long seed = 0,
       rng::type type = DEFAULT_RNG_TYPE)
{
    rng r(type, seed);
    return choose<row_form, T>(x, k, r);
}

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<choose_functor<row_form, T>,
                      typename choose_functor<row_form, T>::ResultType>
choose(const DenseBase<T> &x, size_t k, rng &r)
{
    eigen_assert(k <= x.size());

    using ResultType = typename choose_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(row_form ? 1 : k,
                                   row_form ? k : 1,
                                   choose_functor<row_form, T>(x.derived(),
                                                               k,
                                                               r));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_CHOOSE__ */
