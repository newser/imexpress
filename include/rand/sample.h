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

#ifndef __IEXP_RAND_SAMPLE__
#define __IEXP_RAND_SAMPLE__

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
class sample_functor
{
  public:
    using Scalar = typename T::Scalar;
    using ResultType = typename dense_derive<T,
                                             Scalar,
                                             row_form ? 1 : Dynamic,
                                             row_form ? Dynamic : 1>::type;

    sample_functor(const T &x, size_t k, rng::type type, unsigned long seed)
        : m_result(new Scalar[k])
    {
        rng r(type, seed);
        typename type_eval<T>::type m_x(x.eval());
        gsl_ran_sample(r.gsl(),
                       m_result.get(),
                       k,
                       const_cast<Scalar *>(m_x.data()),
                       m_x.size(),
                       sizeof(Scalar));
    }

    Scalar operator()(Index i) const
    {
        return m_result.get()[i];
    }

  private:
    std::shared_ptr<Scalar> m_result;
};

template <bool row_form = false, typename T = void>
inline CwiseNullaryOp<sample_functor<row_form, T>,
                      typename sample_functor<row_form, T>::ResultType>
sample(const DenseBase<T> &x,
       size_t k,
       unsigned long seed = 0,
       rng::type type = DEFAULT_RNG_TYPE)
{
    using ResultType = typename sample_functor<row_form, T>::ResultType;
    return ResultType::NullaryExpr(k,
                                   sample_functor<row_form, T>(x.derived(),
                                                               k,
                                                               type,
                                                               seed));
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_SAMPLE__ */
