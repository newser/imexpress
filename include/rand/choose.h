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

template <typename T>
class choose_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    choose_functor(const T &x,
                   const size_t k,
                   unsigned long seed,
                   rng_type type)
        : m_result(k)
    {
        rng r(type, seed);
        typename type_eval<T>::type m_x(x.eval());
        gsl_ran_choose(r.gsl(),
                       m_result.data(),
                       k,
                       m_x.data(),
                       m_x.size(),
                       sizeof(typename T::Scalar));
    }

    const typename T::Scalar &operator()(Index i) const
    {
        return m_result[i];
    }

  private:
    ArrayType m_result;
};

template <typename T>
inline CwiseNullaryOp<choose_functor<T>, typename choose_functor<T>::ArrayType>
choose(const ArrayBase<T> &x,
       const size_t k,
       unsigned long seed = 0,
       rng_type type = DEFAULT_RNG)
{
    eigen_assert(IS_VEC(x));
    eigen_assert(k <= x.size());

    using ArrayType = typename choose_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(k,
                                  choose_functor<T>(x.derived(),
                                                    k,
                                                    seed,
                                                    type));
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
