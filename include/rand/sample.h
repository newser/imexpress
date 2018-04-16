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

template <typename T>
class sample_functor
{
  public:
    using ArrayType = Array<typename T::Scalar,
                            T::RowsAtCompileTime,
                            T::ColsAtCompileTime,
                            T::Flags & RowMajorBit ? RowMajor : ColMajor,
                            T::MaxRowsAtCompileTime,
                            T::MaxColsAtCompileTime>;

    sample_functor(const T &x, size_t k, unsigned long seed, rng::type type)
        : m_result(k)
    {
        rng r(type, seed);
        typename type_eval<T>::type m_x(x.eval());
        gsl_ran_sample(r.gsl(),
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
inline CwiseNullaryOp<sample_functor<T>, typename sample_functor<T>::ArrayType>
sample(const ArrayBase<T> &x,
       size_t k,
       unsigned long seed = 0,
       rng::type type = DEFAULT_RNG_TYPE)
{
    eigen_assert(IS_VEC(x));

    using ArrayType = typename sample_functor<T>::ArrayType;
    return ArrayType::NullaryExpr(k,
                                  sample_functor<T>(x.derived(),
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

#endif /* __IEXP_RAND_SAMPLE__ */
