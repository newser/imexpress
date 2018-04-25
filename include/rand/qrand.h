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

#ifndef __IEXP_QRAND__
#define __IEXP_QRAND__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/qrng.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

template <typename T>
inline void qrand_impl(qrng &r, T x[])
{
    UNSUPPORTED_TYPE(T);
}

template <>
inline void qrand_impl<double>(qrng &r, double x[])
{
    return r.next(x);
}

template <typename T>
inline auto qrand(DenseBase<T> &x, qrng::type type) -> decltype(x.derived())
{
    unsigned int dim;
    if ((x.rows() == 1) || (x.cols() == 1)) {
        dim = 1;
    } else if (x.IsRowMajor) {
        dim = x.cols();
    } else {
        dim = x.rows();
    }
    qrng r(type, dim);

    return qrand(x, r);
}

template <typename T>
inline auto qrand(DenseBase<T> &x, qrng &r) -> decltype(x.derived())
{
    typename T::Scalar *data = x.derived().data();
    unsigned int dim = r.dim();
    if ((x.rows() == 1) || (x.cols() == 1)) {
        eigen_assert(dim == 1);
        for (Index i = 0; i < x.size(); ++i) {
            qrand_impl<typename T::Scalar>(r, &data[i]);
        }
    } else if (x.IsRowMajor) {
        eigen_assert(dim == x.cols());
        for (Index i = 0; i < x.rows(); ++i) {
            qrand_impl<typename T::Scalar>(r, &data[i * dim]);
        }
    } else {
        eigen_assert(dim == x.rows());
        for (Index i = 0; i < x.cols(); ++i) {
            qrand_impl<typename T::Scalar>(r, &data[i * dim]);
        }
    }
    return x.derived();
}
}

IEXP_NS_END

#endif /* __IEXP_QRAND__ */
