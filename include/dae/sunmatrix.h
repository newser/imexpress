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

#ifndef __IEXP_DAE_SUNMATRIX__
#define __IEXP_DAE_SUNMATRIX__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <sundials/sundials_matrix.h>
#include <sunmatrix/sunmatrix_dense.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename T>
class sunmat_dense
{
  public:
    sunmat_dense(const DenseBase<T> &a)
        : m_eval(a.derived().eval())
    {
        m_dense.M = a.rows();
        m_dense.N = a.cols();
        m_dense.data = m_eval.data();
        m_dense.ldata = a.size();
        m_dense.cols = new realtype *[a.cols()];
        for (sunindextype j = 0; j < a.cols(); ++j) {
            m_dense.cols[j] = m_eval.data() + j * a.rows();
        }

        m_mat.content = &m_dense;
        m_mat.ops = &s_ops;
    }

    ~sunmat_dense()
    {
    }

  private:
    static struct _generic_SUNMatrix_Ops s_ops;

    typename type_eval<typename T::Derived>::type m_eval;
    struct _SUNMatrixContent_Dense m_dense;
    struct _generic_SUNMatrix m_mat;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_SUNMATRIX__ */
