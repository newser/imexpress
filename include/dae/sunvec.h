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
 * You should have received v copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

#ifndef __IEXP_DAE_SUNVEC__
#define __IEXP_DAE_SUNVEC__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>

IEXP_NS_BEGIN

namespace dae {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

extern struct _generic_N_Vector_Ops sunvec_ops;

class sunvec_serial
{
  public:
    template <typename T>
    sunvec_serial(const DenseBase<T> &v,
                  bool copy = true,
                  typename std::enable_if<bool(T::Flags &DirectAccessBit)>::type
                      * = nullptr)
    {
        static_assert(TYPE_IS(typename T::Scalar, realtype),
                      "only support realtype scalar");

        if (copy) {
            copy_vector(v);
        } else {
            m_result = &m_vec;

            m_serial.length = v.size();
            m_serial.own_data = SUNFALSE;
            m_serial.data = const_cast<realtype *>(v.derived().data());

            m_vec.content = &m_serial;
            m_vec.ops = &sunvec_ops;
        }
    }

    template <typename T>
    sunvec_serial(
        const DenseBase<T> &v,
        typename std::enable_if<!bool(T::Flags & DirectAccessBit)>::type * =
            nullptr)
    {
        static_assert(TYPE_IS(typename T::Scalar, realtype),
                      "only support realtype scalar");

        copy_vector(v);
    }

    ~sunvec_serial()
    {
        if (m_result != &m_vec) {
            N_VDestroy(m_result);
        }
    }

    N_Vector n_vector()
    {
        return m_result;
    }

  private:
    N_Vector m_result;
    struct _N_VectorContent_Serial m_serial;
    struct _generic_N_Vector m_vec;

    template <typename T>
    void copy_vector(const DenseBase<T> &v)
    {
        m_result = N_VNew_Serial(v.size());
        IEXP_NOT_NULLPTR(m_result);
        for (sunindextype i = 0; i < v.size(); ++i) {
            NV_Ith_S(m_result, i) = v(i);
        }
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_DAE_SUNVEC__ */
