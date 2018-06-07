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

#ifndef __IEXP_GSLVEC__
#define __IEXP_GSLVEC__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_vector.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename Scalar>
class gslvec;

#define DEFINE_GSLVEC_SCALAR(scalar,                                           \
                             data_scalar,                                      \
                             set_scalar,                                       \
                             t_vector,                                         \
                             t_block,                                          \
                             f_alloc,                                          \
                             f_free,                                           \
                             f_set)                                            \
    template <>                                                                \
    class gslvec<scalar>                                                       \
    {                                                                          \
      public:                                                                  \
        template <typename T>                                                  \
        gslvec(const DenseBase<T> &v,                                          \
               bool copy = true,                                               \
               typename std::enable_if<bool(T::Flags & DirectAccessBit)>::type \
                   * = nullptr)                                                \
        {                                                                      \
            static_assert(TYPE_IS(typename T::Scalar, scalar),                 \
                          "must be same scalar");                              \
                                                                               \
            if (copy) {                                                        \
                copy_vector(v);                                                \
            } else {                                                           \
                m_result = &m_vec;                                             \
                                                                               \
                m_block.size = v.size();                                       \
                m_block.data = const_cast<data_scalar *>(                      \
                    reinterpret_cast<const data_scalar *>(                     \
                        v.derived().data()));                                  \
                                                                               \
                m_vec.size = v.size();                                         \
                m_vec.stride = 1;                                              \
                m_vec.data = m_block.data;                                     \
                m_vec.block = &m_block;                                        \
                m_vec.owner = 0;                                               \
            }                                                                  \
        }                                                                      \
                                                                               \
        template <typename T>                                                  \
        gslvec(const DenseBase<T> &v,                                          \
               typename std::enable_if<                                        \
                   !bool(T::Flags & DirectAccessBit)>::type * = nullptr)       \
        {                                                                      \
            static_assert(TYPE_IS(typename T::Scalar, scalar),                 \
                          "must be same scalar");                              \
                                                                               \
            copy_vector(v);                                                    \
        }                                                                      \
                                                                               \
        ~gslvec()                                                              \
        {                                                                      \
            if (m_result != &m_vec) {                                          \
                f_free(m_result);                                              \
            }                                                                  \
        }                                                                      \
                                                                               \
        t_vector *gsl_vector()                                                 \
        {                                                                      \
            return m_result;                                                   \
        }                                                                      \
                                                                               \
      private:                                                                 \
        ::t_vector *m_result;                                                  \
        t_block m_block;                                                       \
        ::t_vector m_vec;                                                      \
                                                                               \
        template <typename T>                                                  \
        void copy_vector(const DenseBase<T> &v)                                \
        {                                                                      \
            m_result = f_alloc(v.size());                                      \
            IEXP_NOT_NULLPTR(m_result);                                        \
            for (size_t i = 0; i < v.size(); ++i) {                            \
                typename T::Scalar &&val =                                     \
                    const_cast<typename T::Scalar &&>(v[i]);                   \
                f_set(m_result, i, *(set_scalar *)&val);                       \
            }                                                                  \
        }                                                                      \
    };

#define DEFINE_GSLVEC(scalar, t_vector, t_block, f_alloc, f_free, f_set)       \
    DEFINE_GSLVEC_SCALAR(scalar,                                               \
                         scalar,                                               \
                         scalar,                                               \
                         t_vector,                                             \
                         t_block,                                              \
                         f_alloc,                                              \
                         f_free,                                               \
                         f_set)

DEFINE_GSLVEC_SCALAR(std::complex<long double>,
                     long double,
                     gsl_complex_long_double,
                     gsl_vector_complex_long_double,
                     gsl_block_complex_long_double,
                     gsl_vector_complex_long_double_alloc,
                     gsl_vector_complex_long_double_free,
                     gsl_vector_complex_long_double_set)

DEFINE_GSLVEC_SCALAR(std::complex<double>,
                     double,
                     gsl_complex,
                     gsl_vector_complex,
                     gsl_block_complex,
                     gsl_vector_complex_alloc,
                     gsl_vector_complex_free,
                     gsl_vector_complex_set)

DEFINE_GSLVEC_SCALAR(std::complex<float>,
                     float,
                     gsl_complex_float,
                     gsl_vector_complex_float,
                     gsl_block_complex_float,
                     gsl_vector_complex_float_alloc,
                     gsl_vector_complex_float_free,
                     gsl_vector_complex_float_set)

DEFINE_GSLVEC(long double,
              gsl_vector_long_double,
              gsl_block_long_double,
              gsl_vector_long_double_alloc,
              gsl_vector_long_double_free,
              gsl_vector_long_double_set)

DEFINE_GSLVEC(double,
              gsl_vector,
              gsl_block,
              gsl_vector_alloc,
              gsl_vector_free,
              gsl_vector_set)

DEFINE_GSLVEC(float,
              gsl_vector_float,
              gsl_block_float,
              gsl_vector_float_alloc,
              gsl_vector_float_free,
              gsl_vector_float_set)

DEFINE_GSLVEC(unsigned long,
              gsl_vector_ulong,
              gsl_block_ulong,
              gsl_vector_ulong_alloc,
              gsl_vector_ulong_free,
              gsl_vector_ulong_set)

DEFINE_GSLVEC(long,
              gsl_vector_long,
              gsl_block_long,
              gsl_vector_long_alloc,
              gsl_vector_long_free,
              gsl_vector_long_set)

DEFINE_GSLVEC(unsigned int,
              gsl_vector_uint,
              gsl_block_uint,
              gsl_vector_uint_alloc,
              gsl_vector_uint_free,
              gsl_vector_uint_set)

DEFINE_GSLVEC(int,
              gsl_vector_int,
              gsl_block_int,
              gsl_vector_int_alloc,
              gsl_vector_int_free,
              gsl_vector_int_set)

DEFINE_GSLVEC(unsigned short,
              gsl_vector_ushort,
              gsl_block_ushort,
              gsl_vector_ushort_alloc,
              gsl_vector_ushort_free,
              gsl_vector_ushort_set)

DEFINE_GSLVEC(short,
              gsl_vector_short,
              gsl_block_short,
              gsl_vector_short_alloc,
              gsl_vector_short_free,
              gsl_vector_short_set)

DEFINE_GSLVEC(unsigned char,
              gsl_vector_uchar,
              gsl_block_uchar,
              gsl_vector_uchar_alloc,
              gsl_vector_uchar_free,
              gsl_vector_uchar_set)

DEFINE_GSLVEC(char,
              gsl_vector_char,
              gsl_block_char,
              gsl_vector_char_alloc,
              gsl_vector_char_free,
              gsl_vector_char_set)

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_GSLVEC__ */
