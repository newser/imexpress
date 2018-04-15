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

#ifndef __IEXP_COPY__
#define __IEXP_COPY__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <string.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename DST, typename SRC>
struct copy
{
    copy(DST *d, const SRC *s, size_t n)
    {
        for (size_t i = 0; i < n; ++i) {
            d[i] = static_cast<DST>(s[i]);
        }
    }
};

template <typename T>
struct copy<T, T>
{
    copy(T *d, const T *s, size_t n)
    {
        ::memcpy(d, s, sizeof(T) * n);
    }
};

template <typename D, typename S>
void copy_matrix(D *d, size_t n, const DenseBase<S> &s, std::true_type)
{
    eigen_assert(n == s.size());

    // typename type_eval<S>::type m_s(s.derived().eval());
    int idx = 0;
    for (int i = 0; i < s.rows(); ++i) {
        for (int j = 0; j < s.cols(); ++j) {
            d[idx++] = s(i, j);
        }
    }
}

template <typename D, typename S>
void copy_matrix(D *d, size_t n, const DenseBase<S> &s, std::false_type)
{
    eigen_assert(n == s.size());

    // typename type_eval<S>::type m_s(s.derived().eval());
    int idx = 0;
    for (int j = 0; j < s.cols(); ++j) {
        for (int i = 0; i < s.rows(); ++i) {
            d[idx++] = s(i, j);
        }
    }
}

template <typename D, typename S>
void copy_matrix(D *d, size_t n, const DenseBase<S> &s)
{
    copy_matrix(d, n, s, TYPE_BOOL(bool(S::Flags & RowMajorBit))());
}

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_COPY__ */
