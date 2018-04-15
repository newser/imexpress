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

#ifndef __IEXP_BUF_RC__
#define __IEXP_BUF_RC__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

template <typename Scalar, bool row_major>
class buf_rc
{
  public:
    buf_rc(size_t row, size_t col)
        : m_row(row)
        , m_col(col)
        , m_data(new Scalar[row * col])
    {
    }

    Scalar *data() const
    {
        return m_data.get();
    }

    size_t size() const
    {
        return m_row * m_col;
    }

    Scalar &operator[](size_t i) const
    {
        eigen_assert(i < size());
        return m_data.get()[i];
    }

    Scalar &operator()(size_t i) const
    {
        eigen_assert(i < size());
        return m_data.get()[i];
    }

    Scalar &operator()(size_t i, size_t j) const
    {
        eigen_assert((i < m_row) && (j < m_col));
        return at(i, j, TYPE_BOOL(row_major)());
    }

  private:
    size_t m_row, m_col;
    std::shared_ptr<Scalar> m_data;

    Scalar &at(Index i, Index j, std::false_type) const
    {
        // col major
        return m_data.get()[i + j * m_row];
    }

    Scalar &at(Index i, Index j, std::true_type) const
    {
        // row major
        return m_data.get()[i * m_col + j];
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_BUF_RC__ */
