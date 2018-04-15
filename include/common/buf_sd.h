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

#ifndef __IEXP_BUF_SD__
#define __IEXP_BUF_SD__

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

template <typename T, typename Scalar = TP1(T)>
class buf_sd
{
  public:
    buf_sd(size_t n)
        : m_size(n)
        , m_data(m_static)
    {
        if (m_size > (sizeof(m_static) / sizeof(Scalar))) {
            m_data = new Scalar[n];
            IEXP_NOT_NULLPTR(m_data);
        }
    }

    buf_sd(const buf_sd &&x)
        : m_size(x.size())
        , m_data(m_static)
    {
    }

    ~buf_sd()
    {
        if (!is_static()) {
            delete[] m_data;
        }
    }

    Scalar *data() const
    {
        return m_data;
    }

    size_t size() const
    {
        return m_size;
    }

    bool is_static() const
    {
        return m_data == m_static;
    }

  private:
    buf_sd(const buf_sd &) = delete;
    buf_sd &operator=(const buf_sd &) = delete;

    size_t m_size;
    Scalar *m_data;
    Scalar m_static[is_dynamic<T>::value ? 1 : (TP5(T) * TP6(T))];
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_BUF_SD__ */
