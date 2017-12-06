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

#ifndef __IEXP_MULTISET__
#define __IEXP_MULTISET__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_multiset.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class multiset
{
  public:
    using index = size_t;

    enum class init
    {
        NONE,
        FIRST,
        LAST,
    };

    multiset(const index n, const index k, init i = init::NONE)
    {
        m_gm = gsl_multiset_alloc(n, k);
        IEXP_NOT_NULLPTR(m_gm);
        if (i == init::FIRST) {
            gsl_multiset_init_first(m_gm);
        } else if (i == init::LAST) {
            gsl_multiset_init_last(m_gm);
        }
    }

    multiset(const multiset &c)
    {
        m_gm = gsl_multiset_alloc(c.n(), c.k());
        IEXP_NOT_NULLPTR(m_gm);
        gsl_multiset_memcpy(m_gm, c.m_gm);
    }

    multiset(multiset &&c)
    {
        m_gm = c.m_gm;
        c.m_gm = nullptr;
    }

    ~multiset()
    {
        gsl_multiset_free(m_gm);
    }

    multiset &operator=(const multiset &c)
    {
        eigen_assert((c.n() == n()) && (c.k() == k()));
        gsl_multiset_memcpy(m_gm, c.m_gm);
        return *this;
    }

    multiset &operator=(const multiset &&c)
    {
        memcpy(m_gm, c.m_gm, sizeof(gsl_multiset));
        memset(c.m_gm, 0, sizeof(gsl_multiset));
        return *this;
    }

    index &operator[](const index i)
    {
        eigen_assert(i < k());
        return *(&gsl_multiset_data(m_gm)[i]);
    }

    index &operator()(const index i)
    {
        eigen_assert(i < k());
        return *(&gsl_multiset_data(m_gm)[i]);
    }

    void init_first()
    {
        gsl_multiset_init_first(m_gm);
    }

    void init_last()
    {
        gsl_multiset_init_last(m_gm);
    }

    index n() const
    {
        return gsl_multiset_n(m_gm);
    }

    index k() const
    {
        return gsl_multiset_k(m_gm);
    }

    bool next()
    {
        return bool(gsl_multiset_next(m_gm) == GSL_SUCCESS);
    }

    bool next(multiset &c) const
    {
        c = *this;
        return c.next();
    }

    bool prev()
    {
        return bool(gsl_multiset_prev(m_gm) == GSL_SUCCESS);
    }

    bool prev(multiset &c) const
    {
        c = *this;
        return c.prev();
    }

    bool valid() const
    {
        return bool(gsl_multiset_valid(m_gm) == GSL_SUCCESS);
    }

  private:
    gsl_multiset *m_gm;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_MULTISET__ */
