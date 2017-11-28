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

#ifndef __IEXP_COMBINATION__
#define __IEXP_COMBINATION__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <gsl/gsl_combination.h>

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class combination
{
  public:
    typedef size_t index;

    enum init
    {
        none,
        first,
        last,
    };

    combination(const index n, const index k, init i = none)
    {
        m_gc = gsl_combination_alloc(n, k);
        IEXP_NOT_NULL(m_gc);
        if (i == first) {
            gsl_combination_init_first(m_gc);
        } else if (i == last) {
            gsl_combination_init_last(m_gc);
        }
    }

    combination(const combination &c)
    {
        m_gc = gsl_combination_alloc(c.n(), c.k());
        IEXP_NOT_NULL(m_gc);
        gsl_combination_memcpy(m_gc, c.m_gc);
    }

    combination(combination &&c)
    {
        m_gc = c.m_gc;
        c.m_gc = NULL;
    }

    ~combination()
    {
        gsl_combination_free(m_gc);
    }

    combination &operator=(const combination &c)
    {
        eigen_assert((c.n() == n()) && (c.k() == k()));
        gsl_combination_memcpy(m_gc, c.m_gc);
        return *this;
    }

    combination &operator=(const combination &&c)
    {
        memcpy(m_gc, c.m_gc, sizeof(gsl_combination));
        memset(c.m_gc, 0, sizeof(gsl_combination));
        return *this;
    }

    index &operator[](const index i)
    {
        eigen_assert(i < k());
        return *(&gsl_combination_data(m_gc)[i]);
    }

    index &operator()(const index i)
    {
        eigen_assert(i < k());
        return *(&gsl_combination_data(m_gc)[i]);
    }

    void init_first()
    {
        gsl_combination_init_first(m_gc);
    }

    void init_last()
    {
        gsl_combination_init_last(m_gc);
    }

    index n() const
    {
        return gsl_combination_n(m_gc);
    }

    index k() const
    {
        return gsl_combination_k(m_gc);
    }

    bool next()
    {
        return bool(gsl_combination_next(m_gc) == GSL_SUCCESS);
    }

    bool next(combination &c) const
    {
        c = *this;
        return c.next();
    }

    bool prev()
    {
        return bool(gsl_combination_prev(m_gc) == GSL_SUCCESS);
    }

    bool prev(combination &c) const
    {
        c = *this;
        return c.prev();
    }

    bool valid() const
    {
        return bool(gsl_combination_valid(m_gc) == GSL_SUCCESS);
    }

  private:
    gsl_combination *m_gc;
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_COMBINATION__ */
