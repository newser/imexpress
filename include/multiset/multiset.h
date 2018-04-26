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
    enum class init
    {
        NONE,
        FIRST,
        LAST,
    };

    multiset(size_t n, size_t k, init i = init::NONE)
    {
        m_gm = gsl_multiset_alloc(n, k);
        IEXP_NOT_NULLPTR(m_gm);
        if (i == init::FIRST) {
            gsl_multiset_init_first(m_gm);
        } else if (i == init::LAST) {
            gsl_multiset_init_last(m_gm);
        }
    }

    template <typename T>
    multiset(size_t n, const DenseBase<T> &v)
    {
        static_assert(IS_INTEGER(typename T::Scalar), "must be integer type");

        eigen_assert(n >= v.size());
        m_gm = gsl_multiset_alloc(n, v.size());
        IEXP_NOT_NULLPTR(m_gm);

        typename type_eval<T>::type e_v(v.eval());
        copy<size_t, const typename T::Scalar>(gsl_multiset_data(m_gm),
                                               e_v.data(),
                                               e_v.size());
        gsl_multiset_valid(m_gm);
    }

    multiset(const multiset &rhs)
    {
        m_gm = gsl_multiset_alloc(rhs.n(), rhs.k());
        IEXP_NOT_NULLPTR(m_gm);
        gsl_multiset_memcpy(m_gm, rhs.m_gm);
    }

    multiset(multiset &&c)
    {
        m_gm = c.m_gm;
        c.m_gm = nullptr;
    }

    ~multiset()
    {
        if (m_gm != nullptr) {
            gsl_multiset_free(m_gm);
        }
    }

    multiset &operator=(const multiset &rhs)
    {
        eigen_assert((rhs.n() == n()) && (rhs.k() == k()));
        gsl_multiset_memcpy(m_gm, rhs.m_gm);
        return *this;
    }

    multiset &operator=(multiset &&rhs)
    {
        if (m_gm != nullptr) {
            gsl_multiset_free(m_gm);
        }
        m_gm = rhs.m_gm;
        rhs.m_gm = nullptr;
        return *this;
    }

    template <typename T>
    multiset &operator=(const DenseBase<T> &v)
    {
        static_assert(IS_INTEGER(typename T::Scalar), "must be integer type");

        size_t n = gsl_multiset_n(m_gm);
        if (n < v.size()) {
            throw std::invalid_argument("invalid input size");
        }

        size_t k = gsl_multiset_k(m_gm);
        if (k == v.size()) {
            typename type_eval<T>::type e_v(v.eval());
            copy<size_t, const typename T::Scalar>(gsl_multiset_data(m_gm),
                                                   e_v.data(),
                                                   k);
            gsl_multiset_valid(m_gm);
        } else {
            // rvalue assignment
            *this = multiset(n, v);
        }
        return *this;
    }

    size_t &operator[](size_t i) const
    {
        eigen_assert(i < k());
        return *(&gsl_multiset_data(m_gm)[i]);
    }

    size_t &operator()(size_t i) const
    {
        eigen_assert(i < k());
        return *(&gsl_multiset_data(m_gm)[i]);
    }

    multiset &init_first()
    {
        gsl_multiset_init_first(m_gm);
        return *this;
    }

    multiset &init_last()
    {
        gsl_multiset_init_last(m_gm);
        return *this;
    }

    size_t n() const
    {
        return gsl_multiset_n(m_gm);
    }

    size_t k() const
    {
        return gsl_multiset_k(m_gm);
    }

    bool next()
    {
        return bool(gsl_multiset_next(m_gm) == GSL_SUCCESS);
    }

    bool prev()
    {
        return bool(gsl_multiset_prev(m_gm) == GSL_SUCCESS);
    }

    bool valid() const
    {
        try {
            gsl_multiset_valid(m_gm);
            return true;
        } catch (...) {
            return false;
        }
    }

  private:
    gsl_multiset *m_gm;

    // to matrix
    template <typename _Scalar, int _Rows, int _Cols>
    class matrix_functor
    {
      public:
        matrix_functor(const multiset &ms)
            : m_ms(ms)
        {
        }

        _Scalar operator()(Index i, Index j) const
        {
            // always colume major
            return static_cast<_Scalar>(m_ms(i + j * _Rows));
        }

      private:
        const multiset &m_ms;
    };

    // to array
    template <typename _Scalar, int _Rows, int _Cols>
    class array_functor
    {
      public:
        array_functor(const multiset &ms)
            : m_ms(ms)
        {
        }

        _Scalar operator()(Index i, Index j) const
        {
            // always colume major
            return static_cast<_Scalar>(m_ms(i + j * _Rows));
        }

      private:
        const multiset &m_ms;
    };

  public:
    template <typename _Scalar, int _Rows, int _Cols>
    inline CwiseNullaryOp<matrix_functor<_Scalar, _Rows, _Cols>,
                          Matrix<_Scalar, _Rows, _Cols>>
    matrix()
    {
        static_assert(IS_INTEGER(_Scalar), "must be integer type");
        static_assert((_Rows != Dynamic) && (_Cols != Dynamic),
                      "can not be dynamic size");

        size_t k = gsl_multiset_k(m_gm);
        if (k != _Rows * _Cols) {
            throw std::invalid_argument("invalid matrix size");
        }

        using MatrixType = Matrix<_Scalar, _Rows, _Cols>;
        return MatrixType::NullaryExpr(_Rows,
                                       _Cols,
                                       matrix_functor<_Scalar, _Rows, _Cols>(
                                           *this));
    }

    template <typename _Scalar, int _Rows, int _Cols>
    inline CwiseNullaryOp<array_functor<_Scalar, _Rows, _Cols>,
                          Array<_Scalar, _Rows, _Cols>>
    array()
    {
        static_assert(IS_INTEGER(_Scalar), "must be integer type");
        static_assert((_Rows != Dynamic) && (_Cols != Dynamic),
                      "can not be dynamic size");

        size_t k = gsl_multiset_k(m_gm);
        if (k != _Rows * _Cols) {
            throw std::invalid_argument("invalid array size");
        }

        using ArrayType = Array<_Scalar, _Rows, _Cols>;
        return ArrayType::NullaryExpr(_Rows,
                                      _Cols,
                                      array_functor<_Scalar, _Rows, _Cols>(
                                          *this));
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// size_terface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_MULTISET__ */
