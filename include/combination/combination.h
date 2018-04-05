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

#include <common/common.h>
#include <common/copy.h>

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
    enum class init
    {
        NONE,
        FIRST,
        LAST,
    };

    combination(size_t n, size_t k, init i = init::NONE)
    {
        m_gc = gsl_combination_alloc(n, k);
        IEXP_NOT_NULLPTR(m_gc);
        if (i == init::FIRST) {
            gsl_combination_init_first(m_gc);
        } else if (i == init::LAST) {
            gsl_combination_init_last(m_gc);
        }
    }

    template <typename T>
    combination(size_t n, const DenseBase<T> &v)
    {
        static_assert(IS_INTEGER(typename T::Scalar), "must be integer type");

        eigen_assert(n >= v.size());
        m_gc = gsl_combination_alloc(n, v.size());
        IEXP_NOT_NULLPTR(m_gc);

        typename type_eval<T>::type e_v(v.eval());
        copy<size_t, const typename T::Scalar>(gsl_combination_data(m_gc),
                                               e_v.data(),
                                               e_v.size());
        gsl_combination_valid(m_gc);
    }

    combination(const combination &rhs)
    {
        m_gc = gsl_combination_alloc(rhs.n(), rhs.k());
        IEXP_NOT_NULLPTR(m_gc);
        gsl_combination_memcpy(m_gc, rhs.m_gc);
    }

    combination(combination &&rhs)
    {
        m_gc = rhs.m_gc;
        rhs.m_gc = nullptr;
    }

    ~combination()
    {
        if (m_gc != nullptr) {
            gsl_combination_free(m_gc);
        }
    }

    combination &operator=(const combination &rhs)
    {
        eigen_assert((rhs.n() == n()) && (rhs.k() == k()));
        gsl_combination_memcpy(m_gc, rhs.m_gc);
        return *this;
    }

    combination &operator=(combination &&rhs)
    {
        if (m_gc != nullptr) {
            gsl_combination_free(m_gc);
        }
        m_gc = rhs.m_gc;
        rhs.m_gc = nullptr;
        return *this;
    }

    template <typename T>
    combination &operator=(const DenseBase<T> &v)
    {
        static_assert(IS_INTEGER(typename T::Scalar), "must be integer type");

        size_t n = gsl_combination_n(m_gc);
        if (n < v.size()) {
            throw std::invalid_argument("invalid input size");
        }

        size_t k = gsl_combination_k(m_gc);
        if (k == v.size()) {
            typename type_eval<T>::type e_v(v.eval());
            copy<size_t, const typename T::Scalar>(gsl_combination_data(m_gc),
                                                   e_v.data(),
                                                   k);
            gsl_combination_valid(m_gc);
        } else {
            // rvalue assignment
            *this = combination(n, v);
        }
        return *this;
    }

    size_t &operator[](size_t i) const
    {
        eigen_assert(i < k());
        return *(&gsl_combination_data(m_gc)[i]);
    }

    size_t &operator()(size_t i) const
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

    size_t n() const
    {
        return gsl_combination_n(m_gc);
    }

    size_t k() const
    {
        return gsl_combination_k(m_gc);
    }

    bool next()
    {
        return bool(gsl_combination_next(m_gc) == GSL_SUCCESS);
    }

    bool prev()
    {
        return bool(gsl_combination_prev(m_gc) == GSL_SUCCESS);
    }

    bool valid() const
    {
        try {
            gsl_combination_valid(m_gc);
            return true;
        } catch (...) {
            return false;
        }
    }

  private:
    gsl_combination *m_gc;

    // to matrix
    template <typename _Scalar, int _Rows, int _Cols>
    class matrix_functor
    {
      public:
        matrix_functor(const combination &c)
            : m_c(c)
        {
        }

        _Scalar operator()(Index i, Index j) const
        {
            // always colume major
            return static_cast<_Scalar>(m_c(i + j * _Rows));
        }

      private:
        const combination &m_c;
    };

    // to array
    template <typename _Scalar, int _Rows, int _Cols>
    class array_functor
    {
      public:
        array_functor(const combination &c)
            : m_c(c)
        {
        }

        _Scalar operator()(Index i, Index j) const
        {
            // always colume major
            return static_cast<_Scalar>(m_c(i + j * _Rows));
        }

      private:
        const combination &m_c;
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

        size_t k = gsl_combination_k(m_gc);
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

        size_t k = gsl_combination_k(m_gc);
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
// indexerface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_COMBINATION__ */
