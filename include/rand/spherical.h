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

#ifndef __IEXP_RAND_SPHERICAL__
#define __IEXP_RAND_SPHERICAL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/common.h>

#include <rand/rng.h>

#include <gsl/gsl_randist.h>

IEXP_NS_BEGIN

namespace rand {

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

class sph
{
    template <typename Derived>
    class rng_base
    {
      public:
        rng_base(size_t dim, rand::rng::type type, unsigned long seed)
            : m_dim(dim)
            , m_rng(type, seed)
        {
        }

        void seed(unsigned long seed)
        {
            m_rng.seed(seed);
        }

        void next(double x[])
        {
            static_cast<Derived *>(this)->next_impl(x);
        }

        size_t dim() const
        {
            return m_dim;
        }

      protected:
        size_t m_dim;
        rand::rng m_rng;
    };

  public:
    class rng : public rng_base<rng>
    {
      public:
        rng(size_t dim,
            rand::rng::type type = DEFAULT_RNG_TYPE,
            unsigned long seed = 0)
            : rng_base<rng>(dim, type, seed)
        {
        }

        void next_impl(double x[])
        {
            gsl_ran_dir_nd(m_rng.gsl(), m_dim, x);
        }
    };

    class rng2 : public rng_base<rng2>
    {
      public:
        rng2(rand::rng::type type = DEFAULT_RNG_TYPE, unsigned long seed = 0)
            : rng_base<rng2>(2, type, seed)
        {
        }

        void next_impl(double x[])
        {
            gsl_ran_dir_2d_trig_method(m_rng.gsl(), &x[0], &x[1]);
        }
    };

    class rng3 : public rng_base<rng3>
    {
      public:
        rng3(rand::rng::type type = DEFAULT_RNG_TYPE, unsigned long seed = 0)
            : rng_base<rng3>(3, type, seed)
        {
        }

        void next_impl(double x[])
        {
            gsl_ran_dir_3d(m_rng.gsl(), &x[0], &x[1], &x[2]);
        }
    };

    template <typename T>
    static inline auto fill(DenseBase<T> &x,
                            unsigned long seed = 0,
                            rand::rng::type type = DEFAULT_RNG_TYPE)
        -> decltype(x.derived())
    {
        int dim = x.IsRowMajor ? x.cols() : x.rows();
        if (dim == 2) {
            sph::rng2 r(type, seed);
            return fill(x, r);
        } else if (dim == 3) {
            sph::rng3 r(type, seed);
            return fill(x, r);
        } else {
            sph::rng r(dim, type, seed);
            return fill(x, r);
        }
    }

    template <typename T, typename U>
    static inline auto fill(DenseBase<T> &x, rng_base<U> &r)
        -> decltype(x.derived())
    {
        return fill_impl(x, r, TYPE_BOOL(TP4(T) == RowMajor)());
    }

  private:
    sph() = delete;

    template <typename T, typename U>
    static inline auto fill_impl(DenseBase<T> &x,
                                 rng_base<U> &r,
                                 std::true_type) -> decltype(x.derived())
    {
        // row major
        size_t dim = r.dim();
        eigen_assert(x.cols() == dim);

        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");
        double *data = x.derived().data();
        for (Index i = 0; i < x.rows(); ++i) {
            r.next(&data[i * dim]);
        }
        return x.derived();
    }

    template <typename T, typename U>
    static inline auto fill_impl(DenseBase<T> &x,
                                 rng_base<U> &r,
                                 std::false_type) -> decltype(x.derived())
    {
        // col major
        size_t dim = r.dim();
        eigen_assert(x.rows() == dim);

        static_assert(TYPE_IS(typename T::Scalar, double),
                      "only support double scalar");
        double *data = x.derived().data();
        for (Index i = 0; i < x.cols(); ++i) {
            r.next(&data[i * dim]);
        }
        return x.derived();
    }
};

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// indexerface declaration
////////////////////////////////////////////////////////////
}

IEXP_NS_END

#endif /* __IEXP_RAND_SPHERICAL__ */
