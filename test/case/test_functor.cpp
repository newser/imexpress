#include <catch.hpp>
#include <common/common.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math/constant.h>

using namespace iexp;

template <class T>
class test_f_foreach : public functor_foreach<test_f_foreach<T>, T, int>
{
  public:
    test_f_foreach(const T &x)
        : functor_foreach<test_f_foreach<T>, T, int>(x)
    {
    }

    int foreach_impl(int i) const
    {
        return 2 * i;
    }
};

template <class T, class U>
class test_f_foreach_e
    : public functor_foreach_e<test_f_foreach_e<T, U>, T, U, int>
{
  public:
    test_f_foreach_e(const T &x, U &e)
        : functor_foreach_e<test_f_foreach_e<T, U>, T, U, int>(x, e)
    {
    }

    int foreach_e_impl(int x, char &e) const
    {
        e = x;
        return 2 * x;
    }
};

TEST_CASE("functor_foreach")
{
    Matrix<int, 3, 2> m;
    m << 1, 2, 3, 4, 5, 6;
    test_f_foreach<decltype(m)> ff(m);
    REQUIRE(ff(0, 0) == 2);
    REQUIRE(ff(0, 1) == 4);
    REQUIRE(ff(2, 0) == 10);
    REQUIRE(ff(2, 1) == 12);

    Matrix<char, 3, 2> e;
    e.setZero();
    test_f_foreach_e<decltype(m), decltype(e)> ffe(m, e);
    REQUIRE(ffe(0, 0) == 2);
    REQUIRE(ffe(0, 1) == 4);
    REQUIRE(ffe(2, 0) == 10);
    REQUIRE(ffe(2, 1) == 12);
    REQUIRE(e(0, 0) == 1);
    REQUIRE(e(0, 1) == 2);
    REQUIRE(e(2, 0) == 5);
    REQUIRE(e(2, 1) == 6);

    test_f_foreach<decltype(m + m)> ff2(m + m);
    REQUIRE(ff2(0, 0) == 4);
    REQUIRE(ff2(0, 1) == 8);
    REQUIRE(ff2(2, 0) == 20);
    REQUIRE(ff2(2, 1) == 24);

    test_f_foreach_e<decltype(m + m), decltype(e)> ff2e(m + m, e);
    REQUIRE(ff2e(0, 0) == 4);
    REQUIRE(ff2e(0, 1) == 8);
    REQUIRE(ff2e(2, 0) == 20);
    REQUIRE(ff2e(2, 1) == 24);
    REQUIRE(e(0, 0) == 2);
    REQUIRE(e(0, 1) == 4);
    REQUIRE(e(2, 0) == 10);
    REQUIRE(e(2, 1) == 12);
}

template <class T>
class test_f_2d : public functor_m2vnum_2d<test_f_2d<T>, T, int>
{
  public:
    test_f_2d(const T &x)
        : functor_m2vnum_2d<test_f_2d<T>, T, int>(x)
    {
    }

    int m2vnum_impl(int x0, int x1) const
    {
        return 2 * x0 + x1;
    }
};

template <class T, class U>
class test_f_2d_e : public functor_m2vnum_2d_e<test_f_2d_e<T, U>, T, U, int>
{
  public:
    test_f_2d_e(const T &x, U &e)
        : functor_m2vnum_2d_e<test_f_2d_e<T, U>, T, U, int>(x, e)
    {
    }

    int m2vnum_e_impl(int x0, int x1, char &e) const
    {
        e = (char)x1;
        return 2 * x0 + x1;
    }
};

TEST_CASE("functor_2d")
{
    Matrix<int, 2, 3> m;
    m << 1, 2, 3, 4, 5, 6;

    test_f_2d<decltype(m)> ff(m);
    REQUIRE(ff(0) == 6);
    REQUIRE(ff(1) == 9);
    REQUIRE(ff(2) == 12);

    Matrix<char, 1, 3> e;
    e.setZero();
    test_f_2d_e<decltype(m + m), decltype(e)> ffe(m + m, e);
    REQUIRE(ffe(0) == 12);
    REQUIRE(ffe(1) == 18);
    REQUIRE(ffe(2) == 24);
    REQUIRE(e(0) == 8);
    REQUIRE(e(1) == 10);
    REQUIRE(e(2) == 12);

    ////////////////////////

    Matrix<int, 3, 2, RowMajor> m2;
    m2 << 1, 2, 3, 4, 5, 6;

    test_f_2d<decltype(m2)> ff2(m2);
    REQUIRE(ff2(0) == 4);
    REQUIRE(ff2(1) == 10);
    REQUIRE(ff2(2) == 16);

    e.setZero();
    test_f_2d_e<decltype(m2 + m2), decltype(e)> ffe2(m2 + m2, e);
    REQUIRE(ffe2(0) == 8);
    REQUIRE(ffe2(1) == 20);
    REQUIRE(ffe2(2) == 32);
    REQUIRE(e(0) == 4);
    REQUIRE(e(1) == 8);
    REQUIRE(e(2) == 12);
}

template <class T>
class test_f_3d : public functor_m2vnum_3d<test_f_3d<T>, T, int>
{
  public:
    test_f_3d(const T &x)
        : functor_m2vnum_3d<test_f_3d<T>, T, int>(x)
    {
    }

    int m2vnum_impl(int x0, int x1, int x2) const
    {
        return x0 + x1 + x2;
    }
};

template <class T, class U>
class test_f_3d_e : public functor_m2vnum_3d_e<test_f_3d_e<T, U>, T, U, int>
{
  public:
    test_f_3d_e(const T &x, U &e)
        : functor_m2vnum_3d_e<test_f_3d_e<T, U>, T, U, int>(x, e)
    {
    }

    int m2vnum_e_impl(int x0, int x1, int x2, char &e) const
    {
        e = (char)x2;
        return x0 + x1 + x2;
    }
};

TEST_CASE("functor_3d")
{
    Matrix<int, 3, 2> m;
    m << 1, 2, 3, 4, 5, 6;

    test_f_3d<decltype(m)> ff(m);
    REQUIRE(ff(0) == 9);
    REQUIRE(ff(1) == 12);

    Matrix<char, 1, 2> e;
    e.setZero();
    test_f_3d_e<decltype(m + m), decltype(e)> ffe(m + m, e);
    REQUIRE(ffe(0) == 18);
    REQUIRE(ffe(1) == 24);
    REQUIRE(e(0) == 10);
    REQUIRE(e(1) == 12);

    ////////////////////////

    Matrix<int, 2, 3, RowMajor> m2;
    m2 << 1, 2, 3, 4, 5, 6;

    test_f_3d<decltype(m2)> ff2(m2);
    REQUIRE(ff2(0) == 6);
    REQUIRE(ff2(1) == 15);

    e.setZero();
    test_f_3d_e<decltype(m2 + m2), decltype(e)> ffe2(m2 + m2, e);
    REQUIRE(ffe2(0) == 12);
    REQUIRE(ffe2(1) == 30);
    REQUIRE(e(0) == 6);
    REQUIRE(e(1) == 12);
}

template <class T>
class test_f_m2vd : public functor_m2vdim_va<test_f_m2vd<T>, T, short>
{
  public:
    test_f_m2vd(const T &x)
        : functor_m2vdim_va<test_f_m2vd<T>, T, short>(x)
    {
    }

    void m2vdim_va_impl(const unsigned char *data,
                        int num,
                        int dim,
                        short *result) const
    {
        m_data = data;
        m_num = num;
        m_dim = dim;
        m_result = result;
    }

    mutable const unsigned char *m_data;
    mutable int m_num, m_dim;
    mutable short *m_result;
};

TEST_CASE("functor_m2vdim")
{
    Matrix<unsigned char, 3, 2> m;
    m << 1, 2, 3, 4, 5, 6;

    test_f_m2vd<decltype(m)> fmd(m);
    REQUIRE(fmd.m_num == 2);
    REQUIRE(fmd.m_dim == 3);

    test_f_m2vd<decltype(m + m)> fmd2(m + m);
    REQUIRE(fmd2.m_num == 2);
    REQUIRE(fmd2.m_dim == 3);

    ////////////////////////

    Matrix<unsigned char, 3, 2, RowMajor> m2;
    m2 << 1, 2, 3, 4, 5, 6;

    test_f_m2vd<decltype(m2)> fmd3(m2);
    REQUIRE(fmd3.m_num == 3);
    REQUIRE(fmd3.m_dim == 2);

    test_f_m2vd<decltype(m2 + m2)> fmd4(m2 + m2);
    REQUIRE(fmd4.m_num == 3);
    REQUIRE(fmd4.m_dim == 2);
}
