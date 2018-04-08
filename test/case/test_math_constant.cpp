#include <catch.hpp>
#include <common/common.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math/constant.h>

using namespace iexp;

TEST_CASE("error")
{
    bool ok = false;

    try {
        GSL_ERROR("test", GSL_ERANGE);
    } catch (std::exception &e) {
        if (strncmp(e.what(),
                    "output range error:test",
                    sizeof("output range error:test") - 1) == 0) {
            ok = true;
        }
    }
    REQUIRE(ok);
}

TEST_CASE("type")
{
    REQUIRE(IS_MATRIX(Matrix3f));
    REQUIRE(!IS_MATRIX(Array3f));

    REQUIRE(!IS_ARRAY(Matrix3f));
    REQUIRE(IS_ARRAY(Array3f));

    Matrix<char, 2, 3, RowMajor, 2, 3> m1;
    REQUIRE((TYPE_IS(TP1(decltype(m1)), char)));
    REQUIRE(TP2(decltype(m1)) == 2);
    REQUIRE(TP3(decltype(m1)) == 3);
    REQUIRE((TP4(decltype(m1)) == RowMajor));
    REQUIRE(TP5(decltype(m1)) == 2);
    REQUIRE(TP6(decltype(m1)) == 3);
    REQUIRE(!IS_DYNAMIC(decltype(m1)));

    dense_derive<decltype(m1)>::type d1;
    REQUIRE((TYPE_IS(TP1(decltype(d1)), char)));
    REQUIRE(TP2(decltype(d1)) == 2);
    REQUIRE(TP3(decltype(d1)) == 3);
    REQUIRE((TP4(decltype(d1)) == RowMajor));
    REQUIRE(TP5(decltype(d1)) == 2);
    REQUIRE(TP6(decltype(d1)) == 3);

    dense_derive<decltype(m1 + m1)>::type d2;
    REQUIRE((TYPE_IS(TP1(decltype(d2)), char)));
    REQUIRE(TP2(decltype(d2)) == 2);
    REQUIRE(TP3(decltype(d2)) == 3);
    REQUIRE((TP4(decltype(d2)) == RowMajor));
    REQUIRE(TP5(decltype(d2)) == 2);
    REQUIRE(TP6(decltype(d2)) == 3);

    Array<char, Dynamic, Dynamic, 0, 2, 3> a1;
    REQUIRE((TYPE_IS(TP1(decltype(a1)), char)));
    REQUIRE(TP2(decltype(a1)) == Dynamic);
    REQUIRE(TP3(decltype(a1)) == Dynamic);
    REQUIRE((TP4(decltype(a1)) == ColMajor));
    REQUIRE(TP5(decltype(a1)) == 2);
    REQUIRE(TP6(decltype(a1)) == 3);
    REQUIRE(!IS_DYNAMIC(decltype(a1)));

    dense_derive<decltype(a1)>::type e1;
    REQUIRE((TYPE_IS(TP1(decltype(e1)), char)));
    REQUIRE(TP2(decltype(e1)) == Dynamic);
    REQUIRE(TP3(decltype(e1)) == Dynamic);
    REQUIRE((TP4(decltype(e1)) == ColMajor));
    REQUIRE(TP5(decltype(e1)) == 2);
    REQUIRE(TP6(decltype(e1)) == 3);

    dense_derive<decltype(a1 + a1)>::type e2;
    REQUIRE((TYPE_IS(TP1(decltype(e2)), char)));
    REQUIRE(TP2(decltype(e2)) == Dynamic);
    REQUIRE(TP3(decltype(e2)) == Dynamic);
    REQUIRE((TP4(decltype(e2)) == ColMajor));
    REQUIRE(TP5(decltype(e2)) == 2);
    REQUIRE(TP6(decltype(e2)) == 3);

    // a colmajor + rowmajoe
    Array<char, Dynamic, Dynamic, RowMajor, 2, 3> a2;
    dense_derive<decltype(a2 + a1)>::type e3;
    REQUIRE((TP4(decltype(e3)) == ColMajor));
    REQUIRE(!IS_DYNAMIC(decltype(a2)));

    {
        Matrix<int, Dynamic, 1> m;
        REQUIRE(IS_DYNAMIC(decltype(m)));
        Matrix<int, 1, Dynamic> n;
        REQUIRE(IS_DYNAMIC(decltype(n)));
    }

    REQUIRE(IS_INTEGER(int));
}

TEST_CASE("buf")
{
    buf<Matrix<int, 1, 2>> b1(2);
    REQUIRE(b1.size() == 2);
    REQUIRE(b1.is_static());

    buf<Matrix<int, 1, 2>> b2(3);
    REQUIRE(b2.size() == 3);
    REQUIRE(!b2.is_static());

    buf<Matrix<int, Dynamic, 2>> b3(1);
    REQUIRE(b3.size() == 1);
    REQUIRE(b3.is_static());

    buf<Matrix<int, Dynamic, 2>> b4(2);
    REQUIRE(b4.size() == 2);
    REQUIRE(!b4.is_static());
}

TEST_CASE("math_constant")
{
    REQUIRE(IEXP_E == M_E);
    REQUIRE(IEXP_LOG2E == M_LOG2E);
    REQUIRE(IEXP_LOG10E == M_LOG10E);
    REQUIRE(IEXP_SQRT2 == M_SQRT2);
    REQUIRE(IEXP_SQRT1_2 == M_SQRT1_2);
    REQUIRE(IEXP_SQRT3 == M_SQRT3);
    REQUIRE(IEXP_PI == M_PI);
    REQUIRE(IEXP_PI_2 == M_PI_2);
    REQUIRE(IEXP_PI_4 == M_PI_4);
    REQUIRE(IEXP_SQRTPI == M_SQRTPI);
    REQUIRE(IEXP_2_SQRTPI == M_2_SQRTPI);
    REQUIRE(IEXP_1_PI == M_1_PI);
    REQUIRE(IEXP_2_PI == M_2_PI);
    REQUIRE(IEXP_LN10 == M_LN10);
    REQUIRE(IEXP_LN2 == M_LN2);
    REQUIRE(IEXP_LNPI == M_LNPI);
    REQUIRE(IEXP_EULER == M_EULER);
}
