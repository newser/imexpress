#include <catch.hpp>
#include <combination/combination.h>
#include <iostream>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_combination")
{
    combination c(3, 2);
    REQUIRE(c.n() == 3);
    REQUIRE(c.k() == 2);
    c.init_first();
    REQUIRE(c.valid() == true);
    REQUIRE(c[0] == 0);
    REQUIRE(c[1] == 1);
    c[0] = 1;
    c[1] = 2;
    REQUIRE(c.valid() == true);

    combination c2(c);
    REQUIRE(c2.n() == 3);
    REQUIRE(c2.k() == 2);
    REQUIRE(c2.valid() == true);
    REQUIRE(c2(0) == 1);
    REQUIRE(c2(1) == 2);
    c2(0) = 0;
    c2(1) = 2;
    REQUIRE(c.valid() == true);

    combination c3 = c2;
    REQUIRE(c3.n() == 3);
    REQUIRE(c3.k() == 2);
    REQUIRE(c3.valid() == true);
    REQUIRE(c3(0) == 0);
    REQUIRE(c3(1) == 2);
    REQUIRE(c.valid() == true);

    combination c4(std::move(c3));
    REQUIRE(c4.n() == 3);
    REQUIRE(c4.k() == 2);
    REQUIRE(c4.valid() == true);
    REQUIRE(c4(0) == 0);
    REQUIRE(c4(1) == 2);

    c4.init_last();
    REQUIRE(c4(0) == 1);
    REQUIRE(c4(1) == 2);
    REQUIRE(c4.next() == false);
    REQUIRE(c4.prev() == true);
    REQUIRE(c4(0) == 0);
    REQUIRE(c4(1) == 2);
    REQUIRE(c4.prev() == true);
    REQUIRE(c4.prev() == false);
    REQUIRE(c4(0) == 0);
    REQUIRE(c4(1) == 1);
    REQUIRE(c4.next() == true);
    REQUIRE(c4.next() == true);
    REQUIRE(c4(0) == 1);
    REQUIRE(c4(1) == 2);
    REQUIRE(c4.next() == false);
}

TEST_CASE("test_combination2")
{
    Array<int, Dynamic, Dynamic, RowMajor> a(3, 2), a0(3, 2);
    a << 0, 1, 2, 3, 4, 5;
    combination c(7, a);
    REQUIRE(c[0] == 0);
    REQUIRE(c[5] == 5);

    a0.setZero();
    combination c0(7, a + a0);
    REQUIRE(c0[0] == 0);
    REQUIRE(c0[5] == 5);

    Matrix<int, 2, 3> m = c.matrix<int, 2, 3>();
    REQUIRE(m(0, 0) == 0);
    REQUIRE(m(1, 2) == 5);

    a.resize(2, 1);
    a(0, 0) = 1;
    a(1, 0) = 2;
    c = a;
    REQUIRE(c[0] == 1);
    REQUIRE(c[1] == 2);

    // same size
    a(1, 0) = 3;
    c = a;
    REQUIRE(c[1] == 3);

    // too large
    except_begin()
    {
        a.resize(9, 9);
        c = a;
    }
    except_str("invalid input size");

    // matrix

    Matrix<int, Dynamic, Dynamic> b(3, 2), b0(3, 2);
    b << 0, 3, 1, 4, 2, 5;
    combination c2(7, b);
    REQUIRE(c2[0] == 0);
    REQUIRE(c2[5] == 5);

    b0.setZero();
    combination c20(7, b + b0);
    REQUIRE(c20[0] == 0);
    REQUIRE(c20[5] == 5);

    Array<unsigned char, 2, 3> ar = c2.matrix<unsigned char, 2, 3>();
    REQUIRE(ar(0, 0) == 0);
    REQUIRE(ar(1, 2) == 5);

    b.resize(2, 1);
    b(0, 0) = 1;
    b(1, 0) = 2;
    c2 = b;
    REQUIRE(c2[0] == 1);
    REQUIRE(c2[1] == 2);

    b0.resize(2, 1);
    c20 = b + b0;
    REQUIRE(c20[0] == 1);
    REQUIRE(c20[1] == 2);

    // same size
    b(1, 0) = 3;
    c2 = b;
    REQUIRE(c2[1] == 3);

    // too large
    except_begin()
    {
        b.resize(9, 9);
        c2 = b;
    }
    except_str("invalid input size");
}
