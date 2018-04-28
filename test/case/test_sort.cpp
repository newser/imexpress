#include <catch.hpp>
#include <iostream>
#include <sort/largest.h>
#include <sort/smallest.h>
#include <sort/sort.h>
#include <sort/sort_index.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_sort")
{
    ArrayXi i(3), i2(3), i3(3);
    i << 3, 2, 1;
    i3 << 4, 5, 6;
    i2 = sort(i.array());
    REQUIRE(i2(0) == 1);
    REQUIRE(i2(1) == 2);
    REQUIRE(i2(2) == 3);

    i2 = sort(i.array(), i3);
    REQUIRE(i2(0) == 1);
    REQUIRE(i2(1) == 2);
    REQUIRE(i2(2) == 3);
    REQUIRE(i3(0) == 6);
    REQUIRE(i3(1) == 5);
    REQUIRE(i3(2) == 4);

    {
        VectorXd v1(3), v2(3), v3(3);
        v1 << 3, 1, 2;
        v2 << 9, 9, 9;
        v3 = sort(v1 + v2);
        REQUIRE(v3(0) == 10);
        REQUIRE(v3(1) == 11);
        REQUIRE(v3(2) == 12);
    }

    ArrayXf f(3), f2(3), f3(3);
    f << 3, 2, 1;
    f3 << 4, 5, 6;
    f2 = sort(f.array());
    REQUIRE(f2(0) == 1);
    REQUIRE(f2(1) == 2);
    REQUIRE(f2(2) == 3);

    f2 = sort(f.array(), f3);
    REQUIRE(f2(0) == 1);
    REQUIRE(f2(1) == 2);
    REQUIRE(f2(2) == 3);
    REQUIRE(f3(0) == 6);
    REQUIRE(f3(1) == 5);
    REQUIRE(f3(2) == 4);

    ArrayXd d(3), d2(3), d3(3);
    d << 3, 2, 1;
    d3 << 4, 5, 6;
    d2 = sort(d.array());
    REQUIRE(d2(0) == 1);
    REQUIRE(d2(1) == 2);
    REQUIRE(d2(2) == 3);

    d2 = sort(d.array(), d3);
    REQUIRE(d2(0) == 1);
    REQUIRE(d2(1) == 2);
    REQUIRE(d2(2) == 3);
    REQUIRE(d3(0) == 6);
    REQUIRE(d3(1) == 5);
    REQUIRE(d3(2) == 4);

    {
        VectorXd v1(3), v2(3), v3(3), cc(3);
        v1 << 3, 1, 2;
        v2 << 9, 9, 9;
        cc << 0, 1, 2;
        v3 = sort(v1 + v2, cc);
        REQUIRE(v3(0) == 10);
        REQUIRE(v3(1) == 11);
        REQUIRE(v3(2) == 12);
        REQUIRE(cc(0) == 1);
        REQUIRE(cc(1) == 2);
        REQUIRE(cc(2) == 0);
    }
}

TEST_CASE("test_sort_inplace")
{
    ArrayXd d(3), d3(3);

    d << 3, 2, 1;
    sort_inplace(d);
    REQUIRE(d(0) == 1);
    REQUIRE(d(1) == 2);
    REQUIRE(d(2) == 3);

    d << 3, 2, 1;
    d3 << 4, 5, 6;
    sort_inplace(d, d3);
    REQUIRE(d(0) == 1);
    REQUIRE(d(1) == 2);
    REQUIRE(d(2) == 3);
    REQUIRE(d3(0) == 6);
    REQUIRE(d3(1) == 5);
    REQUIRE(d3(2) == 4);
}

TEST_CASE("test_sort_idx")
{
    ArrayXi i(3), i2(3), i3(3);
    i << 3, 2, 1;
    i2 = sort_idx<int>(i);
    REQUIRE(i2(0) == 2);
    REQUIRE(i2(1) == 1);
    REQUIRE(i2(2) == 0);
}

TEST_CASE("test_smallest")
{
    ArrayXi i(4), i2(4), i3(4);
    i << 3, 2, 1, 4;
    i2 = smallest(2, i);
    REQUIRE(i2.size() == 2);
    REQUIRE((((i2(0) == 1) && (i2(1) == 2)) || ((i2(0) == 2) && (i2(1) == 1))));

    // row form
    RowVector2i j = smallest<true>(2, i);
    REQUIRE(j.size() == 2);
    REQUIRE((((j(0) == 1) && (j(1) == 2)) || ((j(0) == 2) && (j(1) == 1))));

    i2 = smallest_index<int>(2, i);
    REQUIRE(i2.size() == 2);
    REQUIRE((((i2(0) == 1) && (i2(1) == 2)) || ((i2(0) == 2) && (i2(1) == 1))));

    // row form
    j = smallest_index<int, true>(2, i);
    REQUIRE(j.size() == 2);
    REQUIRE((((j(0) == 1) && (j(1) == 2)) || ((j(0) == 2) && (j(1) == 1))));
}

TEST_CASE("test_largest")
{
    ArrayXi i(4), i2(4), i3(4);
    i << 3, 2, 1, 4;
    i2 = largest(2, i);
    REQUIRE(i2.size() == 2);
    REQUIRE((((i2(0) == 3) && (i2(1) == 4)) || ((i2(0) == 4) && (i2(1) == 3))));

    // row form
    RowVector2i j = largest<true>(2, i);
    REQUIRE(j.size() == 2);
    REQUIRE((((j(0) == 3) && (j(1) == 4)) || ((j(0) == 4) && (j(1) == 3))));

    i2 = largest_index<int>(2, i);
    REQUIRE(i2.size() == 2);
    REQUIRE((((i2(0) == 0) && (i2(1) == 3)) || ((i2(0) == 3) && (i2(1) == 0))));

    // row form
    j = largest_index<int, true>(2, i);
    REQUIRE(j.size() == 2);
    REQUIRE((((i2(0) == 0) && (i2(1) == 3)) || ((i2(0) == 3) && (i2(1) == 0))));
}
