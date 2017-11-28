#include <catch.hpp>
#include <combination/combination.h>
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
