#include <catch.hpp>
#include <iostream>
#include <multiset/multiset.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_multiset")
{
    multiset c(3, 2);
    REQUIRE(c.n() == 3);
    REQUIRE(c.k() == 2);
    c.init_first();
    REQUIRE(c.valid() == true);
    REQUIRE(c[0] == 0);
    REQUIRE(c[1] == 0);
    c[0] = 1;
    c[1] = 2;
    REQUIRE(c.valid() == true);

    multiset c2(c);
    REQUIRE(c2.n() == 3);
    REQUIRE(c2.k() == 2);
    REQUIRE(c2.valid() == true);
    REQUIRE(c2(0) == 1);
    REQUIRE(c2(1) == 2);
    c2(0) = 0;
    c2(1) = 2;
    REQUIRE(c.valid() == true);

    multiset c3 = c2;
    REQUIRE(c3.n() == 3);
    REQUIRE(c3.k() == 2);
    REQUIRE(c3.valid() == true);
    REQUIRE(c3(0) == 0);
    REQUIRE(c3(1) == 2);
    REQUIRE(c.valid() == true);

    multiset c4(std::move(c3));
    REQUIRE(c4.n() == 3);
    REQUIRE(c4.k() == 2);
    REQUIRE(c4.valid() == true);
    REQUIRE(c4(0) == 0);
    REQUIRE(c4(1) == 2);

    c4.init_last();
    REQUIRE(c4(0) == 2);
    REQUIRE(c4(1) == 2);
    REQUIRE(c4.next() == false);
    REQUIRE(c4.prev() == true);
    REQUIRE(c4(0) == 1);
    REQUIRE(c4(1) == 2);
    REQUIRE(c4.prev() == true);
    c4(0) = c4(1) = 0;
    REQUIRE(c4.prev() == false);
    REQUIRE(c4.next() == true);
    REQUIRE(c4(0) == 0);
    REQUIRE(c4(1) == 1);
}
