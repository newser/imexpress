#include <../test/test_util.h>
#include <catch.hpp>
#include <histogram/hist.h>

using namespace iexp;

void test_empty(hist &h)
{
    REQUIRE(h[0] == 0);
    REQUIRE(h[3] == 0);
    double low, up;
    h.range(0, low, up);
    REQUIRE((low == 1.0 && up == 2.0));
    h.range(3, low, up);
    REQUIRE((low == 4.0 && up == 5.0));
    REQUIRE(h.max() == 5.0);
    REQUIRE(h.min() == 1.0);
    REQUIRE(h.size() == 4);
    REQUIRE(h.find(1.5) == 0);
    REQUIRE(h.find(3.5) == 2);
    REQUIRE(h.find(4.5) == 3);
    REQUIRE(h.max_val() == 0);
    REQUIRE(h.min_val() == 0);
    REQUIRE(h.mean() == 0);
    REQUIRE(h.std() == 0);
    REQUIRE(h.sum() == 0);
}

TEST_CASE("hist_empty")
{
    double r[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    hist h(r, sizeof(r) / sizeof(r[0]));

    test_empty(h);

    hist h2(h);
    test_empty(h2);

    hist h3(std::move(h2));
    test_empty(h3);

    hist h4 = h3;
    test_empty(h3);

    hist h5 = std::move(h4);
    test_empty(h5);
    REQUIRE(h5 == h);

    h += h5;
    test_empty(h);
    h += 0;
    test_empty(h);
    h -= h5;
    test_empty(h);
    h -= 0;
    test_empty(h);
    h *= h5;
    test_empty(h);
    h *= 6.0;
    test_empty(h);
    h /= 6.0;
    test_empty(h);
}

void test_val(hist &h)
{
    REQUIRE(h[0] == 2);
    REQUIRE(h[3] == 1);
    double low, up;
    h.range(0, low, up);
    REQUIRE((low == 0.0 && up == 1.0));
    h.range(3, low, up);
    REQUIRE((low == 3.0 && up == 4.0));
    REQUIRE(h.max() == 10.0);
    REQUIRE(h.min() == 0.0);
    REQUIRE(h.size() == 10);
    REQUIRE(h.find(1.5) == 1);
    REQUIRE(h.find(3.5) == 3);
    REQUIRE(h.find(9.5) == 9);
    REQUIRE(h.max_val() == 2.0);
    REQUIRE(h.max_idx() == 0);
    REQUIRE(h.min_val() == -1.0);
    REQUIRE(h.min_idx() == 9);
    REQUIRE(__D_EQ7(h.mean(), 2.6428571428571428));
    REQUIRE(__D_EQ7(h.std(), 1.8070158058105026));
    REQUIRE(__D_EQ7(h.sum(), 6));
}

TEST_CASE("hist_val")
{
    double r[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    hist h(r, sizeof(r) / sizeof(r[0]));

    h << 0.5 << 1.5 << 2.5 << 3.5 << 4.5 << 5.5;
    h << 0.5;
    h.add(9.9, -1);
    test_val(h);

    hist h2(h);
    test_val(h2);

    hist h3(std::move(h2));
    test_val(h3);

    hist h4 = h3;
    test_val(h3);

    hist h5 = std::move(h4);
    test_val(h5);
    REQUIRE(h5 == h);

    h5.reset();
    h += h5;
    test_val(h);
    h += 0;
    test_val(h);
    h -= h5;
    test_val(h);
    h -= 0;
    test_val(h);

    h5 << 0.1 << 1.1 << 2.1 << 3.1 << 4.1 << 5.1 << 6.1 << 7.1 << 8.1;
    h5.add(9.1, 1.0);
    h *= h5;
    test_val(h);
    h *= 1.0;
    test_val(h);
    h /= h5;
    test_val(h);
    h /= 1.0;
    test_val(h);
}
