#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/dawson.h>
#include <test_util.h>

TEST_CASE("test_dawson")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    m << 1.0e-15, 0.5, 2.0, 1000.0;
    m2 = iexp::sf::dawson(m);
    REQUIRE(__D_EQ(m2(0, 0), 1.0e-15));
    REQUIRE(__D_EQ9(m2(1, 0), 0.4244363835020222959));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30134038892379196603));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0005000002500003750009));

    __draw_1d(dawson, 0, 10, 100, 0, 0);

    m2 = iexp::sf::dawson(m, e);
    REQUIRE(__D_EQ(m2(0, 0), 1.0e-15));
    REQUIRE(__D_EQ9(m2(1, 0), 0.4244363835020222959));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30134038892379196603));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0005000002500003750009));
    gsl_sf_dawson_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_dawson_e(m2(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));
}
