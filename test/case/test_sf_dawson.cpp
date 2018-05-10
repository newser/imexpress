#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/dawson.h>
#include <special/dilog.h>
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

TEST_CASE("test_dilog")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    m << -3.0, -0.001, 0.1, 1100;
    m2 = iexp::sf::dilog(m);
    REQUIRE(__D_EQ9(m2(0, 0), -1.9393754207667089531));
    REQUIRE(__D_EQ9(m2(1, 0), -0.0009997501110486510834));
    REQUIRE(__D_EQ9(m2(2, 0), 0.1026177910993911));
    REQUIRE(__D_EQ9(m2(3, 0), -21.232504073931749553));

    __draw_1d(dilog, 0, 10, 100, 0, 0);

    m2 = iexp::sf::dilog(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -1.9393754207667089531));
    REQUIRE(__D_EQ9(m2(1, 0), -0.0009997501110486510834));
    REQUIRE(__D_EQ9(m2(2, 0), 0.1026177910993911));
    REQUIRE(__D_EQ9(m2(3, 0), -21.232504073931749553));
    gsl_sf_dawson_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_dawson_e(m2(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    iexp::ArrayXcd cm(4), cm2(4), ce(4);
    cm << std::polar(0.99999, M_PI / 2.0), std::polar(0.8, M_PI / 2.0),
        std::polar(0.5, M_PI / 2.0), std::polar(0.99, 3.0 * M_PI / 4.0);
    cm2 = iexp::sf::dilog(cm);
    REQUIRE(__D_EQ3(cm2(0, 0).real(), -0.20561329262779687646));
    REQUIRE(__D_EQ3(cm2(0, 0).imag(), 0.91595774018131512060));
    //    REQUIRE(__D_EQ3(cm2(3, 0).real(), -0.66210902664245926235));
    //    REQUIRE(__D_EQ3(cm2(3, 0).imag(), 0.51995305609998319025));
}
