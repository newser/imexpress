#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/clausen.h>
#include <test_util.h>

TEST_CASE("test_clausen")
{
    iexp::ArrayXd m(5), m2(5), e(5);
    gsl_sf_result r;

    m << IEXP_PI / 20.0, IEXP_PI / 6.0, IEXP_PI / 3.0,
        2.0 * IEXP_PI + IEXP_PI / 3.0, 100.0 * IEXP_PI + IEXP_PI / 3.0;
    m2 = iexp::sf::clausen(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.4478882448133546));
    REQUIRE(__D_EQ9(m2(1, 0), 0.8643791310538927));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0149416064096535));
    REQUIRE(__D_EQ9(m2(3, 0), 1.0149416064096535));
    REQUIRE(__D_EQ9(m2(4, 0), 1.0149416064096535));

    __draw_1d(clausen, 0, 10, 100, 0, 0);

    m2 = iexp::sf::clausen(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.4478882448133546));
    REQUIRE(__D_EQ9(m2(1, 0), 0.8643791310538927));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0149416064096535));
    REQUIRE(__D_EQ9(m2(3, 0), 1.0149416064096535));
    REQUIRE(__D_EQ9(m2(4, 0), 1.0149416064096535));
    gsl_sf_clausen_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_clausen_e(m2(4, 0), &r);
    REQUIRE(__D_EQ9(e(4, 0), r.err));
}
