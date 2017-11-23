#include <catch.hpp>
#include <iostream>
#include <special/cylindrical_bessel.h>
#include <test_util.h>

TEST_CASE("sf_cbessel_j")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // J0
    m << 0.1, 2.0, 100.0, 1.0e+10;
    m2 = iexp::sf::cbessel_j(0, m);
    REQUIRE(__D_EQ7(m2(0, 0), 0.99750156206604003230));
    REQUIRE(__D_EQ7(m2(1, 0), 0.22389077914123566805));
    REQUIRE(__D_EQ7(m2(2, 0), 0.019985850304223122424));
    REQUIRE(__D_EQ7(m2(3, 0), 2.1755917502468917269e-06));

    m2 = iexp::sf::cbessel_j(0, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), 0.99750156206604003230));
    REQUIRE(__D_EQ7(m2(3, 0), 2.1755917502468917269e-06));
    gsl_sf_bessel_j0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_j0_e(m(3, 0), &r);
    REQUIRE(__D_EQ7(e(3, 0), r.err));

    // J1
    m2 = iexp::sf::cbessel_j(1, m);
    REQUIRE(__D_EQ7(m2(0, 0), 0.04993752603624199756));
    REQUIRE(__D_EQ7(m2(1, 0), 0.57672480775687338720));
    REQUIRE(__D_EQ7(m2(2, 0), -0.07714535201411215803));
    REQUIRE(__D_EQ7(m2(3, 0), -7.676508175684157103e-06));

    m2 = iexp::sf::cbessel_j(1, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), 0.04993752603624199756));
    REQUIRE(__D_EQ7(m2(3, 0), -7.676508175684157103e-06));
    gsl_sf_bessel_j0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_j0_e(m(3, 0), &r);
    REQUIRE(__D_EQ7(e(3, 0), r.err));

    // Jn
    m.resize(2);
    m << 900.0, 15000.0;
    m2 = iexp::sf::cbessel_j(2, m);
    REQUIRE(__D_EQ7(m2(0, 0), -0.019974345269680646400));
    REQUIRE(__D_EQ7(m2(1, 0), -0.0020455820181216382666));

    m2 = iexp::sf::cbessel_j(2, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), -0.019974345269680646400));
    REQUIRE(__D_EQ7(m2(1, 0), -0.0020455820181216382666));
    gsl_sf_bessel_j0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_j0_e(m(1, 0), &r);
    REQUIRE(__D_EQ7(e(1, 0), r.err));
}

TEST_CASE("sf_cbessel_y")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // J0
    m << 0.1, 2.0, 256.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(0, m);
    REQUIRE(__D_EQ7(m2(0, 0), -1.5342386513503668441));
    REQUIRE(__D_EQ7(m2(1, 0), 0.5103756726497451196));
    REQUIRE(__D_EQ7(m2(2, 0), -0.03381290171792454909));
    REQUIRE(__D_EQ7(m2(3, 0), 3.657903190017678681e-06));

    m2 = iexp::sf::cbessel_y(0, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), -1.5342386513503668441));
    REQUIRE(__D_EQ7(m2(3, 0), 3.657903190017678681e-06));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(3, 0), &r);
    REQUIRE(__D_EQ7(e(3, 0), r.err));

    // J1
    m << 0.1, 2.0, 100.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(1, m);
    REQUIRE(__D_EQ7(m2(0, 0), -6.45895109470202698800));
    REQUIRE(__D_EQ7(m2(1, 0), -0.10703243154093754689));
    REQUIRE(__D_EQ7(m2(2, 0), -0.020372312002759793305));
    REQUIRE(__D_EQ7(m2(3, 0), 0.000011612249378370766284));

    m2 = iexp::sf::cbessel_y(1, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), -6.45895109470202698800));
    REQUIRE(__D_EQ7(m2(3, 0), 0.000011612249378370766284));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(3, 0), &r);
    REQUIRE(__D_EQ7(e(3, 0), r.err));

    // Jn
    m.resize(2);
    m << 100.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(100, m);
    REQUIRE(__D_EQ7(m2(0, 0), -0.16692141141757650654));
    REQUIRE(__D_EQ7(m2(1, 0), 3.657889671577715808e-06));

    m2 = iexp::sf::cbessel_y(100, m, e);
    REQUIRE(__D_EQ7(m2(0, 0), -0.16692141141757650654));
    REQUIRE(__D_EQ7(m2(1, 0), 3.657889671577715808e-06));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ7(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(1, 0), &r);
    REQUIRE(__D_EQ7(e(1, 0), r.err));
}
