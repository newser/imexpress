#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/hydrogenicR.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_hydrogenicR1")
{
    iexp::ArrayXd z(2), r(2), s(2), e(2);
    gsl_sf_result ret;

    z = 3.0;
    r << 2.0, 10.0;
    s = sf::hydroR(z, r);
    REQUIRE(__D_EQ9(s(0, 0), 0.025759948256148471036));
    REQUIRE(__D_EQ(s(1, 0), 9.724727052062819704e-13));

    s = sf::hydroR(z, r, e);
    REQUIRE(__D_EQ9(s(0, 0), 0.025759948256148471036));
    REQUIRE(__D_EQ(s(1, 0), 9.724727052062819704e-13));
    gsl_sf_hydrogenicR_1_e(z(0, 0), r(0, 0), &ret);
    REQUIRE(__D_EQ9(e(0, 0), ret.err));
    gsl_sf_hydrogenicR_1_e(z(1, 0), r(1, 0), &ret);
    REQUIRE(__D_EQ9(e(1, 0), ret.err));
}

TEST_CASE("test_hydrogenicR")
{
    iexp::ArrayXi n(4), l(4);
    iexp::ArrayXd z(4), r(4), s(4), e(4);
    gsl_sf_result ret;

    n << 4, 4, 100, 100;
    l << 1, 0, 0, 10;
    z = 3.0;
    r << 0.0, 2.0, 2.0, 2.0;
    s = sf::hydroR(n, l, z, r);
    REQUIRE(__D_EQ9(s(0, 0), 0.0));
    REQUIRE(__D_EQ(s(1, 0), -0.03623182256981820062));
    REQUIRE(__D_EQ(s(2, 0), -0.00007938950980052281367));
    REQUIRE(__D_EQ(s(3, 0), 7.112823375353605977e-12));

    s = sf::hydroR(n, l, z, r, e);
    REQUIRE(__D_EQ9(s(0, 0), 0.0));
    REQUIRE(__D_EQ(s(1, 0), -0.03623182256981820062));
    REQUIRE(__D_EQ(s(2, 0), -0.00007938950980052281367));
    REQUIRE(__D_EQ(s(3, 0), 7.112823375353605977e-12));
    gsl_sf_hydrogenicR_e(n(0, 0), l(0, 0), z(0, 0), r(0, 0), &ret);
    REQUIRE(__D_EQ9(e(0, 0), ret.err));
    gsl_sf_hydrogenicR_e(n(3, 0), l(3, 0), z(3, 0), r(3, 0), &ret);
    REQUIRE(__D_EQ9(e(3, 0), ret.err));
}
