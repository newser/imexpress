#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/hydrogenicR.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_hydrogenicR1")
{
    iexp::ArrayXd s(2), e(2);
    gsl_sf_result ret;

    iexp::Array<double, 2, 2, RowMajor> z_r;
    z_r.col(0).setConstant(3.0);
    z_r.col(1) << 2.0, 10.0;
    s = sf::hydroR1(z_r);
    REQUIRE(__D_EQ9(s(0, 0), 0.025759948256148471036));
    REQUIRE(__D_EQ(s(1, 0), 9.724727052062819704e-13));

    s = sf::hydroR1(z_r, e);
    REQUIRE(__D_EQ9(s(0, 0), 0.025759948256148471036));
    REQUIRE(__D_EQ(s(1, 0), 9.724727052062819704e-13));
    gsl_sf_hydrogenicR_1_e(z_r(0, 0), z_r(0, 1), &ret);
    REQUIRE(__D_EQ9(e(0, 0), ret.err));
    gsl_sf_hydrogenicR_1_e(z_r(1, 0), z_r(1, 1), &ret);
    REQUIRE(__D_EQ9(e(1, 0), ret.err));
}

TEST_CASE("test_hydrogenicR")
{
    iexp::ArrayXd s(4), e(4);
    gsl_sf_result ret;

    iexp::Array<double, 4, 4, RowMajor> nlzr;
    nlzr.col(0) << 4, 4, 100, 100;
    nlzr.col(1) << 1, 0, 0, 10;
    nlzr.col(2).setConstant(3.0);
    nlzr.col(3) << 0.0, 2.0, 2.0, 2.0;
    s = sf::hydroR(nlzr);
    REQUIRE(__D_EQ9(s(0, 0), 0.0));
    REQUIRE(__D_EQ(s(1, 0), -0.03623182256981820062));
    REQUIRE(__D_EQ(s(2, 0), -0.00007938950980052281367));
    REQUIRE(__D_EQ(s(3, 0), 7.112823375353605977e-12));

    s = sf::hydroR(nlzr, e);
    REQUIRE(__D_EQ9(s(0, 0), 0.0));
    REQUIRE(__D_EQ(s(1, 0), -0.03623182256981820062));
    REQUIRE(__D_EQ(s(2, 0), -0.00007938950980052281367));
    REQUIRE(__D_EQ(s(3, 0), 7.112823375353605977e-12));
    gsl_sf_hydrogenicR_e(nlzr(0, 0), nlzr(0, 1), nlzr(0, 2), nlzr(0, 3), &ret);
    REQUIRE(__D_EQ9(e(0, 0), ret.err));
    gsl_sf_hydrogenicR_e(nlzr(3, 0), nlzr(3, 1), nlzr(3, 2), nlzr(3, 3), &ret);
    REQUIRE(__D_EQ9(e(3, 0), ret.err));
}
