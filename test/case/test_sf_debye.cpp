#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/debye.h>
#include <test_util.h>

TEST_CASE("test_debye")
{
    iexp::ArrayXd m(3), m2(3), e(3);
    gsl_sf_result r;

    // 1
    m << 0.1, 1.0, 10.0;
    m2 = iexp::sf::debye1(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.975277750004723276));
    REQUIRE(__D_EQ9(m2(1, 0), 0.777504634112248239));
    REQUIRE(__D_EQ9(m2(2, 0), 0.164443465679946027));

    m2 = iexp::sf::debye1(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.975277750004723276));
    REQUIRE(__D_EQ9(m2(1, 0), 0.777504634112248239));
    REQUIRE(__D_EQ9(m2(2, 0), 0.164443465679946027));
    gsl_sf_debye_1_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_1_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // 2
    m2 = iexp::sf::debye2(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.967083287045302664));
    REQUIRE(__D_EQ9(m2(1, 0), 0.70787847562782924));
    REQUIRE(__D_EQ9(m2(2, 0), 0.0479714980201218708));

    m2 = iexp::sf::debye2(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.967083287045302664));
    REQUIRE(__D_EQ9(m2(1, 0), 0.70787847562782924));
    REQUIRE(__D_EQ9(m2(2, 0), 0.0479714980201218708));
    gsl_sf_debye_2_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_2_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // 3
    m2 = iexp::sf::debye3(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.962999940487211048));
    REQUIRE(__D_EQ9(m2(1, 0), 0.674415564077814667));
    REQUIRE(__D_EQ9(m2(2, 0), 0.0192957656903454886));

    m2 = iexp::sf::debye3(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.962999940487211048));
    REQUIRE(__D_EQ9(m2(1, 0), 0.674415564077814667));
    REQUIRE(__D_EQ9(m2(2, 0), 0.0192957656903454886));
    gsl_sf_debye_3_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_3_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // 4
    m << 0.1, 1.0, 10.0;
    m2 = iexp::sf::debye4(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.960555486124335944));
    REQUIRE(__D_EQ9(m2(1, 0), 0.654874068886737049));
    REQUIRE(__D_EQ9(m2(2, 0), 0.00967367556027115896));

    m2 = iexp::sf::debye4(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.960555486124335944));
    REQUIRE(__D_EQ9(m2(1, 0), 0.654874068886737049));
    REQUIRE(__D_EQ9(m2(2, 0), 0.00967367556027115896));
    gsl_sf_debye_4_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_4_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // 5
    m2 = iexp::sf::debye5(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.95892849428310568745));
    REQUIRE(__D_EQ9(m2(1, 0), 0.6421002580217790246));
    REQUIRE(__D_EQ9(m2(2, 0), 0.005701535852992908538));

    m2 = iexp::sf::debye5(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.95892849428310568745));
    REQUIRE(__D_EQ9(m2(1, 0), 0.6421002580217790246));
    REQUIRE(__D_EQ9(m2(2, 0), 0.005701535852992908538));
    gsl_sf_debye_5_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_5_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // 6
    m2 = iexp::sf::debye6(m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.95776777382605465878));
    REQUIRE(__D_EQ9(m2(1, 0), 0.63311142583495107588));
    REQUIRE(__D_EQ9(m2(2, 0), 3.7938493294615955279e-3));

    m2 = iexp::sf::debye6(m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.95776777382605465878));
    REQUIRE(__D_EQ9(m2(1, 0), 0.63311142583495107588));
    REQUIRE(__D_EQ9(m2(2, 0), 3.7938493294615955279e-3));
    gsl_sf_debye_6_e(m2(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_debye_6_e(m2(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));
}
