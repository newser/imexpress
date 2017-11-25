#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/wigner_symbol.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_wigner3j")
{
    iexp::ArrayXi ja(4), jb(4), jc(4), ma(4), mb(4), mc(4);
    iexp::ArrayXd v(4), e(4);
    gsl_sf_result r;

    ja << 0, 1, 2, 4;
    jb << 1, 1, 4, 4;
    jc << 1, 2, 6, 8;
    ma << 0, 1, 0, 0;
    mb << 1, -1, 2, 0;
    mc << -1, 0, -2, 0;
    v = sf::wigner3j(ja, jb, jc, ma, mb, mc);
    REQUIRE(__D_EQ9(v(0, 0), sqrt(1.0 / 2.0)));
    REQUIRE(__D_EQ9(v(1, 0), sqrt(1.0 / 6.0)));
    REQUIRE(__D_EQ9(v(2, 0), sqrt(8.0 / 105.0)));
    REQUIRE(__D_EQ9(v(3, 0), sqrt(2.0 / 35.0)));

    v = sf::wigner3j(ja, jb, jc, ma, mb, mc, e);
    REQUIRE(__D_EQ9(v(0, 0), sqrt(1.0 / 2.0)));
    REQUIRE(__D_EQ9(v(1, 0), sqrt(1.0 / 6.0)));
    REQUIRE(__D_EQ9(v(2, 0), sqrt(8.0 / 105.0)));
    REQUIRE(__D_EQ9(v(3, 0), sqrt(2.0 / 35.0)));
    gsl_sf_coupling_3j_e(0, 1, 1, 0, 1, -1, &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_coupling_3j_e(4, 4, 8, 0, 0, 0, &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));
}

TEST_CASE("test_wigner6j")
{
    iexp::ArrayXi ja(4), jb(4), jc(4), ma(4), mb(4), mc(4);
    iexp::ArrayXd v(4), e(4);
    gsl_sf_result r;

    ja << 2, 4, 4, 4;
    jb << 2, 4, 4, 4;
    jc << 4, 2, 2, 2;
    ma << 2, 4, 4, 2;
    mb << 2, 4, 4, 2;
    mc << 2, 4, 2, 2;
    v = sf::wigner6j(ja, jb, jc, ma, mb, mc);
    REQUIRE(__D_EQ9(v(0, 0), 1.0 / 6.0));
    REQUIRE(__D_EQ9(v(1, 0), -1.0 / 10.0));
    REQUIRE(__D_EQ9(v(2, 0), 1.0 / 6.0));
    REQUIRE(__D_EQ9(v(3, 0), -0.5 / sqrt(5.0)));

    v = sf::wigner6j(ja, jb, jc, ma, mb, mc, e);
    REQUIRE(__D_EQ9(v(0, 0), 1.0 / 6.0));
    REQUIRE(__D_EQ9(v(1, 0), -1.0 / 10.0));
    REQUIRE(__D_EQ9(v(2, 0), 1.0 / 6.0));
    REQUIRE(__D_EQ9(v(3, 0), -0.5 / sqrt(5.0)));
    gsl_sf_coupling_6j_e(2, 2, 4, 2, 2, 2, &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_coupling_6j_e(4, 4, 2, 2, 2, 2, &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));
}

TEST_CASE("test_wigner9j")
{
    iexp::ArrayXi ja(2), jb(2), jc(2), jd(2), je(2), jf(2), jg(2), jh(2), ji(2);
    iexp::ArrayXd v(2), e(2);
    gsl_sf_result r;

    ja << 4, 8;
    jb << 2, 4;
    jc << 4, 10;
    jd << 3, 7;
    je << 3, 3;
    jf << 2, 8;
    jg << 1, 1;
    jh << 1, 1;
    ji << 2, 2;
    v = sf::wigner9j(ja, jb, jc, jd, je, jf, jg, jh, ji);
    REQUIRE(__D_EQ9(v(0, 0), -sqrt(1.0 / 6.0) / 10.0));
    REQUIRE(__D_EQ9(v(1, 0), sqrt(7.0 / 3.0) / 60.0));

    v = sf::wigner9j(ja, jb, jc, jd, je, jf, jg, jh, ji, e);
    REQUIRE(__D_EQ9(v(0, 0), -sqrt(1.0 / 6.0) / 10.0));
    REQUIRE(__D_EQ9(v(1, 0), sqrt(7.0 / 3.0) / 60.0));
    gsl_sf_coupling_9j_e(4, 2, 4, 3, 3, 2, 1, 1, 2, &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_coupling_9j_e(8, 4, 10, 7, 3, 8, 1, 1, 2, &r);
    REQUIRE(__D_EQ9(e(1, 0), r.err));
}
