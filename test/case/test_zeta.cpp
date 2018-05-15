#include <catch.hpp>
#include <iostream>
#include <special/zeta.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_zeta")
{
    Matrix2i p;
    p << -5, -1, 5, 31;

    Matrix2d r = sf::zeta(p), e;
    REQUIRE(__D_EQ9(r(0, 0), -0.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -1.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 1.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 1.0000000004656629065033784));

    r = sf::zeta(p, e);
    REQUIRE(__D_EQ9(r(0, 0), -0.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -1.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 1.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 1.0000000004656629065033784));

    gsl_sf_result re;
    gsl_sf_zeta_int_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_zeta_int_e(31, &re);
    REQUIRE(re.err == e(1, 1));

    // double
    Matrix2d p2;
    p2 << -5, -1, 5, 31;

    r = sf::zeta(p2);
    REQUIRE(__D_EQ9(r(0, 0), -0.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -1.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 1.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 1.0000000004656629065033784));

    r = sf::zeta(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), -0.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -1.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 1.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 1.0000000004656629065033784));

    gsl_sf_zeta_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_zeta_e(31, &re);
    REQUIRE(re.err == e(1, 1));
}

TEST_CASE("test_zetam1")
{
    Matrix2i p;
    p << -5, -1, 5, 31;

    Matrix2d r = sf::zetam1(p), e;
    REQUIRE(__D_EQ9(r(0, 0), -1.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -13.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 0.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 0.0000000004656629065033784));

    r = sf::zetam1(p, e);
    REQUIRE(__D_EQ9(r(0, 0), -1.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -13.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 0.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 0.0000000004656629065033784));

    gsl_sf_result re;
    gsl_sf_zetam1_int_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_zetam1_int_e(31, &re);
    REQUIRE(re.err == e(1, 1));

    // double
    Matrix2d p2;
    p2 << -5, -1, 5, 31;

    r = sf::zetam1(p2);
    REQUIRE(__D_EQ9(r(0, 0), -1.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -13.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 0.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 0.0000000004656629065033784));

    r = sf::zetam1(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), -1.003968253968253968253968));
    REQUIRE(__D_EQ9(r(0, 1), -13.0 / 12.0));
    REQUIRE(__D_EQ9(r(1, 0), 0.0369277551433699263313655));
    REQUIRE(__D_EQ9(r(1, 1), 0.0000000004656629065033784));

    gsl_sf_zetam1_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_zetam1_e(31, &re);
    REQUIRE(re.err == e(1, 1));
}

TEST_CASE("test_eta")
{
    Matrix2i p;
    p << -5, -1, 0, 5;

    Matrix2d r = sf::eta(p), e;
    REQUIRE(__D_EQ9(r(0, 0), 0.25));
    REQUIRE(__D_EQ9(r(0, 1), 0.25));
    REQUIRE(__D_EQ9(r(1, 0), 0.5));
    REQUIRE(__D_EQ9(r(1, 1), 0.9721197704469093059));

    r = sf::eta(p, e);
    REQUIRE(__D_EQ9(r(0, 0), 0.25));
    REQUIRE(__D_EQ9(r(0, 1), 0.25));
    REQUIRE(__D_EQ9(r(1, 0), 0.5));
    REQUIRE(__D_EQ9(r(1, 1), 0.9721197704469093059));

    gsl_sf_result re;
    gsl_sf_eta_int_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_eta_int_e(5, &re);
    REQUIRE(re.err == e(1, 1));

    // double
    Matrix2d p2;
    p2 << -5, -1, 0, 5;

    r = sf::eta(p2);
    REQUIRE(__D_EQ9(r(0, 0), 0.25));
    REQUIRE(__D_EQ9(r(0, 1), 0.25));
    REQUIRE(__D_EQ9(r(1, 0), 0.5));
    REQUIRE(__D_EQ9(r(1, 1), 0.9721197704469093059));

    r = sf::eta(p2, e);
    REQUIRE(__D_EQ9(r(0, 0), 0.25));
    REQUIRE(__D_EQ9(r(0, 1), 0.25));
    REQUIRE(__D_EQ9(r(1, 0), 0.5));
    REQUIRE(__D_EQ9(r(1, 1), 0.9721197704469093059));

    gsl_sf_eta_e(-5, &re);
    REQUIRE(re.err == e(0, 0));
    gsl_sf_eta_e(5, &re);
    REQUIRE(re.err == e(1, 1));
}

TEST_CASE("test_hzeta")
{
    Matrix2d p;
    p << 2, 30, 10, 0.9;

    Vector2d r = sf::hzeta(p), e;
    REQUIRE(__D_EQ9(r(0), 0.1051663356816857461));
    REQUIRE(__D_EQ9(r(1), 2.3589824880264765e+01));

    r = sf::hzeta(p, e);
    REQUIRE(__D_EQ9(r(0), 0.1051663356816857461));
    REQUIRE(__D_EQ9(r(1), 2.3589824880264765e+01));

    gsl_sf_result re;
    gsl_sf_hzeta_e(2, 10, &re);
    REQUIRE(re.err == e(0));
    gsl_sf_hzeta_e(30, 0.9, &re);
    REQUIRE(re.err == e(1));
}
