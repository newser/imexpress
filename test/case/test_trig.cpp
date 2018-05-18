#include <catch.hpp>
#include <iostream>
#include <special/trigonometry.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_sin")
{
    Matrix2d p, r, e;
    p << -10.0, 1, 1000, 1048576.75;
    r = sf::sin(p);
    REQUIRE(__D_EQ9(r(0, 0), 0.5440211108893698134));
    REQUIRE(__D_EQ9(r(1, 1), 0.8851545351115651914));

    r = sf::sin(p, e);
    REQUIRE(__D_EQ9(r(0, 0), 0.5440211108893698134));
    REQUIRE(__D_EQ9(r(1, 1), 0.8851545351115651914));

    gsl_sf_result re;
    gsl_sf_sin_e(-10.0, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_sin_e(1048576.75, &re);
    REQUIRE(e(1, 1) == re.err);
}

TEST_CASE("test_cos")
{
    Matrix2d p, r, e;
    p << -10.0, 1, 1000, 1048576.75;
    r = sf::cos(p);
    REQUIRE(__D_EQ9(r(0, 0), -0.8390715290764524523));
    REQUIRE(__D_EQ9(r(1, 1), 0.4652971620066351799));

    r = sf::cos(p, e);
    REQUIRE(__D_EQ9(r(0, 0), -0.8390715290764524523));
    REQUIRE(__D_EQ9(r(1, 1), 0.4652971620066351799));

    gsl_sf_result re;
    gsl_sf_cos_e(-10.0, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_cos_e(1048576.75, &re);
    REQUIRE(e(1, 1) == re.err);
}

TEST_CASE("test_sinc")
{
    Matrix2d p, r, e;
    p << 1.0 / 1024.0, 1.0 / 2.0, 80.5, 100.5;
    r = sf::sinc(p);
    REQUIRE(__D_EQ9(r(0, 0), 0.9999984312693665404));
    REQUIRE(__D_EQ9(r(1, 1), 0.0031672625490924445));

    r = sf::sinc(p, e);
    REQUIRE(__D_EQ9(r(0, 0), 0.9999984312693665404));
    REQUIRE(__D_EQ9(r(1, 1), 0.0031672625490924445));

    gsl_sf_result re;
    gsl_sf_sinc_e(1.0 / 1024.0, &re);
    REQUIRE(e(0, 0) == re.err);
    gsl_sf_sinc_e(100.5, &re);
    REQUIRE(e(1, 1) == re.err);
}

TEST_CASE("test_hypot")
{
    Matrix2d p;
    Vector2d r, e;
    p << 1, 2, 3, 4;
    r = sf::hypot(p);
    REQUIRE(__D_EQ9(r(0), std::hypot(1, 3)));
    REQUIRE(__D_EQ9(r(1), std::hypot(2, 4)));

    r = sf::hypot(p, e);
    REQUIRE(__D_EQ9(r(0), std::hypot(1, 3)));
    REQUIRE(__D_EQ9(r(1), std::hypot(2, 4)));

    gsl_sf_result re;
    gsl_sf_hypot_e(1, 3, &re);
    REQUIRE(e(0) == re.err);
    gsl_sf_hypot_e(2, 4, &re);
    REQUIRE(e(1) == re.err);
}
