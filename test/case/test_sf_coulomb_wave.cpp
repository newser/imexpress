#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <special/coulomb_wave.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_coulomb_f")
{
    iexp::ArrayXd l(2), eta(2), x(2), e(2), f(2), g(2);
    iexp::ArrayXi k;
    gsl_sf_result gf, gg, gfp, ggp;
    double ef, eg;

    // F
    l = 0.0;
    eta << -1000.0, -50.0;
    x << 1.0, 5.0;
    f = sf::coulw_f(l, eta, x);
    REQUIRE(__D_EQ9(f(0, 0), 9.68222518991341e-02));
    REQUIRE(__D_EQ9(f(1, 0), 1.52236975714236e-01));

    f = sf::coulw_f(l, eta, x, e);
    REQUIRE(__D_EQ9(f(0, 0), 9.68222518991341e-02));
    REQUIRE(__D_EQ9(f(1, 0), 1.52236975714236e-01));
    gsl_sf_coulomb_wave_FG_e(-1000.0,
                             1.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(0, 0), gf.err));
    gsl_sf_coulomb_wave_FG_e(-50.0,
                             5.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(1, 0), gf.err));

    // G
    g = sf::coulw_g(l, eta, x);
    REQUIRE(__D_EQ9(g(0, 0), 1.13936784379472e-01));
    REQUIRE(__D_EQ9(g(1, 0), 4.41680690236251e-01));

    g = sf::coulw_g(l, eta, x, e);
    REQUIRE(__D_EQ9(g(0, 0), 1.13936784379472e-01));
    REQUIRE(__D_EQ9(g(1, 0), 4.41680690236251e-01));
    gsl_sf_coulomb_wave_FG_e(-1000.0,
                             1.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(0, 0), gg.err));
    gsl_sf_coulomb_wave_FG_e(-50.0,
                             5.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(1, 0), gg.err));

    // Fp
    g = sf::coulw_df(l, eta, x);
    REQUIRE(__D_EQ9(g(0, 0), 5.12063396274631e+00));
    REQUIRE(__D_EQ9(g(1, 0), 2.03091041166137e+00));

    g = sf::coulw_df(l, eta, x, e);
    REQUIRE(__D_EQ9(g(0, 0), 5.12063396274631e+00));
    REQUIRE(__D_EQ9(g(1, 0), 2.03091041166137e+00));
    gsl_sf_coulomb_wave_FG_e(-1000.0,
                             1.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(0, 0), gfp.err));
    gsl_sf_coulomb_wave_FG_e(-50.0,
                             5.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(1, 0), gfp.err));

    // Gp
    g = sf::coulw_dg(l, eta, x);
    REQUIRE(__D_EQ9(g(0, 0), -4.30243486522438e+00));
    REQUIRE(__D_EQ9(g(1, 0), -6.76485374766869e-01));

    g = sf::coulw_dg(l, eta, x, e);
    REQUIRE(__D_EQ9(g(0, 0), -4.30243486522438e+00));
    REQUIRE(__D_EQ9(g(1, 0), -6.76485374766869e-01));
    gsl_sf_coulomb_wave_FG_e(-1000.0,
                             1.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(0, 0), gfp.err));
    gsl_sf_coulomb_wave_FG_e(-50.0,
                             5.0,
                             0.0,
                             0,
                             &gf,
                             &gfp,
                             &gg,
                             &ggp,
                             &ef,
                             &eg);
    REQUIRE(__D_EQ9(e(1, 0), gfp.err));
}

TEST_CASE("test_coulomb_fg")
{
    iexp::ArrayXd l(2), eta(2), x(2), e(2), f(2), g(2), df(2), dg(2);
    iexp::ArrayXi k(2);
    gsl_sf_result gf, gg, gfp, ggp;
    double ef, eg;

    // fg
    l = 0.0;
    k = 0;
    eta << -1000.0, -50.0;
    x << 1.0, 5.0;
    sf::coulw_fg(l, k, eta, x, f, g);
    REQUIRE(__D_EQ9(f(0, 0), 9.68222518991341e-02));
    REQUIRE(__D_EQ9(f(1, 0), 1.52236975714236e-01));
    REQUIRE(__D_EQ9(g(0, 0), 1.13936784379472e-01));
    REQUIRE(__D_EQ9(g(1, 0), 4.41680690236251e-01));

    // dfg
    l = 0.0;
    k = 0;
    eta << -1000.0, -50.0;
    x << 1.0, 5.0;
    sf::coulw_dfg(l, k, eta, x, f, g, df, dg);
    REQUIRE(__D_EQ9(f(0, 0), 9.68222518991341e-02));
    REQUIRE(__D_EQ9(f(1, 0), 1.52236975714236e-01));
    REQUIRE(__D_EQ9(g(0, 0), 1.13936784379472e-01));
    REQUIRE(__D_EQ9(g(1, 0), 4.41680690236251e-01));
    REQUIRE(__D_EQ9(df(0, 0), 5.12063396274631e+00));
    REQUIRE(__D_EQ9(df(1, 0), 2.03091041166137e+00));
    REQUIRE(__D_EQ9(dg(0, 0), -4.30243486522438e+00));
    REQUIRE(__D_EQ9(dg(1, 0), -6.76485374766869e-01));
}

TEST_CASE("test_coulomb_cl")
{
}
