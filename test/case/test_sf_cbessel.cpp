#include <catch.hpp>
#include <iostream>
#include <special/cylindrical_bessel.h>
#include <special/modified_cylindrical_bessel.h>
#include <test_util.h>

TEST_CASE("sf_cbessel_j")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // J0
    m << 0.1, 2.0, 100.0, 1.0e+10;
    m2 = iexp::sf::cbessel_j(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.99750156206604003230));
    REQUIRE(__D_EQ9(m2(1, 0), 0.22389077914123566805));
    REQUIRE(__D_EQ9(m2(2, 0), 0.019985850304223122424));
    REQUIRE(__D_EQ9(m2(3, 0), 2.1755917502468917269e-06));

    m2 = iexp::sf::cbessel_j(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.99750156206604003230));
    REQUIRE(__D_EQ9(m2(3, 0), 2.1755917502468917269e-06));
    gsl_sf_bessel_j0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_j0_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    // J1
    m2 = iexp::sf::cbessel_j(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.04993752603624199756));
    REQUIRE(__D_EQ9(m2(1, 0), 0.57672480775687338720));
    REQUIRE(__D_EQ9(m2(2, 0), -0.07714535201411215803));
    REQUIRE(__D_EQ9(m2(3, 0), -7.676508175684157103e-06));

    m2 = iexp::sf::cbessel_j(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.04993752603624199756));
    REQUIRE(__D_EQ9(m2(3, 0), -7.676508175684157103e-06));
    gsl_sf_bessel_j1_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_j1_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    // Jn
    m.resize(2);
    m << 900.0, 15000.0;
    m2 = iexp::sf::cbessel_j(2, m);
    REQUIRE(__D_EQ9(m2(0, 0), -0.019974345269680646400));
    REQUIRE(__D_EQ9(m2(1, 0), -0.0020455820181216382666));

    m2 = iexp::sf::cbessel_j(2, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -0.019974345269680646400));
    REQUIRE(__D_EQ9(m2(1, 0), -0.0020455820181216382666));
    gsl_sf_bessel_Jn_e(2, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_Jn_e(2, m(1, 0), &r);
    REQUIRE(__D_EQ9(e(1, 0), r.err));
}

TEST_CASE("sf_cbessel_y")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // J0
    m << 0.1, 2.0, 256.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), -1.5342386513503668441));
    REQUIRE(__D_EQ9(m2(1, 0), 0.5103756726497451196));
    REQUIRE(__D_EQ9(m2(2, 0), -0.03381290171792454909));
    REQUIRE(__D_EQ9(m2(3, 0), 3.657903190017678681e-06));

    m2 = iexp::sf::cbessel_y(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -1.5342386513503668441));
    REQUIRE(__D_EQ9(m2(3, 0), 3.657903190017678681e-06));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    // J1
    m << 0.1, 2.0, 100.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), -6.45895109470202698800));
    REQUIRE(__D_EQ9(m2(1, 0), -0.10703243154093754689));
    REQUIRE(__D_EQ9(m2(2, 0), -0.020372312002759793305));
    REQUIRE(__D_EQ9(m2(3, 0), 0.000011612249378370766284));

    m2 = iexp::sf::cbessel_y(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -6.45895109470202698800));
    REQUIRE(__D_EQ9(m2(3, 0), 0.000011612249378370766284));
    gsl_sf_bessel_y0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_y0_e(m(3, 0), &r);
    REQUIRE(__D_EQ9(e(3, 0), r.err));

    // Jn
    m.resize(2);
    m << 100.0, 4294967296.0;
    m2 = iexp::sf::cbessel_y(100, m);
    REQUIRE(__D_EQ9(m2(0, 0), -0.16692141141757650654));
    REQUIRE(__D_EQ9(m2(1, 0), 3.657889671577715808e-06));

    m2 = iexp::sf::cbessel_y(100, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), -0.16692141141757650654));
    REQUIRE(__D_EQ9(m2(1, 0), 3.657889671577715808e-06));
    gsl_sf_bessel_Yn_e(100, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_Yn_e(100, m(1, 0), &r);
    REQUIRE(__D_EQ9(e(1, 0), r.err));
}

TEST_CASE("sf_mcbessel_i")
{
    iexp::ArrayXd m(3), m2(3), e(3);
    gsl_sf_result r;

    // I0
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::mcbessel_i(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 1.0025015629340956014));
    REQUIRE(__D_EQ9(m2(1, 0), 2.2795853023360672674));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0737517071310738235e+42));

    m2 = iexp::sf::mcbessel_i(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 1.0025015629340956014));
    REQUIRE(__D_EQ9(m2(1, 0), 2.2795853023360672674));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0737517071310738235e+42));
    gsl_sf_bessel_I0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_I0_e(m(2, 0), &r);
    REQUIRE(__D_EQ7(e(2, 0), r.err)); // ??

    // I1
    m2 = iexp::sf::mcbessel_i(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.05006252604709269211));
    REQUIRE(__D_EQ9(m2(1, 0), 1.59063685463732906340));
    REQUIRE(__D_EQ_Nep(m2(2, 0), 1.0683693903381627E+42, 1000));

    m2 = iexp::sf::mcbessel_i(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.05006252604709269211));
    REQUIRE(__D_EQ9(m2(1, 0), 1.59063685463732906340));
    REQUIRE(__D_EQ_Nep(m2(2, 0), 1.0683693903381627E+42, 1000));
    gsl_sf_bessel_I1_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_I1_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // In
    m2 = iexp::sf::mcbessel_i(4, m);
    REQUIRE(__D_EQ9(m2(0, 0), 2.6054690212996573677e-07));

    m2 = iexp::sf::mcbessel_i(4, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 2.6054690212996573677e-07));
    gsl_sf_bessel_In_e(4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
}

TEST_CASE("sf_mcbessel_k")
{
    iexp::ArrayXd m(3), m2(3), e(3);
    gsl_sf_result r;

    // K0
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::mcbessel_k(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 2.4270690247020166125));
    REQUIRE(__D_EQ9(m2(1, 0), 0.11389387274953343565));
    REQUIRE(__D_EQ9(m2(2, 0), 4.656628229175902019e-45));

    m2 = iexp::sf::mcbessel_k(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 2.4270690247020166125));
    REQUIRE(__D_EQ9(m2(1, 0), 0.11389387274953343565));
    REQUIRE(__D_EQ9(m2(2, 0), 4.656628229175902019e-45));
    gsl_sf_bessel_K0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_K0_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // K1
    m << 0.1, 2.0, 100.0;
    m2 = iexp::sf::mcbessel_k(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 9.853844780870606135));
    REQUIRE(__D_EQ9(m2(1, 0), 0.13986588181652242728));
    REQUIRE(__D_EQ9(m2(2, 0), 4.679853735636909287e-45));

    m2 = iexp::sf::mcbessel_k(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 9.853844780870606135));
    REQUIRE(__D_EQ9(m2(1, 0), 0.13986588181652242728));
    REQUIRE(__D_EQ9(m2(2, 0), 4.679853735636909287e-45));
    gsl_sf_bessel_K0_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_K0_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // Kn
    m.resize(2);
    m << 0.1, 2.0;
    m2 = iexp::sf::mcbessel_k(4, m);
    REQUIRE(__D_EQ9(m2(0, 0), 479600.2497925682849));

    m2 = iexp::sf::mcbessel_k(4, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 479600.2497925682849));
    gsl_sf_bessel_Kn_e(4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
}

TEST_CASE("sf_mcbessel_i_scaled")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // I0
    m << 1e-10, 0.1, 2.0, 65536.0;
    m2 = iexp::sf::mcbessel_i<true>(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.99999999990000000001));
    REQUIRE(__D_EQ9(m2(1, 0), 0.90710092578230109640));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30850832255367103953));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0015583712551952223537));

    m2 = iexp::sf::mcbessel_i<true>(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.99999999990000000001));
    REQUIRE(__D_EQ9(m2(1, 0), 0.90710092578230109640));
    REQUIRE(__D_EQ9(m2(2, 0), 0.30850832255367103953));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0015583712551952223537));
    gsl_sf_bessel_I0_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_I0_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err)); // ??

    // I1
    m << 0.1, 2.0, 100.0, 65536.0;
    m2 = iexp::sf::mcbessel_i<true>(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 0.04529844680880932501));
    REQUIRE(__D_EQ9(m2(1, 0), 0.21526928924893765916));
    REQUIRE(__D_EQ9(m2(2, 0), 0.039744153025130252));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0015583593657207350452));

    m2 = iexp::sf::mcbessel_i<true>(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 0.04529844680880932501));
    REQUIRE(__D_EQ9(m2(1, 0), 0.21526928924893765916));
    REQUIRE(__D_EQ9(m2(2, 0), 0.039744153025130252));
    REQUIRE(__D_EQ9(m2(3, 0), 0.0015583593657207350452));
    gsl_sf_bessel_I1_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_I1_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // In
    m2 = iexp::sf::mcbessel_i<true>(-4, m);
    REQUIRE(__D_EQ9(m2(0, 0), 2.3575258620054605307e-07));

    m2 = iexp::sf::mcbessel_i<true>(-4, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 2.3575258620054605307e-07));
    gsl_sf_bessel_In_scaled_e(-4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
}

TEST_CASE("sf_mcbessel_k_scaled")
{
    iexp::ArrayXd m(4), m2(4), e(4);
    gsl_sf_result r;

    // K0
    m << 0.1, 1.95, 2.0, 100.0;
    m2 = iexp::sf::mcbessel_k<true>(0, m);
    REQUIRE(__D_EQ9(m2(0, 0), 2.6823261022628943831));
    REQUIRE(__D_EQ9(m2(1, 0), 0.8513330938802157074894));
    REQUIRE(__D_EQ9(m2(2, 0), 0.8415682150707714179));
    REQUIRE(__D_EQ9(m2(3, 0), 0.1251756216591265789));

    m2 = iexp::sf::mcbessel_k<true>(0, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 2.6823261022628943831));
    REQUIRE(__D_EQ9(m2(1, 0), 0.8513330938802157074894));
    REQUIRE(__D_EQ9(m2(2, 0), 0.8415682150707714179));
    REQUIRE(__D_EQ9(m2(3, 0), 0.1251756216591265789));
    gsl_sf_bessel_K0_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_K0_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // K1
    m2 = iexp::sf::mcbessel_k<true>(1, m);
    REQUIRE(__D_EQ9(m2(0, 0), 10.890182683049696574));
    REQUIRE(__D_EQ9(m2(1, 0), 1.050086915104152747182));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0334768470686885732));
    REQUIRE(__D_EQ9(m2(3, 0), 0.1257999504795785293));

    m2 = iexp::sf::mcbessel_k<true>(1, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 10.890182683049696574));
    REQUIRE(__D_EQ9(m2(1, 0), 1.050086915104152747182));
    REQUIRE(__D_EQ9(m2(2, 0), 1.0334768470686885732));
    REQUIRE(__D_EQ9(m2(3, 0), 0.1257999504795785293));
    gsl_sf_bessel_K0_scaled_e(m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
    gsl_sf_bessel_K0_scaled_e(m(2, 0), &r);
    REQUIRE(__D_EQ9(e(2, 0), r.err));

    // Kn
    m2 = iexp::sf::mcbessel_k<true>(4, m);
    REQUIRE(__D_EQ9(m2(0, 0), 530040.2483725626207));

    m2 = iexp::sf::mcbessel_k<true>(4, m, e);
    REQUIRE(__D_EQ9(m2(0, 0), 530040.2483725626207));
    gsl_sf_bessel_Kn_scaled_e(4, m(0, 0), &r);
    REQUIRE(__D_EQ9(e(0, 0), r.err));
}
