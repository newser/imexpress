#include <catch.hpp>
#include <iostream>
#include <special/elliptic_integral.h>
#include <special/elliptic_jacobi.h>
#include <test_util.h>

using namespace iexp;
using namespace iexp::sf;

TEST_CASE("test_ellint_k")
{
    Matrix<double, 3, 1> p, r;
    p << 0.99, 0.50, 0.010;
    r = sf::ellint_k(p);
    REQUIRE(__D_EQ9(r(0), 3.3566005233611923760));
    REQUIRE(__D_EQ9(r(1), 1.6857503548125960429));
    REQUIRE(__D_EQ9(r(2), 1.5708355989121522360));

    __draw_1d(ellint_k, 0.01, 0.99, 100, 0, 0);

    Matrix<double, 3, 1> e;
    r = sf::ellint_k(p, e);
    REQUIRE(__D_EQ9(r(0), 3.3566005233611923760));
    REQUIRE(__D_EQ9(r(1), 1.6857503548125960429));
    REQUIRE(__D_EQ9(r(2), 1.5708355989121522360));

    gsl_sf_result re;
    gsl_sf_ellint_Kcomp_e(p(0), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_Kcomp_e(p(2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_ki")
{
    Matrix<double, 3, 2, RowMajor> p;
    p.col(0).setConstant(M_PI / 3.0);
    p.col(1) << 0.99, 0.50, 0.010;

    Matrix<double, 3, 1> r;
    r = sf::ellint_ki(p);
    REQUIRE(__D_EQ9(r(0), 1.3065333392738766762));
    REQUIRE(__D_EQ9(r(1), 1.0895506700518854093));
    REQUIRE(__D_EQ9(r(2), 1.0472129063770918952));

    Matrix<double, 3, 1> e;
    r = sf::ellint_ki(p, e);
    REQUIRE(__D_EQ9(r(0), 1.3065333392738766762));
    REQUIRE(__D_EQ9(r(1), 1.0895506700518854093));
    REQUIRE(__D_EQ9(r(2), 1.0472129063770918952));

    gsl_sf_result re;
    gsl_sf_ellint_F_e(p(0, 0), p(0, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_F_e(p(2, 0), p(2, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_e")
{
    Matrix<double, 3, 1> p, r;
    p << 0.99, 0.50, 0.010;
    r = sf::ellint_e(p);
    REQUIRE(__D_EQ9(r(0), 1.0284758090288040010));
    REQUIRE(__D_EQ9(r(1), 1.4674622093394271555));
    REQUIRE(__D_EQ9(r(2), 1.5707570561503852873));

    __draw_1d(ellint_e, 0.01, 0.99, 100, 0, 0);

    Matrix<double, 3, 1> e;
    r = sf::ellint_e(p, e);
    REQUIRE(__D_EQ9(r(0), 1.0284758090288040010));
    REQUIRE(__D_EQ9(r(1), 1.4674622093394271555));
    REQUIRE(__D_EQ9(r(2), 1.5707570561503852873));

    gsl_sf_result re;
    gsl_sf_ellint_Ecomp_e(p(0), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_Ecomp_e(p(2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_ei")
{
    Matrix<double, 3, 2, RowMajor> p;
    p.col(0).setConstant(M_PI / 3.0);
    p.col(1) << 0.99, 0.50, 0.010;

    Matrix<double, 3, 1> r;
    r = sf::ellint_ei(p);
    REQUIRE(__D_EQ9(r(0), 0.8704819220377943536));
    REQUIRE(__D_EQ9(r(1), 1.0075555551444720293));
    REQUIRE(__D_EQ9(r(2), 1.0471821963889481104));

    Matrix<double, 3, 1> e;
    r = sf::ellint_ei(p, e);
    REQUIRE(__D_EQ9(r(0), 0.8704819220377943536));
    REQUIRE(__D_EQ9(r(1), 1.0075555551444720293));
    REQUIRE(__D_EQ9(r(2), 1.0471821963889481104));

    gsl_sf_result re;
    gsl_sf_ellint_E_e(p(0, 0), p(0, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_E_e(p(2, 0), p(2, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_p")
{
    Matrix<double, 3, 2, RowMajor> p;
    p.col(0) << 0.99, 0.50, 0.010;
    p.col(1).setConstant(0.1);

    Matrix<double, 3, 1> r;
    r = sf::ellint_p(p);
    REQUIRE(__D_EQ9(r(0), 3.13792612351836506315593));
    REQUIRE(__D_EQ9(r(1), 1.60455249360848890075108));
    REQUIRE(__D_EQ9(r(2), 1.49773208536003801277453));

    Matrix<double, 3, 1> e;
    r = sf::ellint_p(p, e);
    REQUIRE(__D_EQ9(r(0), 3.13792612351836506315593));
    REQUIRE(__D_EQ9(r(1), 1.60455249360848890075108));
    REQUIRE(__D_EQ9(r(2), 1.49773208536003801277453));

    gsl_sf_result re;
    gsl_sf_ellint_Pcomp_e(p(0, 0), p(0, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_Pcomp_e(p(2, 0), p(2, 1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_pi")
{
    Matrix<double, 3, 3, RowMajor> p;
    p.col(0).setConstant(M_PI / 3.0);
    p.col(1) << 0.99, 0.50, 0.010;
    p.col(2).setConstant(0.5);

    Matrix<double, 3, 1> r;
    r = sf::ellint_pi(p);
    REQUIRE(__D_EQ9(r(0), 1.1288726598764099882));
    REQUIRE(__D_EQ9(r(1), 0.9570574331323584890));
    REQUIRE(__D_EQ9(r(2), 0.9228868127118118465));

    Matrix<double, 3, 1> e;
    r = sf::ellint_pi(p, e);
    REQUIRE(__D_EQ9(r(0), 1.1288726598764099882));
    REQUIRE(__D_EQ9(r(1), 0.9570574331323584890));
    REQUIRE(__D_EQ9(r(2), 0.9228868127118118465));

    gsl_sf_result re;
    gsl_sf_ellint_P_e(p(0, 0), p(0, 1), p(0, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_P_e(p(2, 0), p(2, 1), p(2, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_ellint_rc")
{
    Matrix<double, 1, 2> p;
    p << 1.0, 2.0;

    Matrix<double, 1, 1> r;
    r = sf::ellint_rc(p);
    REQUIRE(__D_EQ9(r(0), 0.7853981633974482));

    Matrix<double, 1, 1> e;
    r = sf::ellint_rc(p, e);
    REQUIRE(__D_EQ9(r(0), 0.7853981633974482));

    gsl_sf_result re;
    gsl_sf_ellint_RC_e(p(0), p(1), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
}

TEST_CASE("test_ellint_rd")
{
    Matrix<double, 2, 3, RowMajor> p;
    p.row(0) << 5.0e-11, 1.0e-10, 1.0;
    p.row(1) << 1.0, 2.0, 3.0;

    Matrix<double, 2, 1> r;
    r = sf::ellint_rd(p);
    REQUIRE(__D_EQ9(r(0), 34.0932594919337362));
    REQUIRE(__D_EQ9(r(1), 0.2904602810289906));

    Matrix<double, 2, 1> e;
    r = sf::ellint_rd(p, e);
    REQUIRE(__D_EQ9(r(0), 34.0932594919337362));
    REQUIRE(__D_EQ9(r(1), 0.2904602810289906));

    gsl_sf_result re;
    gsl_sf_ellint_RD_e(p(0, 0), p(0, 1), p(0, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_RD_e(p(1, 0), p(1, 1), p(1, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(1), re.err));
}

TEST_CASE("test_ellint_rf")
{
    Matrix<double, 2, 3, RowMajor> p;
    p.row(0) << 5.0e-11, 1.0e-10, 1.0;
    p.row(1) << 1.0, 2.0, 3.0;

    Matrix<double, 2, 1> r;
    r = sf::ellint_rf(p);
    REQUIRE(__D_EQ9(r(0), 12.36441982979439));
    REQUIRE(__D_EQ9(r(1), 0.7269459354689082));

    Matrix<double, 2, 1> e;
    r = sf::ellint_rf(p, e);
    REQUIRE(__D_EQ9(r(0), 12.36441982979439));
    REQUIRE(__D_EQ9(r(1), 0.7269459354689082));

    gsl_sf_result re;
    gsl_sf_ellint_RF_e(p(0, 0), p(0, 1), p(0, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_ellint_RF_e(p(1, 0), p(1, 1), p(1, 2), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(1), re.err));
}

TEST_CASE("test_ellint_rj")
{
    Matrix<double, 1, 4> p;
    p << 2.0, 3.0, 4.0, 5.0;

    Matrix<double, 1, 1> r;
    r = sf::ellint_rj(p);
    REQUIRE(__D_EQ9(r(0), 0.1429757966715675));

    Matrix<double, 1, 1> e;
    r = sf::ellint_rj(p, e);
    REQUIRE(__D_EQ9(r(0), 0.1429757966715675));

    gsl_sf_result re;
    gsl_sf_ellint_RJ_e(p(0), p(1), p(2), p(3), GSL_PREC_DOUBLE, &re);
    REQUIRE(__D_EQ9(e(0), re.err));
}

TEST_CASE("test_elljac")
{
    Matrix<double, 2, 2> p;
    p.col(0) << 0.5, 0.5;
    p.col(1) << 1.0, 0.3;

    Matrix<std::tuple<double, double, double>, 1, 2> r;
    r = elljac(p);
    REQUIRE(__D_EQ9(std::get<0>(r(0)), 0.4707504736556572833));
    REQUIRE(__D_EQ9(std::get<1>(r(0)), 0.8822663948904402865));
    REQUIRE(__D_EQ9(std::get<2>(r(0)), 0.9429724257773856873));

    REQUIRE(__D_EQ9(std::get<0>(r(1)), 0.8187707145344889190));
    REQUIRE(__D_EQ9(std::get<1>(r(1)), 0.5741206467465548795));
    REQUIRE(__D_EQ9(std::get<2>(r(1)), 0.8938033089590823040));
}
