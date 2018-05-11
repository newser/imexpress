#include <catch.hpp>
#include <iostream>
#include <special/exp_integral.h>
#include <test_util.h>

using namespace iexp;
using namespace iexp::sf;

TEST_CASE("test_expint_1")
{
    Vector4d p, r, e;
    r << -1.0, 1.0e-10, 0.1, 1;
    p = expint_1(r);
    REQUIRE(__D_EQ9(p(0), -1.8951178163559367555));
    REQUIRE(__D_EQ9(p(1), 22.448635265138923980));
    REQUIRE(__D_EQ9(p(2), 1.82292395841939066610));
    REQUIRE(__D_EQ9(p(3), 0.21938393439552027368));

    __draw_1d(expint_1, 0.1, 10, 100, 0, 0);

    p = expint_1(r, e);
    REQUIRE(__D_EQ9(p(0), -1.8951178163559367555));
    REQUIRE(__D_EQ9(p(1), 22.448635265138923980));
    REQUIRE(__D_EQ9(p(2), 1.82292395841939066610));
    REQUIRE(__D_EQ9(p(3), 0.21938393439552027368));
    gsl_sf_result re;
    gsl_sf_expint_E1_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_expint_E1_e(r(3), &re);
    REQUIRE(__D_EQ9(e(3), re.err));
}

TEST_CASE("test_expint_2")
{
    Vector4d p, r, e;
    r << -1.0, 0, 1.0 / 4294967296.0, 1.0;
    p = expint_2(r);
    REQUIRE(__D_EQ9(p(0), 0.8231640121031084799));
    REQUIRE(__D_EQ9(p(1), 1.0));
    REQUIRE(__D_EQ9(p(2), 0.9999999947372139168));
    REQUIRE(__D_EQ9(p(3), 0.14849550677592204792));

    __draw_1d(expint_2, 0.1, 10, 100, 0, 0);

    p = expint_2(r, e);
    REQUIRE(__D_EQ9(p(0), 0.8231640121031084799));
    REQUIRE(__D_EQ9(p(1), 1.0));
    REQUIRE(__D_EQ9(p(2), 0.9999999947372139168));
    REQUIRE(__D_EQ9(p(3), 0.14849550677592204792));
    gsl_sf_result re;
    gsl_sf_expint_E2_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_expint_E2_e(r(3), &re);
    REQUIRE(__D_EQ9(e(3), re.err));
}

TEST_CASE("test_expint")
{
    Vector4d p, r, e;
    r << 0, 1.0 / 4294967296.0, 0.1, 10.0;
    p = expint(3, r);
    REQUIRE(__D_EQ9(p(0), 0.5));
    REQUIRE(__D_EQ9(p(1), 0.499999999767169356972));
    REQUIRE(__D_EQ9(p(2), 0.4162914579082787612543));
    REQUIRE(__D_EQ9(p(3), 0.000003548762553084381959981));

    p = expint(3, r, e);
    REQUIRE(__D_EQ9(p(0), 0.5));
    REQUIRE(__D_EQ9(p(1), 0.499999999767169356972));
    REQUIRE(__D_EQ9(p(2), 0.4162914579082787612543));
    REQUIRE(__D_EQ9(p(3), 0.000003548762553084381959981));
    gsl_sf_result re;
    gsl_sf_expint_En_e(3, r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_expint_En_e(3, r(3), &re);
    REQUIRE(__D_EQ9(e(3), re.err));
}

TEST_CASE("test_expint_ei")
{
    Vector3d p, r, e;
    r << -1.0, 1.0 / 4294967296.0, 1.0;
    p = expint_ei(r);
    REQUIRE(__D_EQ9(p(0), -0.21938393439552027368));
    REQUIRE(__D_EQ9(p(1), -21.603494112783886397));
    REQUIRE(__D_EQ9(p(2), 1.8951178163559367555));

    __draw_1d(expint_ei, 0.1, 10, 100, 0, 0);

    p = expint_ei(r, e);
    REQUIRE(__D_EQ9(p(0), -0.21938393439552027368));
    REQUIRE(__D_EQ9(p(1), -21.603494112783886397));
    REQUIRE(__D_EQ9(p(2), 1.8951178163559367555));
    gsl_sf_result re;
    gsl_sf_expint_Ei_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_expint_Ei_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_expint_ei3")
{
    Vector4d p, r, e;
    r << 1.0e-10, 0.1, 1.0, 10.0;
    p = expint_ei3(r);
    REQUIRE(__D_EQ9(p(0), 1.0e-10));
    REQUIRE(__D_EQ9(p(1), 0.09997500714119079665122));
    REQUIRE(__D_EQ9(p(2), 0.80751118213967145285833));
    REQUIRE(__D_EQ9(p(3), 0.89297951156924921121856));

    __draw_1d(expint_ei3, 0.1, 10, 100, 0, 0);

    p = expint_ei3(r, e);
    REQUIRE(__D_EQ9(p(0), 1.0e-10));
    REQUIRE(__D_EQ9(p(1), 0.09997500714119079665122));
    REQUIRE(__D_EQ9(p(2), 0.80751118213967145285833));
    REQUIRE(__D_EQ9(p(3), 0.89297951156924921121856));
    gsl_sf_result re;
    gsl_sf_expint_3_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_expint_3_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_coshint")
{
    Vector4d p, r, e;
    r << -1.0, 1.0 / 4294967296.0, 1.0, 10.0;
    p = coshint(r);
    REQUIRE(__D_EQ9(p(0), 0.8378669409802082409));
    REQUIRE(__D_EQ9(p(1), -21.603494113016717041));
    REQUIRE(__D_EQ9(p(2), 0.8378669409802082409));
    REQUIRE(__D_EQ9(p(3), 1246.1144860424544147));

    __draw_1d(coshint, 0.1, 10, 100, 0, 0);

    p = coshint(r, e);
    REQUIRE(__D_EQ9(p(0), 0.8378669409802082409));
    REQUIRE(__D_EQ9(p(1), -21.603494113016717041));
    REQUIRE(__D_EQ9(p(2), 0.8378669409802082409));
    REQUIRE(__D_EQ9(p(3), 1246.1144860424544147));
    gsl_sf_result re;
    gsl_sf_Chi_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_Chi_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_sinhint")
{
    Vector4d p, r, e;
    r << -1.0, 1.0 / 4294967296.0, 1.0, 10.0;
    p = sinhint(r);
    REQUIRE(__D_EQ9(p(0), -1.0572508753757285146));
    REQUIRE(__D_EQ9(p(1), 2.3283064365386962891e-10));
    REQUIRE(__D_EQ9(p(2), 1.0572508753757285146));
    REQUIRE(__D_EQ9(p(3), 1246.1144901994233444));

    __draw_1d(sinhint, 0.1, 10, 100, 0, 0);

    p = sinhint(r, e);
    REQUIRE(__D_EQ9(p(0), -1.0572508753757285146));
    REQUIRE(__D_EQ9(p(1), 2.3283064365386962891e-10));
    REQUIRE(__D_EQ9(p(2), 1.0572508753757285146));
    REQUIRE(__D_EQ9(p(3), 1246.1144901994233444));
    gsl_sf_result re;
    gsl_sf_Shi_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_Shi_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_sinint")
{
    Vector4d p, r, e;
    r << -1.0, 1.0e-05, 1.0, 10.0;
    p = sinint(r);
    REQUIRE(__D_EQ9(p(0), -0.9460830703671830149));
    REQUIRE(__D_EQ9(p(1), 9.999999999944444444e-06));
    REQUIRE(__D_EQ9(p(2), 0.9460830703671830149));
    REQUIRE(__D_EQ9(p(3), 1.6583475942188740493));

    __draw_1d(sinint, 0.1, 10, 100, 0, 0);

    p = sinint(r, e);
    REQUIRE(__D_EQ9(p(0), -0.9460830703671830149));
    REQUIRE(__D_EQ9(p(1), 9.999999999944444444e-06));
    REQUIRE(__D_EQ9(p(2), 0.9460830703671830149));
    REQUIRE(__D_EQ9(p(3), 1.6583475942188740493));
    gsl_sf_result re;
    gsl_sf_Si_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_Si_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_cosint")
{
    Vector4d p, r, e;
    r << 1.0 / 4294967296.0, 1.0 / 65536.0, 1.0, 10.0;
    p = cosint(r);
    REQUIRE(__D_EQ9(p(0), -21.603494113016717041));
    REQUIRE(__D_EQ9(p(1), -10.513139224115799751));
    REQUIRE(__D_EQ9(p(2), 0.3374039229009681347));
    REQUIRE(__D_EQ9(p(3), -0.04545643300445537263));

    __draw_1d(cosint, 0.1, 10, 100, 0, 0);

    p = cosint(r, e);
    REQUIRE(__D_EQ9(p(0), -21.603494113016717041));
    REQUIRE(__D_EQ9(p(1), -10.513139224115799751));
    REQUIRE(__D_EQ9(p(2), 0.3374039229009681347));
    REQUIRE(__D_EQ9(p(3), -0.04545643300445537263));
    gsl_sf_result re;
    gsl_sf_Ci_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_Ci_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}

TEST_CASE("test_atanint")
{
    Vector4d p, r, e;
    r << 1.0e-10, 1.0e-05, 1.0, 10.0;
    p = cosint(r);
    REQUIRE(__D_EQ9(p(0), 1.0e-10));
    REQUIRE(__D_EQ9(p(1), 9.99999999988888888889e-06));
    REQUIRE(__D_EQ9(p(2), 0.91596559417721901505));
    REQUIRE(__D_EQ9(p(3), 3.71678149306806859029));

    __draw_1d(cosint, 0.1, 10, 100, 0, 0);

    p = cosint(r, e);
    REQUIRE(__D_EQ9(p(0), 1.0e-10));
    REQUIRE(__D_EQ9(p(1), 9.99999999988888888889e-06));
    REQUIRE(__D_EQ9(p(2), 0.91596559417721901505));
    REQUIRE(__D_EQ9(p(3), 3.71678149306806859029));
    gsl_sf_result re;
    gsl_sf_atanint_e(r(0), &re);
    REQUIRE(__D_EQ9(e(0), re.err));
    gsl_sf_atanint_e(r(2), &re);
    REQUIRE(__D_EQ9(e(2), re.err));
}
