#include <catch.hpp>
#include <integral/qng.h>
#include <iostream>
#include <test_util.h>

using namespace iexp;
using namespace std;

static double foo(double x)
{
    return x * x;
}

TEST_CASE("integral_func")
{
    integral::unary_func<double> f([](double x) -> double { return 2 * x; });
    int dummy[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    double v = f.gsl()->function(1.2, f.gsl()->params);
    REQUIRE(v == 2.4);

    // integral::unary_func<double> f2([](float x) -> float { return x; });

    class t_op
    {
      public:
        double operator()(double x)
        {
            return x + 2;
        }
    };
    integral::unary_func<double> f2((t_op()));
    int dummy2[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    v = f2.gsl()->function(1.2, f2.gsl()->params);
    REQUIRE(v == 3.2);

    integral::unary_func<double> f3(foo);
    int dummy3[10]{1, 2, 3, 4, 5, 6, 7, 8, 9};
    v = f3.gsl()->function(1.2, f3.gsl()->params);
    REQUIRE(v == 1.44);
}

TEST_CASE("integral_qng")
{
    int status;
    double result, abserr;
    size_t neval;

    integral::qng q([](double x) { return pow(x, 2.6) * log(1 / x); },
                    1e-1,
                    0.0);
    status = q(0, 1, &result, &abserr, &neval);
    REQUIRE(status == GSL_SUCCESS);
    REQUIRE(__D_EQ9(result, 7.716049379303083211E-02));
    REQUIRE(__D_EQ9(abserr, 9.424302199601294244E-08));
    REQUIRE(neval == 21);

    status = q(1, 0, &result, &abserr, &neval);
    REQUIRE(status == GSL_SUCCESS);
    REQUIRE(__D_EQ9(result, -7.716049379303083211E-02));
    REQUIRE(__D_EQ9(abserr, 9.424302199601294244E-08));
    REQUIRE(neval == 21);

    status = q(0, 1, &result);
    REQUIRE(status == GSL_SUCCESS);
    REQUIRE(__D_EQ9(result, 7.716049379303083211E-02));
}
