#include <catch.hpp>
#include <fmin/fmin.h>
#include <fmin/mfmin.h>
#include <iostream>
#include <math/constant.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_fmin")
{
    for (int i = 0; i < 3; ++i) {
        iexp::fmin m([](double x, void *) { return pow(x, 4) - 1; },
                     -3,
                     17,
                     -1,
                     nullptr,
                     (iexp::fmin::type)i);
        double r = m.find(0.001, 0.001);
        REQUIRE(__D_EQ3(r, 0));

        iexp::fmin m2([](double x, void *) { return cos(x); },
                      0,
                      6,
                      3,
                      nullptr,
                      (iexp::fmin::type)i);
        r = m2.find(0.001, 0.001);
        REQUIRE(__D_EQ3(r, M_PI));
    }
}

double fff(const VectorXd &x, void *opaque)
{
    double u1 = x(0);
    double u2 = x(1);
    double u3 = x(2);
    double u4 = x(3);

    double t1 = u1 * u1 - u2;
    double t2 = u3 * u3 - u4;
    return 100 * t1 * t1 + (1 - u1) * (1 - u1) + 90 * t2 * t2 +
           (1 - u3) * (1 - u3) +
           10.1 * ((1 - u2) * (1 - u2) + (1 - u4) * (1 - u4)) +
           19.8 * (1 - u2) * (1 - u4);
}

TEST_CASE("test_mfmin")
{
    for (int i = 0; i < 3; ++i) {
        VectorXd x(4), step(4);
        x << -3, -1, -3, -1;
        step.setOnes();
        iexp::mfmin m(fff, x, step, nullptr, (iexp::mfmin::type)i);

        VectorXd r = m.find(1e-3, 5000);
        std::cout << fff(r, nullptr) << std::endl;
        REQUIRE(__D_EQ3(fff(r, nullptr), 0));
    }
}
