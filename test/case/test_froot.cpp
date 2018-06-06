#include <catch.hpp>
#include <froot/dfroot.h>
#include <froot/froot.h>
#include <iostream>
#include <math/constant.h>
#include <test_util.h>

using namespace iexp;

TEST_CASE("test_froot")
{
    // 0 root
    for (int i = 0; i < 3; ++i) {
        try {
            // root: 1, 2
            froot
            fr([](double x)
                   -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
               0.1,
               0.99,
               (froot::type)i);
            REQUIRE(false);
        } catch (...) {
        }
    }

    // 1 root
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        froot fr([](double x)
                     -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                 0.1,
                 1.2,
                 (froot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find<false>(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find(1e-6);
        REQUIRE(__D_EQ6(rt, 1));
    }

    // 2 root
    for (int i = 1; i < 3; ++i) {
        try {
            // root: 1, 2
            froot
            fr([](double x)
                   -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
               0.1,
               2.2,
               (froot::type)i);
            REQUIRE(false);
        } catch (...) {
        }
    }

    // multiple root in interval
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        froot fr([](double x)
                     -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                 0.1,
                 3.2,
                 (froot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find<false>(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find(1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));
    }
}

TEST_CASE("test_dfroot")
{
    // 1 root
    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        dfroot fr([](double x)
                      -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                  [](double x) -> double { return 3 * x * x - 12 * x + 11; },
                  [](double x, double &f, double &df) -> void {
                      f = x * x * x - 6 * x * x + 11 * x - 6;
                      df = 3 * x * x - 12 * x + 11;
                  },
                  -0.1,
                  (dfroot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE(__D_EQ6(rt, 1));

        rt = fr.find(1e-6);
        REQUIRE(__D_EQ6(rt, 1));
    }

    for (int i = 0; i < 3; ++i) {
        // root: 1, 2
        dfroot fr([](double x)
                      -> double { return x * x * x - 6 * x * x + 11 * x - 6; },
                  [](double x) -> double { return 3 * x * x - 12 * x + 11; },
                  [](double x, double &f, double &df) -> void {
                      f = x * x * x - 6 * x * x + 11 * x - 6;
                      df = 3 * x * x - 12 * x + 11;
                  },
                  3.2,
                  (dfroot::type)i);

        double rt;

        rt = fr.find(0, 1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));

        rt = fr.find(1e-6);
        REQUIRE((__D_EQ6(rt, 1) || __D_EQ6(rt, 2) || __D_EQ6(rt, 3)));
    }
}
