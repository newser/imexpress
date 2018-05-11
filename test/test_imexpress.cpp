#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <common/init.h>

int main(int argc, char *argv[])
{
    // global setup...
    iexp::init();

    //    const char *__argv[] = {"", "sf_airy", "sf_airy_scaled",
    //    "sf_airy_deriv", "sf_airy_deriv_scaled"};
    const char *__argv[] = {"", "test_ellint_rf"};
    argc = sizeof(__argv) / sizeof(const char *);
    argv = (char **)__argv;
    int result = Catch::Session().run(argc, argv);

    // global clean-up...
    iexp::exit();

    return result;
}
