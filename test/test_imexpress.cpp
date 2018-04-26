#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <common/init.h>

int main(int argc, char *argv[])
{
    // global setup...
    iexp::init();

    const char *__argv[2] = {"", "test_sph2_rand"};
    argc = sizeof(__argv) / sizeof(const char *);
    argv = (char **)__argv;
    int result = Catch::Session().run(argc, argv);

    // global clean-up...
    iexp::exit();

    return result;
}
