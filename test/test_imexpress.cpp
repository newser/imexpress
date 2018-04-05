#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

#include <common/init.h>

int main(int argc, char *argv[])
{
    // global setup...
    iexp::init();

    int result = Catch::Session().run(argc, argv);

    // global clean-up...
    iexp::exit();

    return result;
}
