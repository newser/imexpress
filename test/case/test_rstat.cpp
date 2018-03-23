#include <../test/test_util.h>
#include <catch.hpp>
#include <rstat/rstat.h>

using namespace iexp;

TEST_CASE("rstat")
{
    rstat rs;
    REQUIRE(rs.size() == 0);

    rs.add(17.2);
    rs.add(18.1);
    rs.add(16.5);
    rs.add(18.3);
    rs.add(12.6);
    REQUIRE(rs.size() == 5);
    REQUIRE(__D_EQ_IN(rs.mean(), 16.54, 1e-9));
    REQUIRE(__D_EQ_IN(rs.var(), 5.373, 1e-9));
    REQUIRE(__D_EQ_IN(rs.max(), 18.3, 1e-9));
    REQUIRE(__D_EQ_IN(rs.min(), 12.6, 1e-9));
    REQUIRE(__D_EQ_IN(rs.median(), 16.5, 1e-9));
    REQUIRE(__D_EQ_IN(rs.std(), 2.3179732526498227, 1e-9));
    REQUIRE(__D_EQ_IN(rs.rms(), 16.669433103738111, 1e-9));
    REQUIRE(__D_EQ_IN(rs.std_mean(), 1.0366291525902596, 1e-9));
    REQUIRE(__D_EQ_IN(rs.skewness(), -0.82905750003696543, 1e-9));
    REQUIRE(__D_EQ_IN(rs.kurtosis(), -1.2217029020861698, 1e-9));

    rs.reset();
    REQUIRE(rs.size() == 0);
}
