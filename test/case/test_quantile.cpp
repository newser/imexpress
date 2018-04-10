#include <catch.hpp>
#include <iostream>
#include <rstat/quantile.h>
#include <sort/sort.h>
#include <stats/quantile.h>

using namespace iexp;

TEST_CASE("quantile")
{
    ArrayXd a = ArrayXd::Random(100);

    quantile q(0.66);
    for (int i = 0; i < a.size(); ++i) {
        q.add(a[i]);
    }

    double v1 = stats::quantile(sort(a), 0.66);
    double v2 = q.get();
}
