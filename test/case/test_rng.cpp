#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <rand/rand.h>
#include <rand/rng.h>
#include <test_util.h>

using namespace std;
using namespace iexp;

TEST_CASE("test_rng")
{
    rand::rng r(rand::BOROSH13, 1);
    REQUIRE(strcmp("borosh13", r.name()) == 0);
#if 0
    cout << "ulong: " << r.uniform_ulong() << endl;
    cout << "ulong < 2: " << r.uniform_ulong(2) << endl;
    cout << "double: " << r.uniform_double() << endl;
    cout << "pos double: " << r.uniform_pos_double() << endl;
#endif

    rand::rng r2(rand::ZUF, 1);
    REQUIRE(strcmp("zuf", r2.name()) == 0);
#if 0
    cout << "ulong: " << r2.uniform_ulong() << endl;
    cout << "ulong < 2: " << r2.uniform_ulong(2) << endl;
    cout << "double: " << r2.uniform_double() << endl;
    cout << "pos double: " << r2.uniform_pos_double() << endl;
#endif

    rand::rng r3;
    REQUIRE(strcmp("mt19937", r3.name()) == 0);
}

TEST_CASE("test_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::rand(v);
    REQUIRE(&vr == &v);

    v2 = rand::rand(v) + rand::rand(v);

    v2 = rand::rand(v, 999) + rand::rand(v, 12345, rand::ZUF);

    iexp::Matrix3i w(3, 3), w2(3, 3);
    rand::rng r;
    w2 = rand::rand(w, r) + rand::rand(w, 12345, rand::ZUF) + rand::rand(w, r);
    // cout << w2;

    iexp::Matrix3i &wr = rand::rand(w);
    REQUIRE(&wr == &w);
}
