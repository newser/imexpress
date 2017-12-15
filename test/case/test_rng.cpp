#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <rand/normal.h>
#include <rand/normal_tail.h>
#include <rand/qrand.h>
#include <rand/qrng.h>
#include <rand/rand.h>
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

    rand::rng r4(std::move(r3));
    rand::rng r5 = std::move(r4);
    r5 = rand::rng();

    rand::rng r6;
    r6 = r5;
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

TEST_CASE("test_qrng")
{
    double x[10];

    rand::qrng r(rand::NIEDERREITER_2, 2);
    REQUIRE(strcmp("niederreiter-base-2", r.name()) == 0);
    r.next(x);

    rand::qrng r2(rand::REVERSEHALTON, 3);
    REQUIRE(strcmp("reversehalton", r2.name()) == 0);
    r.next(x);

    rand::qrng r3(rand::SOBOL, 10);
    REQUIRE(strcmp("sobol", r3.name()) == 0);

    rand::qrng r4(std::move(r3));
    rand::qrng r5 = std::move(r4);
    r5 = rand::qrng(rand::SOBOL, 10);

    rand::qrng r6(rand::SOBOL, 10);
    r6 = r5;
}

TEST_CASE("test_qrand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::qrand(v, rand::SOBOL);
    REQUIRE(&vr == &v);

    v2 = rand::qrand(v, rand::SOBOL) + rand::qrand(v, rand::HALTON);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::qrand(w, rand::REVERSEHALTON);
    w2 = rand::qrand(w, rand::REVERSEHALTON) + rand::qrand(w, rand::HALTON) +
         rand::qrand(w, rand::SOBOL);

    iexp::MatrixXd &wr = rand::qrand(w, rand::HALTON);
    REQUIRE(&wr == &w);

    iexp::MatrixXcd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::qrand(cw, rand::REVERSEHALTON);
    cw2 = rand::qrand(cw, rand::REVERSEHALTON) + rand::qrand(cw, rand::HALTON) +
          rand::qrand(cw, rand::SOBOL);

    iexp::MatrixXcd &cwr = rand::qrand(cw, rand::HALTON);
    REQUIRE(&cwr == &cw);
}

TEST_CASE("test_normal_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::norm_rand(v);
    REQUIRE(&vr == &v);

    v2 =
        rand::norm_rand(v, 2.0) + rand::norm_rand(v, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::norm_rand(w);
    w2 = rand::norm_rand(w) + rand::norm_rand(w, 2.0) +
         rand::norm_rand(w, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &wr = rand::norm_rand(w);
    REQUIRE(&wr == &w);

    iexp::MatrixXd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::norm_rand(cw, 2.0);
    cw2 = rand::norm_rand(cw, 2.0) +
          rand::norm_rand(cw, 3.0, 1234, rand::BOROSH13) +
          rand::norm_rand(cw, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &cwr = rand::norm_rand(cw, 3.0, 1234, rand::BOROSH13);
    REQUIRE(&cwr == &cw);

#ifdef IEXP_MGL2
    iexp::VectorXd vv(100);
    rand::norm_rand(vv);
    mglData y(100);
    y.Link(vv.data(), vv.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 2);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("normal_rand.png");
#endif
}

TEST_CASE("test_normal_tail_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::normt_rand(v, 12);
    REQUIRE(&vr == &v);

    v2 = rand::normt_rand(v, 35, 2.0) +
         rand::normt_rand(v, 23, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::normt_rand(w, 1);
    w2 = rand::normt_rand(w, 2) + rand::normt_rand(w, 3, 2.0) +
         rand::normt_rand(w, 4, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &wr = rand::normt_rand(w, 5);
    REQUIRE(&wr == &w);

    iexp::MatrixXd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::normt_rand(cw, 6, 2.0);
    cw2 = rand::normt_rand(cw, 7, 2.0) +
          rand::normt_rand(cw, 8, 3.0, 1234, rand::BOROSH13) +
          rand::normt_rand(cw, 9, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &cwr = rand::normt_rand(cw, 10, 3.0, 1234, rand::BOROSH13);
    REQUIRE(&cwr == &cw);
}
