#include <catch.hpp>
#include <iostream>
#include <math/constant.h>
#include <rand/beta.h>
#include <rand/bi_gauss.h>
#include <rand/binomial.h>
#include <rand/cauchy.h>
#include <rand/chisq.h>
#include <rand/choose.h>
#include <rand/dirichlet.h>
#include <rand/discrete.h>
#include <rand/exp.h>
#include <rand/expow.h>
#include <rand/f.h>
#include <rand/flat.h>
#include <rand/gamma.h>
#include <rand/gauss.h>
#include <rand/gauss_tail.h>
#include <rand/geometric.h>
#include <rand/gumbel1.h>
#include <rand/gumbel2.h>
#include <rand/hyper_geometric.h>
#include <rand/landau.h>
#include <rand/landau.h>
#include <rand/laplace.h>
#include <rand/levy.h>
#include <rand/levy_skew.h>
#include <rand/lgnorm.h>
#include <rand/log.h>
#include <rand/logistic.h>
#include <rand/mul_gauss.h>
#include <rand/mul_nomial.h>
#include <rand/neg_binomial.h>
#include <rand/pareto.h>
#include <rand/pascal.h>
#include <rand/poisson.h>
#include <rand/qrand.h>
#include <rand/qrng.h>
#include <rand/rand.h>
#include <rand/rayleigh.h>
#include <rand/rayleigh_tail.h>
#include <rand/sample.h>
#include <rand/shuffle.h>
#include <rand/spherical.h>
#include <rand/t.h>
#include <rand/weibull.h>
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

    iexp::VectorXd &vr = rand::gauss_rand(v);
    REQUIRE(&vr == &v);

    v2 = rand::gauss_rand(v, 2.0) +
         rand::gauss_rand(v, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gauss_rand(w);
    w2 = rand::gauss_rand(w) + rand::gauss_rand(w, 2.0) +
         rand::gauss_rand(w, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &wr = rand::gauss_rand(w);
    REQUIRE(&wr == &w);

    iexp::MatrixXd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::gauss_rand(cw, 2.0);
    cw2 = rand::gauss_rand(cw, 2.0) +
          rand::gauss_rand(cw, 3.0, 1234, rand::BOROSH13) +
          rand::gauss_rand(cw, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &cwr = rand::gauss_rand(cw, 3.0, 1234, rand::BOROSH13);
    REQUIRE(&cwr == &cw);

#if 0 // #ifdef IEXP_MGL2
    iexp::VectorXd vv(100);
    rand::gauss_rand(vv);
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

    iexp::VectorXd &vr = rand::gausst_rand(v, 12);
    REQUIRE(&vr == &v);

    v2 = rand::gausst_rand(v, 35, 2.0) +
         rand::gausst_rand(v, 23, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gausst_rand(w, 1);
    w2 = rand::gausst_rand(w, 2) + rand::gausst_rand(w, 3, 2.0) +
         rand::gausst_rand(w, 4, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &wr = rand::gausst_rand(w, 5);
    REQUIRE(&wr == &w);

    iexp::MatrixXd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::gausst_rand(cw, 6, 2.0);
    cw2 = rand::gausst_rand(cw, 7, 2.0) +
          rand::gausst_rand(cw, 8, 3.0, 1234, rand::BOROSH13) +
          rand::gausst_rand(cw, 9, 3.0, 1234, rand::BOROSH13);

    iexp::MatrixXd &cwr = rand::gausst_rand(cw, 10, 3.0, 1234, rand::BOROSH13);
    REQUIRE(&cwr == &cw);

#if 0 // #ifdef IEXP_MGL2
    iexp::VectorXd vv(100);
    rand::gausst_rand(vv, 3.0);
    mglData y(100);
    y.Link(vv.data(), vv.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("normalt_rand.png");
#endif
}

TEST_CASE("test_bi_gauss")
{
    iexp::VectorXcd v(10), v2(10);

    iexp::VectorXcd &vr = rand::bgauss_rand(v, 1.0, 1.0, 0);
    REQUIRE(&vr == &v);

    v2 = rand::bgauss_rand(v, 1.0, 1.0, 0) + rand::bgauss_rand(v, 1.0, 1.0, 0);

    iexp::MatrixXcd m(10, 8), m2(10, 8);

    iexp::MatrixXcd &mr = rand::bgauss_rand(m, 1.0, 1.0, 0);
    REQUIRE(&mr == &m);

    m2 = rand::bgauss_rand(m, 1.0, 2.0, 0) + rand::bgauss_rand(m, 3.0, 4.0, 0);

#if 0 // #ifdef IEXP_MGL2
    iexp::VectorXcd vv(100);
    rand::bgauss_rand(vv, 1.0, 1.0, 0.9);
    iexp::VectorXd vv1 = vv.real();
    iexp::VectorXd vv2 = vv.imag();

    mglData x(100), y(100);
    x.Link(vv1.data(), vv1.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-3, 3, -3, 3);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("bi_norm_rand.png");
#endif
}

TEST_CASE("test_mul_normal")
{
    double mu[2] = {1, 2};
    double L[4] = {4, 2, 2, 3};
    rand::mgauss_rng r(2, mu, L);

    Matrix<double, 2, 100> result;
    for (Index i = 0; i < 100; ++i) {
        r.next(result.data() + i * 2);
    }

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1 = result.row(0);
    VectorXd v2 = result.row(1);
    
    mglData x(100), y(100);
    x.Link(v1.data(), v1.size());
    y.Link(v2.data(), v2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-3, 3, -3, 3);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("mul_norm_rand.png");
#endif
}

TEST_CASE("test_exp_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::exp_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::exp_rand(v, 2.0) + rand::exp_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::exp_rand(w, 2.0);
    w2 = rand::exp_rand(w, 3.3) + rand::exp_rand(w, 4.4) +
         rand::exp_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::exp_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::exp_rand(v1, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -1, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("exp_rand.png");
#endif
}

TEST_CASE("test_laplace_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::laplace_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::laplace_rand(v, 2.0) + rand::laplace_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::laplace_rand(w, 2.0);
    w2 = rand::laplace_rand(w, 3.3) + rand::laplace_rand(w, 4.4) +
         rand::laplace_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::laplace_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::laplace_rand(v1, 1.0);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -10, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("laplace_rand.png");
#endif
}

TEST_CASE("test_expow_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::expow_rand(v, 1.0, 2.0);
    REQUIRE(&vr == &v);

    v2 = rand::expow_rand(v, 2.0, 100) + rand::expow_rand(v, 3.0, 4.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::expow_rand(w, 2.0, 3.3);
    w2 = rand::expow_rand(w, 3.3, 3.3) + rand::expow_rand(w, 4.4, 3.3) +
         rand::expow_rand(w, 5.5, 3.3);

    iexp::MatrixXd &wr = rand::expow_rand(w, 99, 3.3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::expow_rand(v1, 1.0, 2.5);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -10, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("expow_rand.png");
#endif
}

TEST_CASE("test_cauchy_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::cauchy_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::cauchy_rand(v, 2.0) + rand::cauchy_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::cauchy_rand(w, 2.0);
    w2 = rand::cauchy_rand(w, 3.3) + rand::cauchy_rand(w, 4.4) +
         rand::cauchy_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::cauchy_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::cauchy_rand(v1, 1.0);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("cauchy_rand.png");
#endif
}

TEST_CASE("test_rayleigh_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::rayl_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::rayl_rand(v, 2.0) + rand::rayl_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::rayl_rand(w, 2.0);
    w2 = rand::rayl_rand(w, 3.3) + rand::rayl_rand(w, 4.4) +
         rand::rayl_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::rayl_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::rayl_rand(v1, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("rayleigh_rand.png");
#endif
}

TEST_CASE("test_rayleigh_tail_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::raylt_rand(v, 1.0, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::raylt_rand(v, 1.0, 2.0) + rand::raylt_rand(v, 1.0, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::raylt_rand(w, 1.0, 2.0);
    w2 = rand::raylt_rand(w, 1.0, 3.3) + rand::raylt_rand(w, 1.0, 4.4) +
         rand::raylt_rand(w, 1.0, 5.5);

    iexp::MatrixXd &wr = rand::raylt_rand(w, 1.0, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::raylt_rand(v1, 1.5, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("rayleigh_tail_rand.png");
#endif
}

TEST_CASE("test_landau_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::landau_rand(v);
    REQUIRE(&vr == &v);

    v2 = rand::landau_rand(v) + rand::landau_rand(v);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::landau_rand(w);
    w2 = rand::landau_rand(w) + rand::landau_rand(w) + rand::landau_rand(w);

    iexp::MatrixXd &wr = rand::landau_rand(w);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::landau_rand(v1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("landau_rand.png");
#endif
}

TEST_CASE("test_levy_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::levy_rand(v, 1, 1);
    REQUIRE(&vr == &v);

    v2 = rand::levy_rand(v, 2, 3) + rand::levy_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::levy_rand(w, 2, 3);
    w2 = rand::levy_rand(w, 2, 3) + rand::levy_rand(w, 2, 3) +
         rand::levy_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::levy_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::levy_rand(v1, 1, 1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("levy_rand.png");
#endif
}

TEST_CASE("test_levy_skew_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::levysk_rand(v, 1, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::levysk_rand(v, 2, 3, 2) + rand::levysk_rand(v, 2, 4, 2);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::levysk_rand(w, 2, 3, 2);
    w2 = rand::levysk_rand(w, 2, 3, 2) + rand::levysk_rand(w, 2, 3, 2) +
         rand::levysk_rand(w, 2, 3, 2);

    iexp::MatrixXd &wr = rand::levysk_rand(w, 2, 3, 2);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::levysk_rand(v1, 1, 1, 1);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -5, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("levysk_rand.png");
#endif
}

TEST_CASE("test_flat_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::flat_rand(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::flat_rand(v, 2, 3) + rand::flat_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::flat_rand(w, 2, 3);
    w2 = rand::flat_rand(w, 2, 3) + rand::flat_rand(w, 2, 3) +
         rand::flat_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::flat_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::flat_rand(v1, 1, 5);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 6);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("flat_rand.png");
#endif
}

TEST_CASE("test_gamma_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::gamma_rand(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::gamma_rand(v, 2, 3) + rand::gamma_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gamma_rand(w, 2, 3);
    w2 = rand::gamma_rand(w, 2, 3) + rand::gamma_rand(w, 2, 3) +
         rand::flat_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::gamma_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::gamma_rand(v1, 2, 1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 6);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("gamma_rand.png");
#endif
}

TEST_CASE("test_beta_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::beta_rand(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::beta_rand(v, 2, 3) + rand::beta_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::beta_rand(w, 2, 3);
    w2 = rand::beta_rand(w, 2, 3) + rand::beta_rand(w, 2, 3) +
         rand::flat_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::beta_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::beta_rand(v1, 2, 2);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 1);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("beta_rand.png");
#endif
}

TEST_CASE("test_lgnorm_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::lgnorm_rand(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::lgnorm_rand(v, 2, 3) + rand::lgnorm_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::lgnorm_rand(w, 2, 3);
    w2 = rand::lgnorm_rand(w, 2, 3) + rand::lgnorm_rand(w, 2, 3) +
         rand::flat_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::lgnorm_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::lgnorm_rand(v1, 0, 1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 1);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("lgnorm_rand.png");
#endif
}

TEST_CASE("test_chisq_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::chisq_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::chisq_rand(v, 2.0) + rand::chisq_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::chisq_rand(w, 2.0);
    w2 = rand::chisq_rand(w, 3.3) + rand::chisq_rand(w, 4.4) +
         rand::chisq_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::chisq_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef Ichisq_MGL2
    VectorXd v1(100);
    rand::chisq_rand(v1, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -1, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("chisq_rand.png");
#endif
}

TEST_CASE("test_f_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::f_rand(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::f_rand(v, 2, 3) + rand::f_rand(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::f_rand(w, 2, 3);
    w2 = rand::f_rand(w, 2, 3) + rand::f_rand(w, 2, 3) +
         rand::flat_rand(w, 2, 3);

    iexp::MatrixXd &wr = rand::f_rand(w, 2, 3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::f_rand(v1, 1, 1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 1);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("f_rand.png");
#endif
}

TEST_CASE("test_t_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::t_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::t_rand(v, 2.0) + rand::t_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::t_rand(w, 2.0);
    w2 = rand::t_rand(w, 3.3) + rand::t_rand(w, 4.4) + rand::t_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::t_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef It_MGL2
    VectorXd v1(100);
    rand::t_rand(v1, 5.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -4, 4);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("t_rand.png");
#endif
}

TEST_CASE("test_lgst_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::lgst_rand(v, 1.0);
    REQUIRE(&vr == &v);

    v2 = rand::lgst_rand(v, 2.0) + rand::lgst_rand(v, 3.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::lgst_rand(w, 2.0);
    w2 = rand::lgst_rand(w, 3.3) + rand::lgst_rand(w, 4.4) +
         rand::lgst_rand(w, 5.5);

    iexp::MatrixXd &wr = rand::lgst_rand(w, 99);
    REQUIRE(&wr == &w);

#if 0 // #ifdef It_MGL2
    VectorXd v1(100);
    rand::lgst_rand(v1, 1);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, -4, 4);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("lgst_rand.png");
#endif
}

TEST_CASE("test_pareto_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::pareto_rand(v, 1.0, 2.0);
    REQUIRE(&vr == &v);

    v2 = rand::pareto_rand(v, 2.0, 100) + rand::pareto_rand(v, 3.0, 4.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::pareto_rand(w, 2.0, 3.3);
    w2 = rand::pareto_rand(w, 3.3, 3.3) + rand::pareto_rand(w, 4.4, 3.3) +
         rand::pareto_rand(w, 5.5, 3.3);

    iexp::MatrixXd &wr = rand::pareto_rand(w, 99, 3.3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::pareto_rand(v1, 3.0, 2.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("pareto_rand.png");
#endif
}

TEST_CASE("test_sph2_rand")
{
    iexp::VectorXcd v(10), v2(10);

    iexp::VectorXcd &vr = rand::sph2_rand(v);
    REQUIRE(&vr == &v);

    v2 = rand::sph2_rand(v) + rand::sph2_rand(v);

    iexp::MatrixXcd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::sph2_rand<true>(w);
    w2 = rand::sph2_rand(w) + rand::sph2_rand(w) + rand::sph2_rand(w);

    iexp::MatrixXcd &wr = rand::sph2_rand(w);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXcd v1(100);
    rand::sph2_rand(v1);
    VectorXd r = v1.real();
    VectorXd i = v1.imag();
    
    mglData x(100), y(100);
    x.Link(r.data(), r.size());
    y.Link(i.data(), i.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-1.5, 1.5, -1.5, 1.5);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("sph2_rand.png");
#endif
}

TEST_CASE("test_wbl_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::wbl_rand(v, 1.0, 2.0);
    REQUIRE(&vr == &v);

    v2 = rand::wbl_rand(v, 2.0, 100) + rand::wbl_rand(v, 3.0, 4.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::wbl_rand(w, 2.0, 3.3);
    w2 = rand::wbl_rand(w, 3.3, 3.3) + rand::wbl_rand(w, 4.4, 3.3) +
         rand::wbl_rand(w, 5.5, 3.3);

    iexp::MatrixXd &wr = rand::wbl_rand(w, 99, 3.3);
    REQUIRE(&wr == &w);

#if o // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::wbl_rand(v1, 1.0, 2.0);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 3);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("wbl_rand.png");
#endif
}

TEST_CASE("test_gbl1_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::gbl1_rand(v, 1.0, 2.0);
    REQUIRE(&vr == &v);

    v2 = rand::gbl1_rand(v, 2.0, 100) + rand::gbl1_rand(v, 3.0, 4.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gbl1_rand(w, 2.0, 3.3);
    w2 = rand::gbl1_rand(w, 3.3, 3.3) + rand::gbl1_rand(w, 4.4, 3.3) +
         rand::gbl1_rand(w, 5.5, 3.3);

    iexp::MatrixXd &wr = rand::gbl1_rand(w, 99, 3.3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::gbl1_rand(v1, 1.0, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("gbl1_rand.png");
#endif
}

TEST_CASE("test_gbl2_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::gbl2_rand(v, 1.0, 2.0);
    REQUIRE(&vr == &v);

    v2 = rand::gbl2_rand(v, 2.0, 100) + rand::gbl2_rand(v, 3.0, 4.0);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gbl2_rand(w, 2.0, 3.3);
    w2 = rand::gbl2_rand(w, 3.3, 3.3) + rand::gbl2_rand(w, 4.4, 3.3) +
         rand::gbl2_rand(w, 5.5, 3.3);

    iexp::MatrixXd &wr = rand::gbl2_rand(w, 99, 3.3);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::gbl2_rand(v1, 1.0, 1.0);
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("gbl2_rand.png");
#endif
}

TEST_CASE("test_drch_rand")
{
    double alpha[2] = {2.5, 5.5};
    iexp::VectorXcd v(10), v2(10);

    iexp::VectorXcd &vr = rand::drch_rand(v, 2, alpha);
    REQUIRE(&vr == &v);

    v2 = rand::drch_rand(v, 2, alpha) + rand::drch_rand(v, 2, alpha);

    iexp::MatrixXcd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::drch_rand(w, 2, alpha);
    w2 = rand::drch_rand(w, 2, alpha) + rand::drch_rand(w, 2, alpha) +
         rand::drch_rand(w, 2, alpha);

    iexp::MatrixXcd &wr = rand::drch_rand(w, 2, alpha);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXcd v1(100);
    rand::drch_rand(v1, 2, alpha);
    VectorXd r = v1.real();
    VectorXd i = v1.imag();
    
    mglData x(100), y(100);
    x.Link(r.data(), r.size());
    y.Link(i.data(), i.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-1.5, 1.5, -1.5, 1.5);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("drch_rand.png");
#endif
}

TEST_CASE("test_discrete_rand")
{
    const double p[4] = {0.1, 0.2, 0.3, 0.4};
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::discrete_rand(v, 4, p);
    REQUIRE(&vr == &v);

    v2 = rand::discrete_rand(v, 4, p) + rand::discrete_rand(v, 4, p);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::discrete_rand(w, 4, p);
    w2 = rand::discrete_rand(w, 4, p) + rand::discrete_rand(w, 4, p) +
         rand::discrete_rand(w, 4, p);

    iexp::MatrixXi &wr = rand::discrete_rand(w, 4, p);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::discrete_rand(vi, 4, p);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("discrete_rand.png");
#endif
}

TEST_CASE("test_poiss_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::poiss_rand(v, 2.5);
    REQUIRE(&vr == &v);

    v2 = rand::poiss_rand(v, 2.5) + rand::poiss_rand(v, 3.3);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::poiss_rand(w, 4.4);
    w2 = rand::poiss_rand(w, 5.5) + rand::poiss_rand(w, 7.7) +
         rand::poiss_rand(w, 6.6);

    iexp::MatrixXi &wr = rand::poiss_rand(w, 8.8);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::poiss_rand(vi, 2.5);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("poiss_rand.png");
#endif
}

TEST_CASE("test_bnom_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::bnom_rand(v, 0.5, 9);
    REQUIRE(&vr == &v);

    v2 = rand::bnom_rand(v, 0.5, 9) + rand::bnom_rand(v, 0.5, 9);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::bnom_rand(w, 0.5, 9);
    w2 = rand::bnom_rand(w, 0.5, 9) + rand::bnom_rand(w, 0.5, 9) +
         rand::bnom_rand(w, 0.5, 9);

    iexp::MatrixXi &wr = rand::bnom_rand(w, 0.5, 9);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::bnom_rand(vi, 0.5, 9);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("bnom_rand.png");
#endif
}

TEST_CASE("test_nbnom_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::nbnom_rand(v, 0.5, 9);
    REQUIRE(&vr == &v);

    v2 = rand::nbnom_rand(v, 0.5, 9) + rand::nbnom_rand(v, 0.5, 9);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::nbnom_rand(w, 0.5, 9);
    w2 = rand::nbnom_rand(w, 0.5, 9) + rand::nbnom_rand(w, 0.5, 9) +
         rand::nbnom_rand(w, 0.5, 9);

    iexp::MatrixXi &wr = rand::nbnom_rand(w, 0.5, 9);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::nbnom_rand(vi, 0.5, 9);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("nbnom_rand.png");
#endif
}

TEST_CASE("test_multinomial")
{
    double p[2] = {1, 2};
    rand::mnom_rng r(2, p, 6);

    Matrix<unsigned int, 2, 100> result;
    for (Index i = 0; i < 100; ++i) {
        r.next(result.data() + i * 2);
    }

#if 0 // #ifdef IEXP_MGL2
    VectorXd v1 = result.row(0).cast<double>();
    VectorXd v2 = result.row(1).cast<double>();

    mglData x(100), y(100);
    x.Link(v1.data(), v1.size());
    y.Link(v2.data(), v2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 6, 0, 6);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("multinomial.png");
#endif
}

TEST_CASE("test_pascal_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::pascal_rand(v, 0.5, 9);
    REQUIRE(&vr == &v);

    v2 = rand::pascal_rand(v, 0.5, 9) + rand::pascal_rand(v, 0.5, 9);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::pascal_rand(w, 0.5, 9);
    w2 = rand::pascal_rand(w, 0.5, 9) + rand::pascal_rand(w, 0.5, 9) +
         rand::pascal_rand(w, 0.5, 9);

    iexp::MatrixXi &wr = rand::pascal_rand(w, 0.5, 9);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::pascal_rand(vi, 0.5, 3);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 10);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("pascal_rand.png");
#endif
}

TEST_CASE("test_geo_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::geo_rand(v, 0.5);
    REQUIRE(&vr == &v);

    v2 = rand::geo_rand(v, 0.5) + rand::geo_rand(v, 0.5);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::geo_rand(w, 0.5);
    w2 = rand::geo_rand(w, 0.5) + rand::geo_rand(w, 0.5) +
         rand::geo_rand(w, 0.5);

    iexp::MatrixXi &wr = rand::geo_rand(w, 0.5);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::geo_rand(vi, 0.5);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("geo_rand.png");
#endif
}

TEST_CASE("test_hgeo_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::hgeo_rand(v, 3, 10, 5);
    REQUIRE(&vr == &v);

    v2 = rand::hgeo_rand(v, 3, 10, 5) + rand::hgeo_rand(v, 3, 10, 5);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::hgeo_rand(w, 3, 10, 5);
    w2 = rand::hgeo_rand(w, 3, 10, 5) + rand::hgeo_rand(w, 3, 10, 5) +
         rand::hgeo_rand(w, 3, 10, 5);

    iexp::MatrixXi &wr = rand::hgeo_rand(w, 3, 10, 5);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::hgeo_rand(vi, 5, 20, 3);
    VectorXd v1 = vi.cast<double>();

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("hgeo_rand.png");
#endif
}

TEST_CASE("test_log_rand")
{
    iexp::VectorXi v(10), v2(10);

    iexp::VectorXi &vr = rand::log_rand(v, 0.7);
    REQUIRE(&vr == &v);

    v2 = rand::log_rand(v, 0.7) + rand::log_rand(v, 0.7);

    iexp::MatrixXi w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::log_rand(w, 0.7);
    w2 = rand::log_rand(w, 0.7) + rand::log_rand(w, 0.7) +
         rand::log_rand(w, 0.7);

    iexp::MatrixXi &wr = rand::log_rand(w, 0.7);
    REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
    VectorXi vi(100);
    rand::log_rand(vi, 0.7);
    VectorXd v1 = vi.cast<double>();
    
    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 5);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("log_rand.png");
#endif
}

TEST_CASE("test_shuffle")
{
    iexp::VectorXi v(10), v2(10);

    v = VectorXi::LinSpaced(10, 0, 9);
    iexp::VectorXi &vr = rand::shuffle(v);
    REQUIRE(&vr == &v);

    v2 = rand::shuffle(v, 0.7) + rand::shuffle(v);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::shuffle(w);
    w2 = rand::shuffle(w) + rand::shuffle(w) + rand::shuffle(w);

    iexp::MatrixXd &wr = rand::shuffle(w);
    REQUIRE(&wr == &w);
}

TEST_CASE("test_choose")
{
    iexp::VectorXi v(10), v2(10);

    v = VectorXi::LinSpaced(10, 0, 9);
    iexp::VectorXi vr = rand::choose(v.array(), 5).matrix();
    REQUIRE(vr.size() == 5);

    v2 = rand::choose(rand::shuffle(v).array() + rand::shuffle(v).array(), 7);
    REQUIRE(v2.size() == 7);
}

TEST_CASE("test_sample")
{
    iexp::VectorXd v(10), v2(10);

    v = VectorXd::LinSpaced(10, 0, 9);
    iexp::VectorXd vr = rand::sample(v.array(), 5).matrix();
    REQUIRE(vr.size() == 5);

    v2 = rand::sample(rand::shuffle(v).array() + rand::shuffle(v).array(), 7);
    REQUIRE(v2.size() == 7);
}
