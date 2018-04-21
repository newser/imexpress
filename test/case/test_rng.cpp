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
    rand::rng r(rand::rng::type::BOROSH13, 1);
    REQUIRE(strcmp("borosh13", r.name()) == 0);
#if 0
    cout << "ulong: " << r.uniform_ulong() << endl;
    cout << "ulong < 2: " << r.uniform_ulong(2) << endl;
    cout << "double: " << r.uniform_double() << endl;
    cout << "pos double: " << r.uniform_pos_double() << endl;
#endif

    rand::rng r2(rand::rng::type::ZUF, 1);
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

    v2 = rand::rand(v, 999) + rand::rand(v, 12345, rand::rng::type::ZUF);

    iexp::Matrix3i w(3, 3), w2(3, 3);
    rand::rng r;
    w2 = rand::rand(w, r) + rand::rand(w, 12345, rand::rng::type::ZUF) +
         rand::rand(w, r);
    // cout << w2;

    iexp::Matrix3i &wr = rand::rand(w);
    REQUIRE(&wr == &w);
}

TEST_CASE("test_qrng")
{
    double x[10];

    rand::qrng r(rand::qrng::type::NIEDERREITER_2, 2);
    REQUIRE(strcmp("niederreiter-base-2", r.name()) == 0);
    r.next(x);

    rand::qrng r2(rand::qrng::type::REVERSEHALTON, 3);
    REQUIRE(strcmp("reversehalton", r2.name()) == 0);
    r.next(x);

    rand::qrng r3(rand::qrng::type::SOBOL, 10);
    REQUIRE(strcmp("sobol", r3.name()) == 0);

    rand::qrng r4(std::move(r3));
    rand::qrng r5 = std::move(r4);
    r5 = rand::qrng(rand::qrng::type::SOBOL, 10);

    rand::qrng r6(rand::qrng::type::SOBOL, 10);
    r6 = r5;
}

TEST_CASE("test_qrand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::qrand(v, rand::qrng::type::SOBOL);
    REQUIRE(&vr == &v);

    v2 = rand::qrand(v, rand::qrng::type::SOBOL) +
         rand::qrand(v, rand::qrng::type::HALTON);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::qrand(w, rand::qrng::type::REVERSEHALTON);
    w2 = rand::qrand(w, rand::qrng::type::REVERSEHALTON) +
         rand::qrand(w, rand::qrng::type::HALTON) +
         rand::qrand(w, rand::qrng::type::SOBOL);

    iexp::MatrixXd &wr = rand::qrand(w, rand::qrng::type::HALTON);
    REQUIRE(&wr == &w);

    iexp::MatrixXcd cw(3, 4), cw2(3, 4);
    cw.fill(9.9999);
    cw2 = rand::qrand(cw, rand::qrng::type::REVERSEHALTON);
    cw2 = rand::qrand(cw, rand::qrng::type::REVERSEHALTON) +
          rand::qrand(cw, rand::qrng::type::HALTON) +
          rand::qrand(cw, rand::qrng::type::SOBOL);

    iexp::MatrixXcd &cwr = rand::qrand(cw, rand::qrng::type::HALTON);
    REQUIRE(&cwr == &cw);
}

TEST_CASE("test_normal_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::gauss::fill(v);
        REQUIRE(&vr == &v);

        v2 = rand::gauss::fill(v, 2.0) +
             rand::gauss::fill(v, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::gauss::fill(w);
        w2 = rand::gauss::fill(w) + rand::gauss::fill(w, 2.0) +
             rand::gauss::fill(w, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd &wr = rand::gauss::fill(w);
        REQUIRE(&wr == &w);

        iexp::MatrixXd cw(3, 4), cw2(3, 4);
        cw.fill(9.9999);
        cw2 = rand::gauss::fill(cw, 2.0);
        cw2 = rand::gauss::fill(cw, 2.0) +
              rand::gauss::fill(cw, 3.0, 1234, rand::rng::type::BOROSH13) +
              rand::gauss::fill(cw, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd &cwr =
            rand::gauss::fill(cw, 3.0, 1234, rand::rng::type::BOROSH13);
        REQUIRE(&cwr == &cw);

#if 1 // #ifdef IEXP_MGL2
        iexp::VectorXd vv(100);
        rand::gauss::fill(vv);
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

    {
        VectorXd v = VectorXd::LinSpaced(10, -5, 5);
        VectorXd v2 = rand::gauss::pdf(v.array(), 2.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::gauss::pdf(m.array(), 2.0);

        // test compile
        v2 =
            rand::gauss::pdf(v.array(), 2.0) + rand::gauss::pdf(v.array(), 2.0);
        m2 = rand::gauss::pdf(m.array() + m.array(), 2.0);

#if 1 // #ifdef IEXP_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::gauss::pdf(vv.array(), 2.0);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("gauss_pdf.png");
#endif
    }
}

TEST_CASE("test_normal_tail_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::gausst::fill(v, 12);
        REQUIRE(&vr == &v);

        v2 = rand::gausst::fill(v, 35, 2.0) +
             rand::gausst::fill(v, 23, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::gausst::fill(w, 1);
        w2 = rand::gausst::fill(w, 2) + rand::gausst::fill(w, 3, 2.0) +
             rand::gausst::fill(w, 4, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd &wr = rand::gausst::fill(w, 5);
        REQUIRE(&wr == &w);

        iexp::MatrixXd cw(3, 4), cw2(3, 4);
        cw.fill(9.9999);
        cw2 = rand::gausst::fill(cw, 6, 2.0);
        cw2 = rand::gausst::fill(cw, 7, 2.0) +
              rand::gausst::fill(cw, 8, 3.0, 1234, rand::rng::type::BOROSH13) +
              rand::gausst::fill(cw, 9, 3.0, 1234, rand::rng::type::BOROSH13);

        iexp::MatrixXd &cwr =
            rand::gausst::fill(cw, 10, 3.0, 1234, rand::rng::type::BOROSH13);
        REQUIRE(&cwr == &cw);

#if 1 // #ifdef IEXP_MGL2
        iexp::VectorXd vv(100);
        rand::gausst::fill(vv, 3.0);
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

    {
        VectorXd v = VectorXd::LinSpaced(10, -5, 5);
        VectorXd v2 = rand::gausst::pdf(v.array(), 1.5, 2.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::gausst::pdf(m.array(), 1.5, 2.0);

        // test compile
        v2 = rand::gausst::pdf(v.array(), 1.5, 2.0) +
             rand::gausst::pdf(v.array(), 1.5, 2.0);
        m2 = rand::gausst::pdf(m.array() + m.array(), 1.5, 2.0);

#if 1 // #ifdef IEXP_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, 0, 5);
        VectorXd vv2 = rand::gausst::pdf(vv.array(), 1.5, 1.0);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 2);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("gausst_pdf.png");
#endif
    }
}

TEST_CASE("test_bi_gauss")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::bgauss::fill(v, 1.0, 1.0, 0);
        REQUIRE(&vr == &v);

        v2 = rand::bgauss::fill(v, 1.0, 1.0, 0) +
             rand::bgauss::fill(v, 1.0, 1.0, 0);

        iexp::MatrixXd m(10, 8), m2(10, 8);

        iexp::MatrixXd &mr = rand::bgauss::fill(m, 1.0, 1.0, 0);
        REQUIRE(&mr == &m);

        m2 = rand::bgauss::fill(m, 1.0, 2.0, 0) +
             rand::bgauss::fill(m, 3.0, 4.0, 0);

        rand::bgauss::fill(v, v2, 1.0, 1.0, 0);

#if 0 // #ifdef IEXP_MGL2
        iexp::VectorXd vv1(100), vv2(100);
        rand::bgauss::fill(vv1, vv2, 1.0, 1.0, 0.9);

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

    {
        VectorXd vv = VectorXd::LinSpaced(10, -2, 2);
        double z[100];

        rand::bgauss::dist bg(1, 1, 0.9);
        for (int i = 0; i < vv.size(); ++i) {
            for (int j = 0; j < vv.size(); ++j) {
                z[i * vv.size() + j] = bg.pdf(vv[i], vv[j]);
            }
        }

#if 0 // #ifdef IEXP_MGL2
        mglData zz(100);
        zz.Create(10, 10);
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                zz.a[i * 10 + j] = bg.pdf(vv[i], vv[j]);
            }
        }

        mglGraph gr;
        // gr.SetOrigin(0, 0);
        // gr.SetRanges(-2, 2, -2, 2, -2, 2);
        gr.Light(true);
        gr.Rotate(50, 60);
        gr.Box();
        gr.Surf(zz);
        gr.WriteFrame("bgauss_pdf.png");
#endif
    }
}

TEST_CASE("test_mul_normal")
{
    {
        double mu[2] = {1, 2};
        double L[4] = {4, 2, 2, 3};
        rand::mgauss::rng r(2, mu, L);

        Matrix<double, 2, 100> result;
#if 0
        for (Index i = 0; i < 100; ++i) {
            r.next(result.data() + i * 2);
        }
#endif
        rand::mgauss::fill(result, 2, mu, L);

#if 1 // #ifdef IEXP_MGL2
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

    {
        double x[2] = {0};
        double mu[2] = {1, 2};
        double cov[4] = {4, 2, 2, 3};

        rand::mgauss::dist mg(2, mu, cov);
        double p = mg.pdf(x);
        REQUIRE(__D_EQ9(p, 0.028294217120391));

        p = mg.lnpdf(x);
        REQUIRE(__D_EQ9(p, log(0.028294217120391)));

        rand::mgauss::rng r(2, mu, cov);
        MatrixXd m(2, 100);
        rand::mgauss::fill(m, r);

        double e_mu[2], e_cov[4];
        rand::mgauss::dist::mean(m.cols(), m.rows(), m.data(), e_mu);
        rand::mgauss::dist::cov(m.cols(), m.rows(), m.data(), e_cov);
//        cout << e_mu[0] << " " << e_mu[1] << endl;
//         cout << e_cov[0] << " " << e_cov[1] << endl;
//         cout << e_cov[2] << " " << e_cov[3] << endl;

#if 0
        for (int i = 0; i < m.cols(); ++i) {
            double s[2];
            r.next(s);
            m(0, i) = s[0];
            m(1, i) = s[1];
        }
#endif
        rand::mgauss::fill(m, r);
        Vector2d v = rand::mgauss::mean(m.array());
        cout << v << endl;

        Matrix2d rcov = rand::mgauss::cov(m);
        cout << rcov << endl;
    }
}

TEST_CASE("test_exp_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::exp::fill(v, 1.0);
        REQUIRE(&vr == &v);

        v2 = rand::exp::fill(v, 2.0) + rand::exp::fill(v, 3.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::exp::fill(w, 2.0);
        w2 = rand::exp::fill(w, 3.3) + rand::exp::fill(w, 4.4) +
             rand::exp::fill(w, 5.5);

        iexp::MatrixXd &wr = rand::exp::fill(w, 99);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::exp::fill(v1, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::exp::pdf(v.array(), 2.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::exp::pdf(m.array(), 2.0);

        // test compile
        v2 = rand::exp::pdf(v.array(), 2.0) + rand::exp::pdf(v.array(), 2.0);
        m2 = rand::exp::pdf(m.array() + m.array(), 2.0);

#if 1 // #ifdef IEXP_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::exp::pdf(vv.array(), 1.0);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("exp_pdf.png");
#endif
    }
}

TEST_CASE("test_laplace_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::laplace::fill(v, 1.0);
        REQUIRE(&vr == &v);

        v2 = rand::laplace::fill(v, 2.0) + rand::laplace::fill(v, 3.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::laplace::fill(w, 2.0);
        w2 = rand::laplace::fill(w, 3.3) + rand::laplace::fill(w, 4.4) +
             rand::laplace::fill(w, 5.5);

        iexp::MatrixXd &wr = rand::laplace::fill(w, 99);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::laplace::fill(v1, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::laplace::pdf(v.array(), 2.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::laplace::pdf(m.array(), 2.0);

        // test compile
        v2 = rand::laplace::pdf(v.array(), 2.0) +
             rand::laplace::pdf(v.array(), 2.0);
        m2 = rand::laplace::pdf(m.array() + m.array(), 2.0);

#if 1 // #ifdef Ilaplace_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::laplace::pdf(vv.array(), 1.0);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("laplace_pdf.png");
#endif
    }
}

TEST_CASE("test_expow_fill")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::expow::fill(v, 1.0, 2.0);
        REQUIRE(&vr == &v);

        v2 = rand::expow::fill(v, 2.0, 100) + rand::expow::fill(v, 3.0, 4.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::expow::fill(w, 2.0, 3.3);
        w2 = rand::expow::fill(w, 3.3, 3.3) + rand::expow::fill(w, 4.4, 3.3) +
             rand::expow::fill(w, 5.5, 3.3);

        iexp::MatrixXd &wr = rand::expow::fill(w, 99, 3.3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::expow::fill(v1, 1.0, 2.5);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::expow::pdf(v.array(), 2.0, 3.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::expow::pdf(m.array(), 2.0, 3.0);

        // test compile
        v2 = rand::expow::pdf(v.array(), 2.0, 3.0) +
             rand::expow::pdf(v.array(), 2.0, 3.0);
        m2 = rand::expow::pdf(m.array() + m.array(), 2.0, 3.0);

#if 1 // #ifdef Iexpow_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -3, 3);
        VectorXd vv2 = rand::expow::pdf(vv.array(), 1.0, 2.5);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-3, 3, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("expow_pdf.png");
#endif
    }
}

TEST_CASE("test_cauchy_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::cauchy::fill(v, 1.0);
        REQUIRE(&vr == &v);

        v2 = rand::cauchy::fill(v, 2.0) + rand::cauchy::fill(v, 3.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::cauchy::fill(w, 2.0);
        w2 = rand::cauchy::fill(w, 3.3) + rand::cauchy::fill(w, 4.4) +
             rand::cauchy::fill(w, 5.5);

        iexp::MatrixXd &wr = rand::cauchy::fill(w, 99);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::cauchy::fill(v1, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::cauchy::pdf(v.array(), 2.0);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::cauchy::pdf(m.array(), 2.0);

        // test compile
        v2 = rand::cauchy::pdf(v.array(), 2.0) +
             rand::cauchy::pdf(v.array(), 2.0);
        m2 = rand::cauchy::pdf(m.array() + m.array(), 2.0);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::cauchy::pdf(vv.array(), 1.0);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("cauchy_pdf.png");
#endif
    }
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
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::landau::fill(v);
        REQUIRE(&vr == &v);

        v2 = rand::landau::fill(v) + rand::landau::fill(v);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::landau::fill(w);
        w2 = rand::landau::fill(w) + rand::landau::fill(w) +
             rand::landau::fill(w);

        iexp::MatrixXd &wr = rand::landau::fill(w);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::landau::fill(v1);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::landau::pdf(v.array());

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::landau::pdf(m.array());

        // test compile
        v2 = rand::landau::pdf(v.array()) + rand::landau::pdf(v.array());
        m2 = rand::landau::pdf(m.array() + m.array());

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::landau::pdf(vv.array());

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 0.2);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("laudau_pdf.png");
#endif
    }
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
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::flat::fill(v, 1, 2);
        REQUIRE(&vr == &v);

        v2 = rand::flat::fill(v, 2, 3) + rand::flat::fill(v, 2, 4);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::flat::fill(w, 2, 3);
        w2 = rand::flat::fill(w, 2, 3) + rand::flat::fill(w, 2, 3) +
             rand::flat::fill(w, 2, 3);

        iexp::MatrixXd &wr = rand::flat::fill(w, 2, 3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::flat::fill(v1, 1, 5);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::flat::pdf(v.array(), 0, 1);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::flat::pdf(m.array(), 0, 1);

        // test compile
        v2 =
            rand::flat::pdf(v.array(), 0, 1) + rand::flat::pdf(v.array(), 0, 1);
        m2 = rand::flat::pdf(m.array() + m.array(), 0, 1);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::flat::pdf(vv.array(), 0.5, 2.5);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("flat_pdf.png");
#endif
    }
}

TEST_CASE("test_gamma_rand")
{
    iexp::VectorXd v(10), v2(10);

    iexp::VectorXd &vr = rand::gamma::fill(v, 1, 2);
    REQUIRE(&vr == &v);

    v2 = rand::gamma::fill(v, 2, 3) + rand::gamma::fill(v, 2, 4);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::gamma::fill(w, 2, 3);
    w2 = rand::gamma::fill(w, 2, 3) + rand::gamma::fill(w, 2, 3) +
         rand::flat::fill(w, 2, 3);

    iexp::MatrixXd &wr = rand::gamma::fill(w, 2, 3);
    REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
    VectorXd v1(100);
    rand::gamma::fill(v1, 2, 1);

    mglData y(100);
    y.Link(v1.data(), v1.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 100, 0, 6);
    gr.Axis();
    gr.Plot(y, "+");
    gr.WriteFrame("gamma_rand.png");
#endif
    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::gamma::pdf(v.array(), 1, 1);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::gamma::pdf(m.array(), 1, 1);

        // test compile
        v2 = rand::gamma::pdf(v.array(), 1, 1) +
             rand::gamma::pdf(v.array(), 1, 1);
        m2 = rand::gamma::pdf(m.array() + m.array(), 1, 1);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::gamma::pdf(vv.array(), 1, 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("gamma_pdf.png");
#endif
    }
}

TEST_CASE("test_beta")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::beta::fill(v, 1, 2);
        REQUIRE(&vr == &v);

        v2 = rand::beta::fill(v, 2, 3) + rand::beta::fill(v, 2, 4);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::beta::fill(w, 2, 3);
        w2 = rand::beta::fill(w, 2, 3) + rand::beta::fill(w, 2, 3) +
             rand::flat_rand(w, 2, 3);

        iexp::MatrixXd &wr = rand::beta::fill(w, 2, 3);
        REQUIRE(&wr == &w);

#if 0 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::beta::fill(v1, 2, 2);

        mglData y(100);
        y.Link(v1.data(), v1.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 100, 0, 1);
        gr.Axis();
        gr.Plot(y, "+");
        gr.WriteFrame("beta_fill.png");
#endif
    }

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::beta::pdf(v.array(), 2, 3);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::beta::pdf(m.array(), 2, 3);

        // test compile
        v2 = rand::beta::pdf(v, 2, 3) + rand::beta::pdf(v, 2, 3);
        m2 = rand::beta::pdf(m.array() + m.array(), 2, 3);

#if 1 // #ifdef IEXP_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::beta::pdf(vv.array(), 2, 2);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 2);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("beta_pdf.png");
#endif
    }
}

TEST_CASE("test_lgnorm_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::lgnorm::fill(v, 1, 2);
        REQUIRE(&vr == &v);

        v2 = rand::lgnorm::fill(v, 2, 3) + rand::lgnorm::fill(v, 2, 4);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::lgnorm::fill(w, 2, 3);
        w2 = rand::lgnorm::fill(w, 2, 3) + rand::lgnorm::fill(w, 2, 3) +
             rand::flat_rand(w, 2, 3);

        iexp::MatrixXd &wr = rand::lgnorm::fill(w, 2, 3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::lgnorm::fill(v1, 0, 1);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::lgnorm::pdf(v.array(), 0, 1);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::lgnorm::pdf(m.array(), 0, 1);

        // test compile
        v2 = rand::lgnorm::pdf(v.array(), 0, 1) +
             rand::lgnorm::pdf(v.array(), 0, 1);
        m2 = rand::lgnorm::pdf(m.array() + m.array(), 0, 1);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::lgnorm::pdf(vv.array(), 0, 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("lgnorm_pdf.png");
#endif
    }
}

TEST_CASE("test_chisq_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::chisq::fill(v, 1.0);
        REQUIRE(&vr == &v);

        v2 = rand::chisq::fill(v, 2.0) + rand::chisq::fill(v, 3.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::chisq::fill(w, 2.0);
        w2 = rand::chisq::fill(w, 3.3) + rand::chisq::fill(w, 4.4) +
             rand::chisq::fill(w, 5.5);

        iexp::MatrixXd &wr = rand::chisq::fill(w, 99);
        REQUIRE(&wr == &w);

#if 1 // #ifdef Ichisq_MGL2
        VectorXd v1(100);
        rand::chisq::fill(v1, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::chisq::pdf(v.array(), 1);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::chisq::pdf(m.array(), 1);

        // test compile
        v2 = rand::chisq::pdf(v.array(), 1) + rand::chisq::pdf(v.array(), 1);
        m2 = rand::chisq::pdf(m.array() + m.array(), 1);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::chisq::pdf(vv.array(), 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("chisq_pdf.png");
#endif
    }
}

TEST_CASE("test_f_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::f::fill(v, 1, 2);
        REQUIRE(&vr == &v);

        v2 = rand::f::fill(v, 2, 3) + rand::f::fill(v, 2, 4);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::f::fill(w, 2, 3);
        w2 = rand::f::fill(w, 2, 3) + rand::f::fill(w, 2, 3) +
             rand::f::fill(w, 2, 3);

        iexp::MatrixXd &wr = rand::f::fill(w, 2, 3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::f::fill(v1, 1, 1);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::f::pdf(v.array(), 1, 2);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::f::pdf(m.array(), 1, 2);

        // test compile
        v2 = rand::f::pdf(v.array(), 1, 2) + rand::f::pdf(v.array(), 1, 2);
        m2 = rand::f::pdf(m.array() + m.array(), 1, 2);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::f::pdf(vv.array(), 3, 2);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("f_pdf.png");
#endif
    }
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
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::lgst::fill(v, 1.0);
        REQUIRE(&vr == &v);

        v2 = rand::lgst::fill(v, 2.0) + rand::lgst::fill(v, 3.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::lgst::fill(w, 2.0);
        w2 = rand::lgst::fill(w, 3.3) + rand::lgst::fill(w, 4.4) +
             rand::lgst::fill(w, 5.5);

        iexp::MatrixXd &wr = rand::lgst::fill(w, 99);
        REQUIRE(&wr == &w);

#if 1 // #ifdef It_MGL2
        VectorXd v1(100);
        rand::lgst::fill(v1, 1);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::lgst::pdf(v.array(), 2);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::lgst::pdf(m.array(), 2);

        // test compile
        v2 = rand::lgst::pdf(v.array(), 2) + rand::lgst::pdf(v.array(), 2);
        m2 = rand::lgst::pdf(m.array() + m.array(), 2);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::lgst::pdf(vv.array(), 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 0.5);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("lgst_pdf.png");
#endif
    }
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
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::gbl1::fill(v, 1.0, 2.0);
        REQUIRE(&vr == &v);

        v2 = rand::gbl1::fill(v, 2.0, 100) + rand::gbl1::fill(v, 3.0, 4.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::gbl1::fill(w, 2.0, 3.3);
        w2 = rand::gbl1::fill(w, 3.3, 3.3) + rand::gbl1::fill(w, 4.4, 3.3) +
             rand::gbl1::fill(w, 5.5, 3.3);

        iexp::MatrixXd &wr = rand::gbl1::fill(w, 99, 3.3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::gbl1::fill(v1, 1.0, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::gbl1::pdf(v.array(), 2, 3);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::gbl1::pdf(m.array(), 2, 3);

        // test compile
        v2 =
            rand::gbl1::pdf(v.array(), 2, 3) + rand::gbl1::pdf(v.array(), 2, 3);
        m2 = rand::gbl1::pdf(m.array() + m.array(), 2, 3);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::gbl1::pdf(vv.array(), 1, 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(-5, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("gbl1_pdf.png");
#endif
    }
}

TEST_CASE("test_gbl2_rand")
{
    {
        iexp::VectorXd v(10), v2(10);

        iexp::VectorXd &vr = rand::gbl2::fill(v, 1.0, 2.0);
        REQUIRE(&vr == &v);

        v2 = rand::gbl2::fill(v, 2.0, 100) + rand::gbl2::fill(v, 3.0, 4.0);

        iexp::MatrixXd w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::gbl2::fill(w, 2.0, 3.3);
        w2 = rand::gbl2::fill(w, 3.3, 3.3) + rand::gbl2::fill(w, 4.4, 3.3) +
             rand::gbl2::fill(w, 5.5, 3.3);

        iexp::MatrixXd &wr = rand::gbl2::fill(w, 99, 3.3);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXd v1(100);
        rand::gbl2::fill(v1, 1.0, 1.0);

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

    {
        VectorXd v = VectorXd::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::gbl2::pdf(v.array(), 2, 3);

        Matrix2Xd m = Matrix2Xd::Random(2, 10);
        Matrix2Xd m2 = rand::gbl2::pdf(m.array(), 2, 3);

        // test compile
        v2 =
            rand::gbl2::pdf(v.array(), 2, 3) + rand::gbl2::pdf(v.array(), 2, 3);
        m2 = rand::gbl2::pdf(m.array() + m.array(), 2, 3);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
        VectorXd vv2 = rand::gbl2::pdf(vv.array(), 1, 1);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 5, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("gbl2_pdf.png");
#endif
    }
}

TEST_CASE("test_drch_rand")
{
    double alpha[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.5};
    iexp::VectorXd v(12), v2(12);

    iexp::VectorXd &vr = rand::drch_rand(v, 6, alpha);
    REQUIRE(&vr == &v);

    v2 = rand::drch_rand(v, 6, alpha) + rand::drch_rand(v, 6, alpha);

    iexp::MatrixXd w(3, 4), w2(3, 4);
    w.fill(9.9999);
    w2 = rand::drch_rand(w, 6, alpha);
    w2 = rand::drch_rand(w, 6, alpha) + rand::drch_rand(w, 6, alpha) +
         rand::drch_rand(w, 6, alpha);

    iexp::MatrixXd &wr = rand::drch_rand(w, 6, alpha);
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
    {
        const double p[4] = {0.1, 0.2, 0.3, 0.4};
        iexp::VectorXi v(10), v2(10);

        iexp::VectorXi &vr = rand::discrete::fill(v, 4, p);
        REQUIRE(&vr == &v);

        v2 = rand::discrete::fill(v, 4, p) + rand::discrete::fill(v, 4, p);

        iexp::MatrixXi w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::discrete::fill(w, 4, p);
        w2 = rand::discrete::fill(w, 4, p) + rand::discrete::fill(w, 4, p) +
             rand::discrete::fill(w, 4, p);

        iexp::MatrixXi &wr = rand::discrete::fill(w, 4, p);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXi vi(100);
        rand::discrete::fill(vi, 4, p);
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

    {
        const double p[4] = {0.1, 0.2, 0.3, 0.4};

        VectorXi v = VectorXi::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::discrete::pdf(v.array(), 4, p);

        Matrix2Xi m = Matrix2Xi::Random(2, 10);
        Matrix2Xd m2 = rand::discrete::pdf(m.array(), 4, p);

        // test compile
        v2 = rand::discrete::pdf(v.array(), 4, p) +
             rand::discrete::pdf(v.array(), 4, p);
        m2 = rand::discrete::pdf(m.array() + m.array(), 4, p);
    }
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
    {
        iexp::VectorXi v(10), v2(10);

        iexp::VectorXi &vr = rand::bnom::fill(v, 0.5, 9);
        REQUIRE(&vr == &v);

        v2 = rand::bnom::fill(v, 0.5, 9) + rand::bnom::fill(v, 0.5, 9);

        iexp::MatrixXi w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::bnom::fill(w, 0.5, 9);
        w2 = rand::bnom::fill(w, 0.5, 9) + rand::bnom::fill(w, 0.5, 9) +
             rand::bnom::fill(w, 0.5, 9);

        iexp::MatrixXi &wr = rand::bnom::fill(w, 0.5, 9);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXi vi(100);
        rand::bnom::fill(vi, 0.5, 9);
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

    {
        VectorXi v = VectorXi::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::bnom::pdf(v.array(), 0.5, 2);

        Matrix2Xi m = Matrix2Xi::Random(2, 10);
        Matrix2Xd m2 = rand::bnom::pdf(m.array(), 0.5, 2);

        // test compile
        v2 = rand::bnom::pdf(v.array(), 0.5, 2) +
             rand::bnom::pdf(v.array(), 0.5, 2);
        m2 = rand::bnom::pdf(m.array() + m.array(), 0.5, 2);

#if 1 // #ifdef IEXP_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
        VectorXd vv2 = rand::bnom::pdf(vv.cast<unsigned int>().array(), 0.5, 9);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 10, 0, 0.5);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("bnom_pdf.png");
#endif
    }
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
    {
        iexp::VectorXi v(10), v2(10);

        iexp::VectorXi &vr = rand::geo::fill(v, 0.5);
        REQUIRE(&vr == &v);

        v2 = rand::geo::fill(v, 0.5) + rand::geo::fill(v, 0.5);

        iexp::MatrixXi w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::geo::fill(w, 0.5);
        w2 = rand::geo::fill(w, 0.5) + rand::geo::fill(w, 0.5) +
             rand::geo::fill(w, 0.5);

        iexp::MatrixXi &wr = rand::geo::fill(w, 0.5);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXi vi(100);
        rand::geo::fill(vi, 0.5);
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

    {
        VectorXi v = VectorXi::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::geo::pdf(v.array(), 0.5);

        Matrix2Xi m = Matrix2Xi::Random(2, 10);
        Matrix2Xd m2 = rand::geo::pdf(m.array(), 0.5);

        // test compile
        v2 = rand::geo::pdf(v.array(), 0.5) + rand::geo::pdf(v.array(), 0.5);
        m2 = rand::geo::pdf(m.array() + m.array(), 0.5);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
        VectorXd vv2 = rand::geo::pdf(vv.cast<unsigned int>().array(), 0.5);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 10, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("geometric_pdf.png");
#endif
    }
}

TEST_CASE("test_hgeo_rand")
{
    {
        iexp::VectorXi v(10), v2(10);

        iexp::VectorXi &vr = rand::hgeo::fill(v, 3, 10, 5);
        REQUIRE(&vr == &v);

        v2 = rand::hgeo::fill(v, 3, 10, 5) + rand::hgeo::fill(v, 3, 10, 5);

        iexp::MatrixXi w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::hgeo::fill(w, 3, 10, 5);
        w2 = rand::hgeo::fill(w, 3, 10, 5) + rand::hgeo::fill(w, 3, 10, 5) +
             rand::hgeo::fill(w, 3, 10, 5);

        iexp::MatrixXi &wr = rand::hgeo::fill(w, 3, 10, 5);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXi vi(100);
        rand::hgeo::fill(vi, 5, 20, 3);
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

    {
        VectorXi v = VectorXi::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::hgeo::pdf(v.array(), 1, 2, 3);

        Matrix2Xi m = Matrix2Xi::Random(2, 10);
        Matrix2Xd m2 = rand::hgeo::pdf(m.array(), 1, 2, 3);

        // test compile
        v2 = rand::hgeo::pdf(v.array(), 1, 2, 3) +
             rand::hgeo::pdf(v.array(), 1, 2, 3);
        m2 = rand::hgeo::pdf(m.array() + m.array(), 1, 2, 3);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
        VectorXd vv2 =
            rand::hgeo::pdf(vv.cast<unsigned int>().array(), 5, 20, 3);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 10, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("hgeo_pdf.png");
#endif
    }
}

TEST_CASE("test_log_rand")
{
    {
        iexp::VectorXi v(10), v2(10);

        iexp::VectorXi &vr = rand::log::fill(v, 0.7);
        REQUIRE(&vr == &v);

        v2 = rand::log::fill(v, 0.7) + rand::log::fill(v, 0.7);

        iexp::MatrixXi w(3, 4), w2(3, 4);
        w.fill(9.9999);
        w2 = rand::log::fill(w, 0.7);
        w2 = rand::log::fill(w, 0.7) + rand::log::fill(w, 0.7) +
             rand::log::fill(w, 0.7);

        iexp::MatrixXi &wr = rand::log::fill(w, 0.7);
        REQUIRE(&wr == &w);

#if 1 // #ifdef IEXP_MGL2
        VectorXi vi(100);
        rand::log::fill(vi, 0.7);
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

    {
        VectorXi v = VectorXi::LinSpaced(10, 0, 3);
        VectorXd v2 = rand::log::pdf(v.array(), 0.5);

        Matrix2Xi m = Matrix2Xi::Random(2, 10);
        Matrix2Xd m2 = rand::log::pdf(m.array(), 0.5);

        // test compile
        v2 = rand::log::pdf(v.array(), 0.5) + rand::log::pdf(v.array(), 0.5);
        m2 = rand::log::pdf(m.array() + m.array(), 0.5);

#if 1 // #ifdef Icauchy_MGL2
        VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
        VectorXd vv2 = rand::log::pdf(vv.cast<unsigned int>().array(), 0.7);

        mglData x(100), y(100);
        x.Link(vv.data(), vv.size());
        y.Link(vv2.data(), vv2.size());
        mglGraph gr;
        gr.SetOrigin(0, 0);
        gr.SetRanges(0, 10, 0, 1);
        gr.Axis();
        gr.Plot(x, y, "+");
        gr.WriteFrame("log_pdf.png");
#endif
    }
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

    iexp::RowVectorXi vr2 = rand::choose<true>(v, 9).matrix();
    REQUIRE(vr2.size() == 9);

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
