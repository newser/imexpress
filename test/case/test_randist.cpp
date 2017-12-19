#include <catch.hpp>
#include <iostream>
#include <rand/mul_gauss.h>
#include <randist/bi_gauss.h>
#include <randist/bi_gauss.h>
#include <randist/gauss.h>
#include <randist/gauss_tail.h>
#include <randist/mul_gauss.h>
#include <test_util.h>

using namespace iexp;
using namespace std;

TEST_CASE("randist_gauss")
{
    VectorXd v = VectorXd::LinSpaced(10, -5, 5);
    VectorXd v2 = rdist::gauss_pdf(v.array(), 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::gauss_pdf(m.array(), 2.0);

    // test compile
    v2 = rdist::gauss_pdf(v.array(), 2.0) + rdist::gauss_pdf(v.array(), 2.0);
    m2 = rdist::gauss_pdf(m.array() + m.array(), 2.0);

#if 0 // #ifdef IEXP_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::gauss_pdf(vv.array(), 2.0);

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

TEST_CASE("randist_gausst")
{
    VectorXd v = VectorXd::LinSpaced(10, -5, 5);
    VectorXd v2 = rdist::gausst_pdf(v.array(), 1.5, 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::gausst_pdf(m.array(), 1.5, 2.0);

    // test compile
    v2 = rdist::gausst_pdf(v.array(), 1.5, 2.0) +
         rdist::gausst_pdf(v.array(), 1.5, 2.0);
    m2 = rdist::gausst_pdf(m.array() + m.array(), 1.5, 2.0);

#if 0 // #ifdef IEXP_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 5);
    VectorXd vv2 = rdist::gausst_pdf(vv.array(), 1.5, 1.0);

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

TEST_CASE("randist_bgauss")
{
    VectorXd vv = VectorXd::LinSpaced(10, -2, 2);
    double z[100];

    rdist::bgauss bg(1, 1, 0.9);
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
    //gr.SetOrigin(0, 0);
    //gr.SetRanges(-2, 2, -2, 2, -2, 2);
    gr.Light(true);
    gr.Rotate(50, 60);
    gr.Box();
    gr.Surf(zz);
    gr.WriteFrame("bgauss_pdf.png");
#endif
}

TEST_CASE("randist_mgauss")
{
    double x[2] = {0};
    double mu[2] = {1, 2};
    double cov[4] = {4, 2, 2, 3};

    rdist::mgauss mg(2, mu, cov);
    double p = mg.pdf(x);
    REQUIRE(__D_EQ9(p, 0.028294217120391));

    p = mg.lnpdf(x);
    REQUIRE(__D_EQ9(p, log(0.028294217120391)));

    rand::mgauss_rng r(2, mu, cov);
    double sample[10 * 2], e_mu[2], e_cov[4];
    for (int i = 0; i < 10; ++i) {
        r.next(&sample[i * 2]);
    }
    rdist::mgauss::mean(10, 2, sample, e_mu);
    // cout << e_mu[0] << " " << e_mu[1] << endl;
    rdist::mgauss::cov(10, 2, sample, e_cov);
    // cout << e_cov[0] << " " << e_cov[1] << endl;
    // cout << e_cov[2] << " " << e_cov[3] << endl;

    MatrixXd m(2, 100);
    for (int i = 0; i < m.cols(); ++i) {
        double s[2];
        r.next(s);
        m(0, i) = s[0];
        m(1, i) = s[1];
    }
    Vector2d v = rdist::mgauss_mean(m.array());
    // cout << v << endl;

    Matrix<double, 100, 2, RowMajor> rm;
    for (int i = 0; i < rm.rows(); ++i) {
        double s[2];
        r.next(s);
        rm(i, 0) = s[0];
        rm(i, 1) = s[1];
    }
    Matrix2d rcov = rdist::mgauss_cov(rm.array());
    // cout << rcov << endl;
}
