#include <catch.hpp>
#include <iostream>
#include <rand/mul_gauss.h>
#include <randist/beta.h>
#include <randist/bi_gauss.h>
#include <randist/binomial.h>
#include <randist/cauchy.h>
#include <randist/chisq.h>
#include <randist/discrete.h>
#include <randist/exp.h>
#include <randist/expow.h>
#include <randist/f.h>
#include <randist/flat.h>
#include <randist/gamma.h>
#include <randist/gauss.h>
#include <randist/gauss_tail.h>
#include <randist/geometric.h>
#include <randist/gumbel1.h>
#include <randist/gumbel2.h>
#include <randist/hyper_geometric.h>
#include <randist/landau.h>
#include <randist/laplace.h>
#include <randist/lgnorm.h>
#include <randist/log.h>
#include <randist/logistic.h>
#include <randist/mul_gauss.h>
#include <randist/mul_nomial.h>
#include <randist/neg_binomial.h>
#include <randist/pareto.h>
#include <randist/pascal.h>
#include <randist/poisson.h>
#include <randist/rayleigh.h>
#include <randist/rayleigh_tail.h>
#include <randist/t.h>
#include <randist/weibull.h>
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

TEST_CASE("randist_exp")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::exp_pdf(v.array(), 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::exp_pdf(m.array(), 2.0);

    // test compile
    v2 = rdist::exp_pdf(v.array(), 2.0) + rdist::exp_pdf(v.array(), 2.0);
    m2 = rdist::exp_pdf(m.array() + m.array(), 2.0);

#if 0 // #ifdef IEXP_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::exp_pdf(vv.array(), 1.0);
    
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

TEST_CASE("randist_laplace")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::laplace_pdf(v.array(), 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::laplace_pdf(m.array(), 2.0);

    // test compile
    v2 =
        rdist::laplace_pdf(v.array(), 2.0) + rdist::laplace_pdf(v.array(), 2.0);
    m2 = rdist::laplace_pdf(m.array() + m.array(), 2.0);

#if 0 // #ifdef Ilaplace_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::laplace_pdf(vv.array(), 1.0);

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

TEST_CASE("randist_expow")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::expow_pdf(v.array(), 2.0, 3.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::expow_pdf(m.array(), 2.0, 3.0);

    // test compile
    v2 = rdist::expow_pdf(v.array(), 2.0, 3.0) +
         rdist::expow_pdf(v.array(), 2.0, 3.0);
    m2 = rdist::expow_pdf(m.array() + m.array(), 2.0, 3.0);

#if 0 // #ifdef Iexpow_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -3, 3);
    VectorXd vv2 = rdist::expow_pdf(vv.array(), 1.0, 2.5);
    
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

TEST_CASE("randist_cauchy")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::cauchy_pdf(v.array(), 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::cauchy_pdf(m.array(), 2.0);

    // test compile
    v2 = rdist::cauchy_pdf(v.array(), 2.0) + rdist::cauchy_pdf(v.array(), 2.0);
    m2 = rdist::cauchy_pdf(m.array() + m.array(), 2.0);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::cauchy_pdf(vv.array(), 1.0);
    
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

TEST_CASE("randist_rayleigh")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::rayl_pdf(v.array(), 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::rayl_pdf(m.array(), 2.0);

    // test compile
    v2 = rdist::rayl_pdf(v.array(), 2.0) + rdist::rayl_pdf(v.array(), 2.0);
    m2 = rdist::rayl_pdf(m.array() + m.array(), 2.0);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::rayl_pdf(vv.array(), 2.0);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-5, 5, 0, 1);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("rayl_pdf.png");
#endif
}

TEST_CASE("randist_rayleigh_tail")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::raylt_pdf(v.array(), 1.0, 2.0);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::raylt_pdf(m.array(), 1.0, 2.0);

    // test compile
    v2 = rdist::raylt_pdf(v.array(), 1.0, 2.0) +
         rdist::raylt_pdf(v.array(), 1.0, 2.0);
    m2 = rdist::raylt_pdf(m.array() + m.array(), 1.0, 2.0);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::raylt_pdf(vv.array(), 0.5, 2.0);

    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-5, 5, 0, 1);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("raylt_pdf.png");
#endif
}

TEST_CASE("randist_laudau")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::laudau_pdf(v.array());

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::laudau_pdf(m.array());

    // test compile
    v2 = rdist::laudau_pdf(v.array()) + rdist::laudau_pdf(v.array());
    m2 = rdist::laudau_pdf(m.array() + m.array());

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::laudau_pdf(vv.array());
    
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

TEST_CASE("randist_gamma")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::gamma_pdf(v.array(), 1, 1);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::gamma_pdf(m.array(), 1, 1);

    // test compile
    v2 = rdist::gamma_pdf(v.array(), 1, 1) + rdist::gamma_pdf(v.array(), 1, 1);
    m2 = rdist::gamma_pdf(m.array() + m.array(), 1, 1);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::gamma_pdf(vv.array(), 1, 1);
    
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

TEST_CASE("randist_flat")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::flat_pdf(v.array(), 0, 1);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::flat_pdf(m.array(), 0, 1);

    // test compile
    v2 = rdist::flat_pdf(v.array(), 0, 1) + rdist::flat_pdf(v.array(), 0, 1);
    m2 = rdist::flat_pdf(m.array() + m.array(), 0, 1);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::flat_pdf(vv.array(), 0.5, 2.5);
    
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

TEST_CASE("randist_lognormal")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::lgnorm_pdf(v.array(), 0, 1);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::lgnorm_pdf(m.array(), 0, 1);

    // test compile
    v2 =
        rdist::lgnorm_pdf(v.array(), 0, 1) + rdist::lgnorm_pdf(v.array(), 0, 1);
    m2 = rdist::lgnorm_pdf(m.array() + m.array(), 0, 1);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::lgnorm_pdf(vv.array(), 0, 1);
    
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

TEST_CASE("randist_chisq")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::chisq_pdf(v.array(), 1);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::chisq_pdf(m.array(), 1);

    // test compile
    v2 = rdist::chisq_pdf(v.array(), 1) + rdist::chisq_pdf(v.array(), 1);
    m2 = rdist::chisq_pdf(m.array() + m.array(), 1);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::chisq_pdf(vv.array(), 1);
    
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

TEST_CASE("randist_f")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::f_pdf(v.array(), 1, 2);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::f_pdf(m.array(), 1, 2);

    // test compile
    v2 = rdist::f_pdf(v.array(), 1, 2) + rdist::f_pdf(v.array(), 1, 2);
    m2 = rdist::f_pdf(m.array() + m.array(), 1, 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::f_pdf(vv.array(), 3, 2);
    
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

TEST_CASE("randist_t")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::t_pdf(v.array(), 2);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::t_pdf(m.array(), 2);

    // test compile
    v2 = rdist::t_pdf(v.array(), 2) + rdist::t_pdf(v.array(), 2);
    m2 = rdist::t_pdf(m.array() + m.array(), 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::t_pdf(vv.array(), 5);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(-5, 5, 0, 1);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("t_pdf.png");
#endif
}

TEST_CASE("randist_beta")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::beta_pdf(v.array(), 2, 3);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::beta_pdf(m.array(), 2, 3);

    // test compile
    v2 = rdist::beta_pdf(v.array(), 2, 3) + rdist::beta_pdf(v.array(), 2, 3);
    m2 = rdist::beta_pdf(m.array() + m.array(), 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::beta_pdf(vv.array(), 2, 2);
    
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

TEST_CASE("randist_logistic")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::lgst_pdf(v.array(), 2);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::lgst_pdf(m.array(), 2);

    // test compile
    v2 = rdist::lgst_pdf(v.array(), 2) + rdist::lgst_pdf(v.array(), 2);
    m2 = rdist::lgst_pdf(m.array() + m.array(), 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::lgst_pdf(vv.array(), 1);
    
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

TEST_CASE("randist_pareto")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::pareto_pdf(v.array(), 2, 3);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::pareto_pdf(m.array(), 2, 3);

    // test compile
    v2 =
        rdist::pareto_pdf(v.array(), 2, 3) + rdist::pareto_pdf(v.array(), 2, 3);
    m2 = rdist::pareto_pdf(m.array() + m.array(), 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::pareto_pdf(vv.array(), 3, 2);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 5, 0, 2);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("pareto_pdf.png");
#endif
}

TEST_CASE("randist_weibull")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::wbl_pdf(v.array(), 2, 3);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::wbl_pdf(m.array(), 2, 3);

    // test compile
    v2 = rdist::wbl_pdf(v.array(), 2, 3) + rdist::wbl_pdf(v.array(), 2, 3);
    m2 = rdist::wbl_pdf(m.array() + m.array(), 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::wbl_pdf(vv.array(), 1, 2);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 5, 0, 1);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("wbl_pdf.png");
#endif
}

TEST_CASE("randist_gumbel1")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::gbl1_pdf(v.array(), 2, 3);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::gbl1_pdf(m.array(), 2, 3);

    // test compile
    v2 = rdist::gbl1_pdf(v.array(), 2, 3) + rdist::gbl1_pdf(v.array(), 2, 3);
    m2 = rdist::gbl1_pdf(m.array() + m.array(), 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::gbl1_pdf(vv.array(), 1, 1);
    
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

TEST_CASE("randist_gumbel2")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::gbl2_pdf(v.array(), 2, 3);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::gbl2_pdf(m.array(), 2, 3);

    // test compile
    v2 = rdist::gbl2_pdf(v.array(), 2, 3) + rdist::gbl2_pdf(v.array(), 2, 3);
    m2 = rdist::gbl2_pdf(m.array() + m.array(), 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, -5, 5);
    VectorXd vv2 = rdist::gbl2_pdf(vv.array(), 1, 1);

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

TEST_CASE("randist_discrete")
{
    const double p[4] = {0.1, 0.2, 0.3, 0.4};

    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::discrete_pdf(v.array(), 4, p);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::discrete_pdf(m.array(), 4, p);

    // test compile
    v2 = rdist::discrete_pdf(v.array(), 4, p) +
         rdist::discrete_pdf(v.array(), 4, p);
    m2 = rdist::discrete_pdf(m.array() + m.array(), 4, p);
}

TEST_CASE("randist_poisson")
{
    VectorXd v = VectorXd::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::poiss_pdf(v.array(), 2);

    Matrix2Xd m = Matrix2Xd::Random(2, 10);
    Matrix2Xd m2 = rdist::poiss_pdf(m.array(), 2);

    // test compile
    v2 = rdist::poiss_pdf(v.array(), 2) + rdist::poiss_pdf(v.array(), 2);
    m2 = rdist::poiss_pdf(m.array() + m.array(), 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::poiss_pdf(vv.array(), 2.5);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 10, 0, 0.5);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("poiss_pdf.png");
#endif
}

TEST_CASE("randist_binomial")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::bnom_pdf(v.array(), 0.5, 2);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::bnom_pdf(m.array(), 0.5, 2);

    // test compile
    v2 =
        rdist::bnom_pdf(v.array(), 0.5, 2) + rdist::bnom_pdf(v.array(), 0.5, 2);
    m2 = rdist::bnom_pdf(m.array() + m.array(), 0.5, 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::bnom_pdf(vv.cast<unsigned int>().array(), 0.5, 9);
    
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

TEST_CASE("randist_multinomial")
{
    const double p[] = {0.1, 0.2, 0.3, 0.4};
    rdist::mnom mn(4, p);

    const unsigned int x[4] = {1, 2, 3, 4};
    double p1 = mn.pdf(x);
    double p2 = mn.lnpdf(x);
}

TEST_CASE("randist_neg_binomial")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::nbnom_pdf(v.array(), 0.5, 2);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::nbnom_pdf(m.array(), 0.5, 2);

    // test compile
    v2 = rdist::nbnom_pdf(v.array(), 0.5, 2) +
         rdist::bnom_pdf(v.array(), 0.5, 2);
    m2 = rdist::nbnom_pdf(m.array() + m.array(), 0.5, 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::nbnom_pdf(vv.cast<unsigned int>().array(), 0.5, 3.5);

    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 10, 0, 0.5);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("nbnom_pdf.png");
#endif
}

TEST_CASE("randist_pascal")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::pascal_pdf(v.array(), 0.5, 2);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::pascal_pdf(m.array(), 0.5, 2);

    // test compile
    v2 = rdist::pascal_pdf(v.array(), 0.5, 2) +
         rdist::pascal_pdf(v.array(), 0.5, 2);
    m2 = rdist::pascal_pdf(m.array() + m.array(), 0.5, 2);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::pascal_pdf(vv.cast<unsigned int>().array(), 0.5, 3);
    
    mglData x(100), y(100);
    x.Link(vv.data(), vv.size());
    y.Link(vv2.data(), vv2.size());
    mglGraph gr;
    gr.SetOrigin(0, 0);
    gr.SetRanges(0, 10, 0, 0.5);
    gr.Axis();
    gr.Plot(x, y, "+");
    gr.WriteFrame("pascal_pdf.png");
#endif
}

TEST_CASE("randist_geometric")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::geo_pdf(v.array(), 0.5);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::geo_pdf(m.array(), 0.5);

    // test compile
    v2 = rdist::geo_pdf(v.array(), 0.5) + rdist::geo_pdf(v.array(), 0.5);
    m2 = rdist::geo_pdf(m.array() + m.array(), 0.5);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::geo_pdf(vv.cast<unsigned int>().array(), 0.5);
    
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

TEST_CASE("randist_hyper_geometric")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::hgeo_pdf(v.array(), 1, 2, 3);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::hgeo_pdf(m.array(), 1, 2, 3);

    // test compile
    v2 = rdist::hgeo_pdf(v.array(), 1, 2, 3) + rdist::geo_pdf(v.array(), 0.5);
    m2 = rdist::hgeo_pdf(m.array() + m.array(), 1, 2, 3);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::hgeo_pdf(vv.cast<unsigned int>().array(), 5, 20, 3);

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

TEST_CASE("randist_log")
{
    VectorXi v = VectorXi::LinSpaced(10, 0, 3);
    VectorXd v2 = rdist::log_pdf(v.array(), 0.5);

    Matrix2Xi m = Matrix2Xi::Random(2, 10);
    Matrix2Xd m2 = rdist::log_pdf(m.array(), 0.5);

    // test compile
    v2 = rdist::log_pdf(v.array(), 0.5) + rdist::geo_pdf(v.array(), 0.5);
    m2 = rdist::log_pdf(m.array() + m.array(), 0.5);

#if 0 // #ifdef Icauchy_MGL2
    VectorXd vv = VectorXd::LinSpaced(100, 0, 10);
    VectorXd vv2 = rdist::log_pdf(vv.cast<unsigned int>().array(), 0.7);
    
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
