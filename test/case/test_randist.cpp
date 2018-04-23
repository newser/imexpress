#include <catch.hpp>
#include <iostream>
#include <rand/beta.h>
#include <rand/binomial.h>
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
         rand::bnom::pdf(v.array(), 0.5, 2);
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
