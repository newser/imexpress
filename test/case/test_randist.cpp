#include <catch.hpp>
#include <iostream>
#include <rand/beta.h>
#include <rand/binomial.h>
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
