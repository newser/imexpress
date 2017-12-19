#include <catch.hpp>
#include <iostream>
#include <randist/gauss.h>
#include <randist/gauss_tail.h>
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
