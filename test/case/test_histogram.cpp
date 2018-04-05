#include <../test/test_util.h>
#include <catch.hpp>
#include <histogram/hist.h>
#include <histogram/hist2.h>
#include <histogram/hist2pdf.h>
#include <histogram/histpdf.h>

using namespace iexp;

void test_empty(hist &h)
{
    REQUIRE(h[0] == 0);
    REQUIRE(h[3] == 0);
    double low, up;
    h.range(0, low, up);
    REQUIRE((low == 1.0 && up == 2.0));
    h.range(3, low, up);
    REQUIRE((low == 4.0 && up == 5.0));
    REQUIRE(h.max() == 5.0);
    REQUIRE(h.min() == 1.0);
    REQUIRE(h.size() == 4);
    REQUIRE(h.find(1.5) == 0);
    REQUIRE(h.find(3.5) == 2);
    REQUIRE(h.find(4.5) == 3);
    REQUIRE(h.max_val() == 0);
    REQUIRE(h.min_val() == 0);
    REQUIRE(h.mean() == 0);
    REQUIRE(h.std() == 0);
    REQUIRE(h.sum() == 0);
}

TEST_CASE("hist_empty")
{
    double r[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    hist h(r, sizeof(r) / sizeof(r[0]));

    test_empty(h);

    hist h2(h);
    test_empty(h2);

    hist h3(std::move(h2));
    test_empty(h3);

    hist h4 = h3;
    test_empty(h4);

    hist h5 = std::move(h4);
    test_empty(h5);
    REQUIRE(h5 == h);

    hist h6(4, 1.0, 5.0);
    test_empty(h6);

    h += h5;
    test_empty(h);
    h += 0;
    test_empty(h);
    h -= h5;
    test_empty(h);
    h -= 0;
    test_empty(h);
    h *= h5;
    test_empty(h);
    h *= 6.0;
    test_empty(h);
    h /= 6.0;
    test_empty(h);
}

void test_val(hist &h)
{
    REQUIRE(h[0] == 2);
    REQUIRE(h[3] == 1);
    double low, up;
    h.range(0, low, up);
    REQUIRE((low == 0.0 && up == 1.0));
    h.range(3, low, up);
    REQUIRE((low == 3.0 && up == 4.0));
    REQUIRE(h.max() == 10.0);
    REQUIRE(h.min() == 0.0);
    REQUIRE(h.size() == 10);
    REQUIRE(h.find(1.5) == 1);
    REQUIRE(h.find(3.5) == 3);
    REQUIRE(h.find(9.5) == 9);
    REQUIRE(h.max_val() == 2.0);
    REQUIRE(h.max_idx() == 0);
    REQUIRE(h.min_val() == -1.0);
    REQUIRE(h.min_idx() == 9);
    REQUIRE(__D_EQ7(h.mean(), 2.6428571428571428));
    REQUIRE(__D_EQ7(h.std(), 1.8070158058105026));
    REQUIRE(__D_EQ7(h.sum(), 6));
}

TEST_CASE("hist_val")
{
    //    double r[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    //    hist h(r, sizeof(r) / sizeof(r[0]));
    hist h({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});

    h << 0.5 << 1.5 << 2.5 << 3.5 << 4.5 << 5.5;
    h << 0.5;
    h.add(9.9, -1);
    test_val(h);

    hist h2(h);
    test_val(h2);

    hist h3(std::move(h2));
    test_val(h3);

    hist h4 = h3;
    test_val(h4);

    hist h5 = std::move(h4);
    test_val(h5);
    REQUIRE(h5 == h);

    hist h6(10, 0.0, 10.0);
    h6 << 0.5 << 1.5 << 2.5 << 3.5 << 4.5 << 5.5;
    h6 << 0.5;
    h6.add(9.9, -1);
    test_val(h6);

    h5.reset();
    h += h5;
    test_val(h);
    h += 0;
    test_val(h);
    h -= h5;
    test_val(h);
    h -= 0;
    test_val(h);

    h5 << 0.1 << 1.1 << 2.1 << 3.1 << 4.1 << 5.1 << 6.1 << 7.1 << 8.1;
    h5.add(9.1, 1.0);
    h *= h5;
    test_val(h);
    h *= 1.0;
    test_val(h);
    h /= h5;
    test_val(h);
    h /= 1.0;
    test_val(h);

    h.reset();
    histpdf hp(h);
    hp.next();
    hp.next();

    VectorXd v(10);
    hp.next(v);
    // hp.next(v + v);
}

TEST_CASE("hist_vec")
{
    // matrix
    {
        Matrix<double, 1, 11> m;
        m << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
        hist h(m);

        VectorXd x(6);
        x << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
        h << x;

        h << 0.5;
        h.add(9.9, -1);
        test_val(h);

        size_t i = h.find(100);
        REQUIRE(i == -1);

        hist h2(m + m);
    }

    // array
    {
        Array<double, 11, 1> a;
        a << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
        hist h(a);

        VectorXd x(6), w(6);
        x << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
        w.setOnes();
        h.add(x, w);

        h << 0.5;
        h.add(9.9, -1);
        test_val(h);

        size_t i = h.find(100);
        REQUIRE(i == -1);

        hist h2(a + a);
    }
}

void test_empty2(hist2 &h)
{
    REQUIRE(h.get(0, 0) == 0);
    REQUIRE(h.get(8, 2) == 0);
    REQUIRE(h.get(9, 4) == 0);

    double low, up;
    h.xrange(0, low, up);
    REQUIRE((low == 0.0 && up == 1.0));
    h.xrange(9, low, up);
    REQUIRE((low == 9.0 && up == 10.0));
    h.yrange(0, low, up);
    REQUIRE((low == 90.0 && up == 91.0));
    h.yrange(4, low, up);
    REQUIRE((low == 94.0 && up == 95.0));

    REQUIRE(h.xmax() == 10.0);
    REQUIRE(h.xmin() == 0.0);
    REQUIRE(h.xsize() == 10);
    REQUIRE(h.ymax() == 95.0);
    REQUIRE(h.ymin() == 90.0);
    REQUIRE(h.ysize() == 5);

    size_t i, j;
    h.find(0.5, 90.1, i, j);
    REQUIRE((i == 0 && j == 0));
    h.find(9.5, 94.1, i, j);
    REQUIRE((i == 9 && j == 4));
    h.find(4.5, 93.1, i, j);
    REQUIRE((i == 4 && j == 3));

    REQUIRE(h.max_val() == 0);
    REQUIRE(h.min_val() == 0);
    REQUIRE(h.xmean() == 0);
    REQUIRE(h.ymean() == 0);
    REQUIRE(h.xstd() == 0);
    REQUIRE(h.ystd() == 0);
    REQUIRE(h.sum() == 0);
    REQUIRE(h.cov() == 0);
}

TEST_CASE("hist2_empty")
{
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    double y[] = {90.0, 91.0, 92.0, 93.0, 94.0, 95.0};
    hist2 h(x, sizeof(x) / sizeof(x[0]), y, sizeof(y) / sizeof(y[0]));

    test_empty2(h);

    hist2 h2(h);
    test_empty2(h2);

    hist2 h3(std::move(h2));
    test_empty2(h3);

    hist2 h4 = h3;
    test_empty2(h4);

    hist2 h5 = std::move(h4);
    test_empty2(h5);
    REQUIRE(h5 == h);

    hist2 h6(10, 0.0, 10.0, 5, 90.0, 95.0);
    test_empty2(h6);

    h += h5;
    test_empty2(h);
    h += 0;
    test_empty2(h);
    h -= h5;
    test_empty2(h);
    h -= 0;
    test_empty2(h);
    h *= h5;
    test_empty2(h);
    h *= 6.0;
    test_empty2(h);
    h /= 6.0;
    test_empty2(h);
}

void test_val2(hist2 &h)
{
    REQUIRE(h.get(0, 0) == 1);
    REQUIRE(h.get(0, 1) == 0);
    REQUIRE(h.get(9, 4) == 1);
    REQUIRE(h.get(9, 3) == 0);

    double low, up;
    h.xrange(0, low, up);
    REQUIRE((low == 0.0 && up == 1.0));
    h.xrange(9, low, up);
    REQUIRE((low == 9.0 && up == 10.0));
    h.yrange(0, low, up);
    REQUIRE((low == 90.0 && up == 91.0));
    h.yrange(4, low, up);
    REQUIRE((low == 94.0 && up == 95.0));

    REQUIRE(h.xmax() == 10.0);
    REQUIRE(h.xmin() == 0.0);
    REQUIRE(h.xsize() == 10);
    REQUIRE(h.ymax() == 95.0);
    REQUIRE(h.ymin() == 90.0);
    REQUIRE(h.ysize() == 5);

    size_t i, j;
    h.find(0.5, 90.1, i, j);
    REQUIRE((i == 0 && j == 0));
    h.find(9.5, 94.1, i, j);
    REQUIRE((i == 9 && j == 4));
    h.find(4.5, 93.1, i, j);
    REQUIRE((i == 4 && j == 3));

    REQUIRE(h.max_val() == 3.2);
    h.max_idx(i, j);
    REQUIRE((i == 5 && j == 2));

    REQUIRE(h.min_val() == -2.0);
    h.min_idx(i, j);
    REQUIRE((i == 4 && j == 2));

    REQUIRE(__D_EQ7(h.xmean(), 5.1428571428571432));
    REQUIRE(__D_EQ7(h.xstd(), 2.7152254012497448));
    REQUIRE(__D_EQ7(h.ymean(), 92.5));
    REQUIRE(__D_EQ7(h.ystd(), 1.3363062095621219));
    REQUIRE(__D_EQ7(h.sum(), 9.1999999999999993));
    REQUIRE(__D_EQ7(h.cov(), 3.5714285714285712));
}

TEST_CASE("hist2_val")
{
    //    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    //    double y[] = {90.0, 91.0, 92.0, 93.0, 94.0, 95.0};
    //    hist2 h(x, sizeof(x) / sizeof(x[0]), y, sizeof(y) / sizeof(y[0]));
    hist2 h({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0},
            {90.0, 91.0, 92.0, 93.0, 94.0, 95.0});

    h.add(0.1, 90.1);
    h.add(1.1, 90.9);
    h.add(2.1, 91.1);
    h.add(3.1, 91.9);
    h.add(4.1, 92.1);
    h.add(4.1, 92.1, -3.0);
    h.add(5.1, 92.9);
    h.add(5.1, 92.9, 2.2);
    h.add(6.1, 93.1);
    h.add(7.1, 93.9);
    h.add(8.1, 94.1);
    h.add(9.1, 94.9);
    test_val2(h);

    hist2 h2(h);
    test_val2(h2);

    hist2 h3(std::move(h2));
    test_val2(h3);

    hist2 h4 = h3;
    test_val2(h4);

    hist2 h5 = std::move(h4);
    test_val2(h5);
    REQUIRE(h5 == h);

    hist2 h6(10, 0.0, 10.0, 5, 90.0, 95.0);
    h6.add(0.1, 90.1);
    h6.add(1.1, 90.9);
    h6.add(2.1, 91.1);
    h6.add(3.1, 91.9);
    h6.add(4.1, 92.1);
    h6.add(4.1, 92.1, -3.0);
    h6.add(5.1, 92.9);
    h6.add(5.1, 92.9, 2.2);
    h6.add(6.1, 93.1);
    h6.add(7.1, 93.9);
    h6.add(8.1, 94.1);
    h6.add(9.1, 94.9);
    test_val2(h6);

    h5.reset();
    h += h5;
    test_val2(h);
    h += 0;
    test_val2(h);
    h -= h5;
    test_val2(h);
    h -= 0;
    test_val2(h);

    int i, j;
    for (i = 0; i < 10; ++i) {
        for (j = 0; j < 5; ++j) {
            h5.add(i + 0.2, j + 90.2);
        }
    }
    h *= h5;
    test_val2(h);
    h *= 1.0;
    test_val2(h);
    h /= h5;
    test_val2(h);
    h /= 1.0;
    test_val2(h);

    h.reset();
    hist2pdf hp(h);
    double x1, y1;
    hp.next(x1, y1);
    hp.next(x1, y1);

    Matrix<double, 2, 3> mx;
    Array<double, 2, 3> ay;
    hp.next(mx, ay);
}

TEST_CASE("hist2_vec")
{
    // matrix
    {
        Matrix<double, 1, 11> m;
        m << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
        Array<double, 6, 1> my;
        my << 90.0, 91.0, 92.0, 93.0, 94.0, 95.0;
        hist2 h(m, my);

        VectorXd x(10), y(10);
        x << 0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1;
        y << 90.1, 90.9, 91.1, 91.9, 92.1, 92.9, 93.1, 93.9, 94.1, 94.9;
        h.add(x, y);
        h.add(4.1, 92.1, -3.0);
        h.add(5.1, 92.9, 2.2);
        test_val2(h);

        size_t i, j;
        REQUIRE(!h.find(9, 99, i, j));

        hist h2(m + m);
    }

    // array
    {
        Array<double, 11, 1> a;
        a << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0;
        Matrix<double, 6, 1> ay;
        ay << 90.0, 91.0, 92.0, 93.0, 94.0, 95.0;
        hist2 h(a, ay);

        VectorXd x(10);
        ArrayXd y(10);
        x << 0.1, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1, 9.1;
        y << 90.1, 90.9, 91.1, 91.9, 92.1, 92.9, 93.1, 93.9, 94.1, 94.9;
        h.add(x, y);
        h.add(4.1, 92.1, -3.0);
        h.add(5.1, 92.9, 2.2);
        test_val2(h);

        hist h2(a + a);
    }
}
