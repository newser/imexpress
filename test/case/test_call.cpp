#include <Eigen/Dense>
#include <catch.hpp>
#include <gsl/gsl_sf_bessel.h>

#include <iostream>

using namespace Eigen;

int call_gsl()
{
    double x = 5.0;
    double y = gsl_sf_bessel_J0(x);
    printf("J0(%g) = %.18e\n", x, y);
    return 0;
}

TEST_CASE("call_gsl")
{
    REQUIRE(call_gsl() == 0);
}

int call_eigen()
{
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    std::cout << "Here is the matrix m:\n" << m << std::endl;
    VectorXd v(2);
    v(0) = 4;
    v(1) = v(0) - 1;
    std::cout << "Here is the vector v:\n" << v << std::endl;
    return 0;
}

TEST_CASE("call_eigen")
{
    REQUIRE(call_eigen() == 0);
}
