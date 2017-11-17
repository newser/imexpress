#include <iostream>

#include <Eigen/Core>
#include <gsl/gsl_const.h>

#include <common/common.h>

int main()
{
    std::cout << GSL_CONST_NUM_ZETTA << std::endl;

    Eigen::MatrixXd m(1, 2);
    std::cout << m << std::endl;

    Eigen::Array<int, 2, 2> a;
    // std::cout << Eigen::isnan(a) << Eigen::pos;

    return 0;
}
