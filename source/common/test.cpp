#include <iostream>

#include <gsl/gsl_const.h>
#include <Eigen/Core>

int main()
{
    std::cout << GSL_CONST_NUM_ZETTA << std::endl;
    
    Eigen::MatrixXd m(1, 2);
    std::cout << m << std::endl;

    return 0;
}
