#include <catch.hpp>
#include <dae/dense_ode.h>
#include <dae/dense_ode.h>
#include <dae/func_ode.h>
#include <dae/function.h>
#include <dae/linsol.h>
#include <dae/ode.h>
#include <dae/sunmat.h>
#include <dae/sunvec.h>
#include <iostream>
#include <test_util.h>

using namespace iexp;
using namespace iexp::dae;

TEST_CASE("test_is_matrix_ar")
{
    Matrix2d m1, m2;
    Array2d a1, a2;

    REQUIRE(is_matrix<decltype(m1)>::value);
    REQUIRE(is_matrix<decltype(m1 + m2)>::value);
    REQUIRE(!is_matrix<decltype(a1)>::value);
    REQUIRE(!is_matrix<decltype(a1 + a2)>::value);

    REQUIRE(!is_array<decltype(m1)>::value);
    REQUIRE(!is_array<decltype(m1 + m2)>::value);
    REQUIRE(is_array<decltype(a1)>::value);
    REQUIRE(is_array<decltype(a1 + a2)>::value);
}

TEST_CASE("test_sunmat")
{
    Matrix2d m1;

    // no alloc
    m1 << 1, 2, 3, 4;
    dae::sunmat_dense smd(m1, false);
    SUNMatrix s = smd.sunmatrix();
    REQUIRE(SM_ROWS_D(s) == m1.rows());
    REQUIRE(SM_COLUMNS_D(s) == m1.cols());
    REQUIRE(SM_LDATA_D(s) == m1.size());
    REQUIRE(SM_DATA_D(s) == m1.data());
    REQUIRE(SM_ELEMENT_D(s, 0, 0) == m1(0, 0));
    REQUIRE(SM_ELEMENT_D(s, 0, 1) == m1(0, 1));
    REQUIRE(SM_ELEMENT_D(s, 1, 0) == m1(1, 0));
    REQUIRE(SM_ELEMENT_D(s, 1, 1) == m1(1, 1));
    SUNMatZero(s);
    REQUIRE(m1(0, 0) == 0);
    REQUIRE(m1(1, 1) == 0);

    // copy
    dae::sunmat_dense smd4(m1, true);
    SUNMatrix s4 = smd4.sunmatrix();
    REQUIRE(SM_ROWS_D(s4) == m1.rows());
    REQUIRE(SM_COLUMNS_D(s4) == m1.cols());
    REQUIRE(SM_LDATA_D(s4) == m1.size());
    REQUIRE(SM_DATA_D(s4) != m1.data());
    REQUIRE(SM_ELEMENT_D(s4, 0, 0) == m1(0, 0));
    REQUIRE(SM_ELEMENT_D(s4, 0, 1) == m1(0, 1));
    REQUIRE(SM_ELEMENT_D(s4, 1, 0) == m1(1, 0));
    REQUIRE(SM_ELEMENT_D(s4, 1, 1) == m1(1, 1));

    // row major
    Matrix<double, 2, 2, RowMajor> rm1;
    rm1 << 4, 5, 9, 1;
    dae::sunmat_dense smd2(rm1);
    SUNMatrix s2 = smd2.sunmatrix();
    REQUIRE(SM_ROWS_D(s2) == rm1.rows());
    REQUIRE(SM_COLUMNS_D(s2) == rm1.cols());
    REQUIRE(SM_LDATA_D(s2) == rm1.size());
    REQUIRE(SM_DATA_D(s2) != rm1.data());
    REQUIRE(SM_ELEMENT_D(s2, 0, 0) == rm1(0, 0));
    REQUIRE(SM_ELEMENT_D(s2, 0, 1) == rm1(0, 1));
    REQUIRE(SM_ELEMENT_D(s2, 1, 0) == rm1(1, 0));
    REQUIRE(SM_ELEMENT_D(s2, 1, 1) == rm1(1, 1));

    // expression
    dae::sunmat_dense smd3(m1 + m1);
    SUNMatrix s3 = smd3.sunmatrix();
    REQUIRE(SM_ROWS_D(s3) == m1.rows());
    REQUIRE(SM_COLUMNS_D(s3) == m1.cols());
    REQUIRE(SM_LDATA_D(s3) == m1.size());
    REQUIRE(SM_DATA_D(s3) != m1.data());
    REQUIRE(SM_ELEMENT_D(s3, 0, 0) == 2 * m1(0, 0));
    REQUIRE(SM_ELEMENT_D(s3, 0, 1) == 2 * m1(0, 1));
    REQUIRE(SM_ELEMENT_D(s3, 1, 0) == 2 * m1(1, 0));
    REQUIRE(SM_ELEMENT_D(s3, 1, 1) == 2 * m1(1, 1));
}

TEST_CASE("test_sunvec")
{
    Vector4d v1;
    v1 << 1, 2, 3, 4;

    dae::sunvec_serial sv1(v1, false);
    N_Vector nv1 = sv1.n_vector();
    REQUIRE(NV_LENGTH_S(nv1) == v1.size());
    REQUIRE(NV_DATA_S(nv1) == v1.data());
    REQUIRE(NV_OWN_DATA_S(nv1) == SUNFALSE);
    REQUIRE(NV_Ith_S(nv1, 0) == v1(0));
    REQUIRE(NV_Ith_S(nv1, 3) == v1(3));

    dae::sunvec_serial sv3(v1, true);
    N_Vector nv3 = sv3.n_vector();
    REQUIRE(NV_LENGTH_S(nv3) == v1.size());
    REQUIRE(NV_DATA_S(nv3) != v1.data());
    REQUIRE(NV_OWN_DATA_S(nv3) == SUNTRUE);
    REQUIRE(NV_Ith_S(nv3, 0) == v1(0));
    REQUIRE(NV_Ith_S(nv3, 3) == v1(3));

    dae::sunvec_serial sv2(v1 + v1);
    N_Vector nv2 = sv2.n_vector();
    REQUIRE(NV_LENGTH_S(nv2) == v1.size());
    REQUIRE(NV_DATA_S(nv2) != v1.data());
    REQUIRE(NV_OWN_DATA_S(nv2) == SUNTRUE);
    REQUIRE(NV_Ith_S(nv2, 0) == 2 * v1(0));
    REQUIRE(NV_Ith_S(nv2, 3) == 2 * v1(3));
}

class test_ls : public dae::linsol<test_ls>
{
  public:
    test_ls()
        : dae::linsol<test_ls>()
    {
    }

    void test_ec(int err_code)
    {
        ls_check(err_code);
    }
};

TEST_CASE("test_linsol_solver")
{
    test_ls ls;

    ls.test_ec(0);
    ls.test_ec(1);

    except_begin()
    {
        ls.test_ec(-1);
    }
    except_str("the memory argument to the function is NULL");

    except_begin()
    {
        ls.test_ec(-9);
    }
    except_str("a singular R matrix was encountered in a QR factorization");
}

class test_dy
{
  public:
    int compute_dy(double t, Map<const VectorXd> &y, Map<VectorXd> &dy)
    {
        for (int i = 0; i < dy.size(); ++i) {
            dy[i] = y[i] + t;
        }
        return 1;
    }
};

TEST_CASE("test_dy")
{
    N_Vector y = N_VNew_Serial(3);
    NV_Ith_S(y, 0) = 1;
    NV_Ith_S(y, 1) = 2;
    NV_Ith_S(y, 2) = 3;
    N_Vector dy = N_VNew_Serial(3);

    test_dy td;
    int ret = dy_func<test_dy>::s_dy(9, y, dy, &td);
    REQUIRE(ret == 1);
    REQUIRE(NV_Ith_S(dy, 0) == 10);
    REQUIRE(NV_Ith_S(dy, 1) == 11);
    REQUIRE(NV_Ith_S(dy, 2) == 12);
}

class test_jac
{
  public:
    int compute_jac(double t,
                    Map<const VectorXd> &y,
                    Map<const VectorXd> &fy,
                    Map<MatrixXd> &jac)
    {
        REQUIRE(y.size() == 3);
        REQUIRE(fy.size() == 3);
        REQUIRE(jac.rows() == 3);
        REQUIRE(jac.cols() == 3);
        for (int i = 0; i < jac.rows(); ++i) {
            for (int j = 0; j < jac.cols(); ++j) {
                jac(i, j) = y(i) * fy(j);
            }
        }
        return 2;
    }
};

TEST_CASE("test_jac")
{
    N_Vector y = N_VNew_Serial(3);
    NV_Ith_S(y, 0) = 1;
    NV_Ith_S(y, 1) = 2;
    NV_Ith_S(y, 2) = 3;

    N_Vector fy = N_VNew_Serial(3);
    NV_Ith_S(fy, 0) = 4;
    NV_Ith_S(fy, 1) = 5;
    NV_Ith_S(fy, 2) = 6;

    SUNMatrix j = SUNDenseMatrix(3, 3);

    test_jac tj;
    jac_func<test_jac>::s_jac(9, y, fy, j, &tj, nullptr, nullptr, nullptr);
    REQUIRE(SM_ELEMENT_D(j, 0, 0) == 4);
    REQUIRE(SM_ELEMENT_D(j, 0, 2) == 6);
    REQUIRE(SM_ELEMENT_D(j, 2, 0) == 12);
    REQUIRE(SM_ELEMENT_D(j, 2, 2) == 18);
}

double sol_val(double t)
{
    return -(0.5 + 2 * t) * exp(-6 * t);
}

TEST_CASE("test_func_ode")
{
    /*
     ivp:
     u''/2 + 6u' + 18u = 0
     u(0) = -1/2
     u'(0) = 1

     solution
     u = -e^(-6t)/2 - 2te^(-6t)

     converted:
     y0 = u
     y1 = u'
     */

    for (int i = 0; i < 2; ++i) {
        VectorXd y(2);
        y << -0.5, 1;

        multistep ms = (i == 0 ? multistep::ADAMS : multistep::BDF);

        func_ode dfo(ms,
                     [](double t,
                        Map<const VectorXd> &y,
                        Map<VectorXd> &dy) -> int {
                         dy[0] = y[1];
                         dy[1] = -12 * y[1] - 36 * y[0];
                         return CV_SUCCESS;
                     },
                     0,
                     y);

        dfo.tolerance(0, 1e-6);

        for (double t = 0.1; t < 1; t += 0.1) {
            double tv = t;
            dfo.go(tv, y);

            double ans = sol_val(tv);
            REQUIRE(__D_EQ6(y[0], ans));
        }
    }
}

TEST_CASE("test_dense_ode")
{
    /*
     ivp:
     u''/2 + 6u' + 18u = 0
     u(0) = -1/2
     u'(0) = 1

     solution
     u = -e^(-6t)/2 - 2te^(-6t)

     converted:
     y0 = u
     y1 = u'
     */

    for (int i = 0; i < 2; ++i) {
        VectorXd y(2);
        y << -0.5, 1;

        multistep ms = (i == 0 ? multistep::ADAMS : multistep::BDF);

        dense_ode dfo(ms,
                      [](double t,
                         Map<const VectorXd> &y,
                         Map<VectorXd> &dy) -> int {
                          dy[0] = y[1];
                          dy[1] = -12 * y[1] - 36 * y[0];
                          return CV_SUCCESS;
                      },
                      0,
                      y);

        dfo.tolerance(0, 1e-6);

        for (double t = 0.1; t < 1; t += 0.1) {
            double tv = t;
            dfo.go(tv, y);

            double ans = sol_val(tv);
            REQUIRE(__D_EQ6(y[0], ans));
        }
    }

    for (int i = 0; i < 2; ++i) {
        VectorXd y(2);
        y << -0.5, 1;

        multistep ms = (i == 0 ? multistep::ADAMS : multistep::BDF);

        dense_ode dfo(ms,
                      [](double t,
                         Map<const VectorXd> &y,
                         Map<VectorXd> &dy) -> int {
                          dy[0] = y[1];
                          dy[1] = -12 * y[1] - 36 * y[0];
                          return CV_SUCCESS;
                      },
                      0,
                      y);

        dfo.tolerance(0, 1e-6);

        dfo.jac([](double t,
                   Map<const VectorXd> &y,
                   Map<const VectorXd> &fy,
                   Map<MatrixXd> &jac) -> int {
            jac(0, 0) = 0;
            jac(0, 1) = 1;
            jac(1, 1) = -36;
            jac(1, 1) = -12;

            return CV_SUCCESS;
        });

        for (double t = 0.1; t < 1; t += 0.1) {
            double tv = t;
            dfo.go(tv, y);

            double ans = sol_val(tv);
            REQUIRE(__D_EQ6(y[0], ans));
        }
    }
}
