#include <catch.hpp>
#include <dae/linsol.h>
#include <dae/sunmat.h>
#include <dae/sunvec.h>
#include <iostream>
#include <test_util.h>

using namespace iexp;

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

class test_ls : public dae::linsol::solver<test_ls>
{
  public:
    test_ls()
        : dae::linsol::solver<test_ls>()
    {
    }

    void test_ec(int err_code)
    {
        check(err_code);
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
