#include <catch.hpp>
#include <dae/linsol.h>
#include <dae/sunmatrix.h>
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

    // no alloc
    m1 << 1, 2, 3, 4;
    dae::sunmat_dense smd(m1);
    SUNMatrix s = smd.sunmat();
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

    // row major
    Matrix<double, 2, 2, RowMajor> rm1;
    rm1 << 4, 5, 9, 1;
    dae::sunmat_dense smd2(rm1);
    SUNMatrix s2 = smd2.sunmat();
    REQUIRE(SM_ROWS_D(s2) == rm1.rows());
    REQUIRE(SM_COLUMNS_D(s2) == rm1.cols());
    REQUIRE(SM_LDATA_D(s2) == rm1.size());
    REQUIRE(SM_DATA_D(s2) != rm1.data());
    REQUIRE(SM_ELEMENT_D(s2, 0, 0) == rm1(0, 0));
    REQUIRE(SM_ELEMENT_D(s2, 0, 1) == rm1(0, 1));
    REQUIRE(SM_ELEMENT_D(s2, 1, 0) == rm1(1, 0));
    REQUIRE(SM_ELEMENT_D(s2, 1, 1) == rm1(1, 1));

    dae::sunmat_dense smd3(m1 + m1);
    SUNMatrix s3 = smd3.sunmat();
    REQUIRE(SM_ROWS_D(s3) == m1.rows());
    REQUIRE(SM_COLUMNS_D(s3) == m1.cols());
    REQUIRE(SM_LDATA_D(s3) == m1.size());
    REQUIRE(SM_DATA_D(s3) != m1.data());
    REQUIRE(SM_ELEMENT_D(s3, 0, 0) == 2 * m1(0, 0));
    REQUIRE(SM_ELEMENT_D(s3, 0, 1) == 2 * m1(0, 1));
    REQUIRE(SM_ELEMENT_D(s3, 1, 0) == 2 * m1(1, 0));
    REQUIRE(SM_ELEMENT_D(s3, 1, 1) == 2 * m1(1, 1));
}

TEST_CASE("test_linsol")
{
}
