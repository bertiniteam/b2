//This file is part of Bertini 2.
//
//products_of_linears_block_test.cpp is free software: you can redistribute it and/or
//modify it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//products_of_linears_block_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

#include <boost/test/unit_test.hpp>

#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/products_of_linears_block.hpp"

BOOST_AUTO_TEST_SUITE(products_of_linears_block_suite)

using namespace bertini;
using bertini::blocks::ProductsOfLinearsBlock;
using bertini::DefaultPrecision;

static_assert(bertini::blocks::is_block_v<ProductsOfLinearsBlock>,
              "ProductsOfLinearsBlock must satisfy the block contract");

// f0 = (2x + 3y + 1)(x - y + 4),  f1 = (x + 1)
// at (x,y) = (1,1):  f0 = 6*4 = 24,  f1 = 2
//   df0/dx = 2*4 + 6*1 = 14,  df0/dy = 3*4 + 6*(-1) = 6
//   df1/dx = 1,               df1/dy = 0
static ProductsOfLinearsBlock MakeTestBlock()
{
	bertini::Mat<complex_mp> f0(2, 3);
	f0 << complex_mp(2), complex_mp(3),  complex_mp(1),
	      complex_mp(1), complex_mp(-1), complex_mp(4);

	bertini::Mat<complex_mp> f1(1, 3);
	f1 << complex_mp(1), complex_mp(0), complex_mp(1);

	std::vector<bertini::Mat<complex_mp>> factors{f0, f1};
	return ProductsOfLinearsBlock(2, std::move(factors));
}

BOOST_AUTO_TEST_CASE(shape)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	BOOST_CHECK_EQUAL(block.NumFunctions(), 2u);
	BOOST_CHECK_EQUAL(block.NumVariables(), 2u);
	BOOST_CHECK(!block.DependsOnPathVariable());
}

BOOST_AUTO_TEST_CASE(eval_double)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Vec<complex_dbl> result(2);
	block.EvalInPlace<complex_dbl>(result, x, complex_dbl(0));

	BOOST_CHECK_CLOSE(result(0).real(), 24.0, 1e-11);
	BOOST_CHECK_CLOSE(result(1).real(), 2.0, 1e-11);
	BOOST_CHECK_SMALL(result(0).imag(), 1e-11);
	BOOST_CHECK_SMALL(result(1).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Mat<complex_dbl> J(2, 2);
	block.JacobianInPlace<complex_dbl>(J, x, complex_dbl(0));

	BOOST_CHECK_CLOSE(J(0, 0).real(), 14.0, 1e-11);
	BOOST_CHECK_CLOSE(J(0, 1).real(),  6.0, 1e-11);
	BOOST_CHECK_CLOSE(J(1, 0).real(),  1.0, 1e-11);
	BOOST_CHECK_SMALL(J(1, 1).real(), 1e-11);
}

BOOST_AUTO_TEST_CASE(eval_mpfr)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	block.Precision(30);

	bertini::Vec<complex_mp> x(2); x << complex_mp(1), complex_mp(1);
	bertini::Vec<complex_mp> result(2);
	block.EvalInPlace<complex_mp>(result, x, complex_mp(0));

	BOOST_CHECK(abs(result(0) - complex_mp(24)) < real_mp("1e-25"));
	BOOST_CHECK(abs(result(1) - complex_mp(2))  < real_mp("1e-25"));
}

BOOST_AUTO_TEST_CASE(jacobian_mpfr)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	block.Precision(30);

	bertini::Vec<complex_mp> x(2); x << complex_mp(1), complex_mp(1);
	bertini::Mat<complex_mp> J(2, 2);
	block.JacobianInPlace<complex_mp>(J, x, complex_mp(0));

	BOOST_CHECK(abs(J(0, 0) - complex_mp(14)) < real_mp("1e-25"));
	BOOST_CHECK(abs(J(0, 1) - complex_mp(6))  < real_mp("1e-25"));
	BOOST_CHECK(abs(J(1, 0) - complex_mp(1))  < real_mp("1e-25"));
	BOOST_CHECK(abs(J(1, 1))                    < real_mp("1e-25"));
}

BOOST_AUTO_TEST_SUITE_END()
