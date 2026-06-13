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
	Mat<mpfr_complex> f0(2, 3);
	f0 << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1),
	      mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);

	Mat<mpfr_complex> f1(1, 3);
	f1 << mpfr_complex(1), mpfr_complex(0), mpfr_complex(1);

	std::vector<Mat<mpfr_complex>> factors{f0, f1};
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

	Vec<dbl> x(2); x << dbl(1), dbl(1);
	Vec<dbl> result(2);
	block.EvalInPlace<dbl>(result, x, dbl(0));

	BOOST_CHECK_CLOSE(result(0).real(), 24.0, 1e-11);
	BOOST_CHECK_CLOSE(result(1).real(), 2.0, 1e-11);
	BOOST_CHECK_SMALL(result(0).imag(), 1e-11);
	BOOST_CHECK_SMALL(result(1).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();

	Vec<dbl> x(2); x << dbl(1), dbl(1);
	Mat<dbl> J(2, 2);
	block.JacobianInPlace<dbl>(J, x, dbl(0));

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

	Vec<mpfr_complex> x(2); x << mpfr_complex(1), mpfr_complex(1);
	Vec<mpfr_complex> result(2);
	block.EvalInPlace<mpfr_complex>(result, x, mpfr_complex(0));

	BOOST_CHECK(abs(result(0) - mpfr_complex(24)) < mpfr_float("1e-25"));
	BOOST_CHECK(abs(result(1) - mpfr_complex(2))  < mpfr_float("1e-25"));
}

BOOST_AUTO_TEST_CASE(jacobian_mpfr)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	block.Precision(30);

	Vec<mpfr_complex> x(2); x << mpfr_complex(1), mpfr_complex(1);
	Mat<mpfr_complex> J(2, 2);
	block.JacobianInPlace<mpfr_complex>(J, x, mpfr_complex(0));

	BOOST_CHECK(abs(J(0, 0) - mpfr_complex(14)) < mpfr_float("1e-25"));
	BOOST_CHECK(abs(J(0, 1) - mpfr_complex(6))  < mpfr_float("1e-25"));
	BOOST_CHECK(abs(J(1, 0) - mpfr_complex(1))  < mpfr_float("1e-25"));
	BOOST_CHECK(abs(J(1, 1))                    < mpfr_float("1e-25"));
}

BOOST_AUTO_TEST_SUITE_END()
