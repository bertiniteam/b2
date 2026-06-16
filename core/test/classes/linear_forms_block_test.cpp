//This file is part of Bertini 2.
//
//linear_forms_block_test.cpp is free software: you can redistribute it and/or
//modify it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//linear_forms_block_test.cpp is distributed in the hope that it will be useful,
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
#include "bertini2/system/blocks/linear_forms_block.hpp"

BOOST_AUTO_TEST_SUITE(linear_forms_block_suite)

using namespace bertini;
using bertini::blocks::LinearFormsBlock;
using bertini::DefaultPrecision;

static_assert(bertini::blocks::is_block_v<LinearFormsBlock>,
              "LinearFormsBlock must satisfy the block contract");

// f0 = 2x + 3y + 1,  f1 = x - y + 4   (augmented rows: [coef_x, coef_y, const])
// at (x,y) = (1,1):  f0 = 6,  f1 = 4
//   Jacobian is constant: [[2, 3], [1, -1]]
static LinearFormsBlock MakeTestBlock()
{
	bertini::Mat<mpfr_complex> M(2, 3);
	M << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1),
	     mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);
	return LinearFormsBlock(2, std::move(M));
}

BOOST_AUTO_TEST_CASE(shape)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	BOOST_CHECK_EQUAL(block.NumFunctions(), 2u);
	BOOST_CHECK_EQUAL(block.NumVariables(), 2u);
	BOOST_CHECK(!block.DependsOnPathVariable());
	BOOST_CHECK(block.HasConstantJacobian());
}

BOOST_AUTO_TEST_CASE(eval_double)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	bertini::Vec<dbl> result(2);
	block.EvalInPlace<dbl>(result, x, dbl(0));

	BOOST_CHECK_CLOSE(result(0).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(result(1).real(), 4.0, 1e-11);
	BOOST_CHECK_SMALL(result(0).imag(), 1e-11);
	BOOST_CHECK_SMALL(result(1).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(eval_mpfr)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	block.Precision(30);

	bertini::Vec<mpfr_complex> x(2); x << mpfr_complex(1), mpfr_complex(1);
	bertini::Vec<mpfr_complex> result(2);
	block.EvalInPlace<mpfr_complex>(result, x, mpfr_complex(0));

	BOOST_CHECK(abs(result(0) - mpfr_complex(6)) < mpfr_float("1e-25"));
	BOOST_CHECK(abs(result(1) - mpfr_complex(4)) < mpfr_float("1e-25"));
}

BOOST_AUTO_TEST_CASE(jacobian_is_constant_coefficient_matrix)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();

	// evaluate the Jacobian at two different points; it must not change.
	bertini::Mat<dbl> J0(2, 2), J1(2, 2);
	bertini::Vec<dbl> p0(2); p0 << dbl(1), dbl(1);
	bertini::Vec<dbl> p1(2); p1 << dbl(-5), dbl(7);
	block.JacobianInPlace<dbl>(J0, p0, dbl(0));
	block.JacobianInPlace<dbl>(J1, p1, dbl(0));

	BOOST_CHECK_CLOSE(J0(0,0).real(), 2.0, 1e-11);
	BOOST_CHECK_CLOSE(J0(0,1).real(), 3.0, 1e-11);
	BOOST_CHECK_CLOSE(J0(1,0).real(), 1.0, 1e-11);
	BOOST_CHECK_CLOSE(J0(1,1).real(), -1.0, 1e-11);
	BOOST_CHECK_SMALL((J0 - J1).norm(), 1e-11);   // constant in x
}

BOOST_AUTO_TEST_CASE(time_derivative_is_zero)
{
	DefaultPrecision(30);
	auto block = MakeTestBlock();
	bertini::Vec<dbl> x(2); x << dbl(3), dbl(4);
	bertini::Vec<dbl> dt(2);
	block.TimeDerivInPlace<dbl>(dt, x, dbl(0));
	BOOST_CHECK_SMALL(dt.norm(), 1e-12);
}

// eval against an independent matrix-vector reference at a generic point, in mpfr.
BOOST_AUTO_TEST_CASE(eval_matches_matrix_vector_reference)
{
	DefaultPrecision(40);
	bertini::Mat<mpfr_complex> M(3, 4);
	M << mpfr_complex(2),  mpfr_complex(-1), mpfr_complex(0),  mpfr_complex(5),
	     mpfr_complex(1),  mpfr_complex(3),  mpfr_complex(-2), mpfr_complex(-1),
	     mpfr_complex(0),  mpfr_complex(4),  mpfr_complex(7),  mpfr_complex(2);
	LinearFormsBlock block(3, M);
	block.Precision(40);

	bertini::Vec<mpfr_complex> x(3); x << mpfr_complex(2), mpfr_complex(-3), mpfr_complex(1);
	bertini::Vec<mpfr_complex> result(3);
	block.EvalInPlace<mpfr_complex>(result, x, mpfr_complex(0));

	bertini::Vec<mpfr_complex> aug(4); aug << x(0), x(1), x(2), mpfr_complex(1);
	bertini::Vec<mpfr_complex> reference = M * aug;
	BOOST_CHECK((result - reference).norm() < mpfr_float("1e-30"));
}

BOOST_AUTO_TEST_SUITE_END()
