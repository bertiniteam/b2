//This file is part of Bertini 2.
//
//blend_block_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//blend_block_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"

BOOST_AUTO_TEST_SUITE(blend_block_suite)

using namespace bertini;
using bertini::blocks::BlendBlock;
using bertini::DefaultPrecision;

using Blend = BlendBlock<System>;

static_assert(bertini::blocks::is_block_v<Blend>,
              "BlendBlock<System> must satisfy the block contract");

// A: f = x + 2y     B: f = x*y     (variable group {x,y}, one function each, no patch)
// H = (1-t)*A + t*B
// at (x,y)=(1,1):  A = 3,  B = 1
//   t = 1/2:  H = 0.5*3 + 0.5*1 = 2
//   dH/dt = -A + B = -2
//   dH/dx = 0.5*[1,2] + 0.5*[y,x]=0.5*[1,1] = [1, 1.5]
static Blend MakeTestBlend()
{
	auto x = node::Variable::Make("x");
	auto y = node::Variable::Make("y");

	auto A = std::make_shared<System>();
	A->AddVariableGroup(VariableGroup{x, y});
	A->AddFunction(x + 2 * y);

	auto B = std::make_shared<System>();
	B->AddVariableGroup(VariableGroup{x, y});
	B->AddFunction(x * y);

	auto t = node::Variable::Make("t");
	std::vector<std::shared_ptr<node::Node>> coeffs{node::Integer::Make(1) - t, t};
	std::vector<std::shared_ptr<const System>> ops{A, B};

	return Blend(t, coeffs, ops);
}

BOOST_AUTO_TEST_CASE(shape)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();
	BOOST_CHECK_EQUAL(blend.NumFunctions(), 1u);
	BOOST_CHECK(blend.DependsOnPathVariable());
}

BOOST_AUTO_TEST_CASE(eval_double)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Vec<complex_dbl> result(1);
	blend.EvalInPlace<complex_dbl>(result, x, complex_dbl(0.5));

	BOOST_CHECK_CLOSE(result(0).real(), 2.0, 1e-11);
	BOOST_CHECK_SMALL(result(0).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Mat<complex_dbl> J(1, 2);
	blend.JacobianInPlace<complex_dbl>(J, x, complex_dbl(0.5));

	BOOST_CHECK_CLOSE(J(0, 0).real(), 1.0, 1e-11);
	BOOST_CHECK_CLOSE(J(0, 1).real(), 1.5, 1e-11);
}

BOOST_AUTO_TEST_CASE(time_derivative_double)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Vec<complex_dbl> dHdt(1);
	blend.TimeDerivInPlace<complex_dbl>(dHdt, x, complex_dbl(0.5));

	BOOST_CHECK_CLOSE(dHdt(0).real(), -2.0, 1e-11);
}

BOOST_AUTO_TEST_CASE(eval_and_time_derivative_mpfr)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();
	blend.Precision(30);

	bertini::Vec<complex_mp> x(2); x << complex_mp(1), complex_mp(1);
	complex_mp t("0.5");

	bertini::Vec<complex_mp> result(1);
	blend.EvalInPlace<complex_mp>(result, x, t);
	BOOST_CHECK(abs(result(0) - complex_mp(2)) < real_mp("1e-25"));

	bertini::Vec<complex_mp> dHdt(1);
	blend.TimeDerivInPlace<complex_mp>(dHdt, x, t);
	BOOST_CHECK(abs(dHdt(0) - complex_mp(-2)) < real_mp("1e-25"));
}

// The coefficients are evaluated through a cached coefficient System (compiled once,
// reused), not by node-level tree evaluation.  Repeated calls --- and a call at a changed
// value --- must keep giving the right answer, exercising the cache across invocations.
BOOST_AUTO_TEST_CASE(repeated_evaluation_reuses_coefficients_consistently)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	bertini::Vec<complex_dbl> result(1);

	// First call builds the cached coefficient system; subsequent calls reuse it.
	blend.EvalInPlace<complex_dbl>(result, x, complex_dbl(0.5));
	BOOST_CHECK_CLOSE(result(0).real(), 2.0, 1e-11);

	blend.EvalInPlace<complex_dbl>(result, x, complex_dbl(0.5));
	BOOST_CHECK_CLOSE(result(0).real(), 2.0, 1e-11);

	// A different path value: H = (1-t)*3 + t*1 = 3 - 2t; at t=0 -> 3.
	blend.EvalInPlace<complex_dbl>(result, x, complex_dbl(0.0));
	BOOST_CHECK_CLOSE(result(0).real(), 3.0, 1e-11);

	// And at t=1 -> 1.
	blend.EvalInPlace<complex_dbl>(result, x, complex_dbl(1.0));
	BOOST_CHECK_CLOSE(result(0).real(), 1.0, 1e-11);
}

BOOST_AUTO_TEST_SUITE_END()
