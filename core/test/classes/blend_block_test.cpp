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

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

#include <boost/test/unit_test.hpp>

#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"
#include "bertini2/system/blocks/products_of_linears_block.hpp"

BOOST_AUTO_TEST_SUITE(blend_block_suite)

using namespace bertini;
using bertini::blocks::BlendBlock;
using bertini::blocks::ProductsOfLinearsBlock;
using bertini::DefaultPrecision;

using Blend = BlendBlock<ProductsOfLinearsBlock>;

static_assert(bertini::blocks::is_block_v<Blend>,
              "BlendBlock<ProductsOfLinearsBlock> must satisfy the block contract");

// A: f = x + 2y + 1     B: g = x + y - 5     (each one function, two variables)
// H = (1-t)*A + t*B
// at (x,y)=(1,1):  A = 4,  B = -3
//   t = 1/2:  H = 0.5*4 + 0.5*(-3) = 0.5
//   dH/dt = -A + B = -7
//   dH/dx = 0.5*[1,2] + 0.5*[1,1] = [1, 1.5]
static Blend MakeTestBlend()
{
	bertini::Mat<mpfr_complex> a(1, 3);
	a << mpfr_complex(1), mpfr_complex(2), mpfr_complex(1);
	ProductsOfLinearsBlock A(2, std::vector<bertini::Mat<mpfr_complex>>{a});

	bertini::Mat<mpfr_complex> b(1, 3);
	b << mpfr_complex(1), mpfr_complex(1), mpfr_complex(-5);
	ProductsOfLinearsBlock B(2, std::vector<bertini::Mat<mpfr_complex>>{b});

	auto t = node::Variable::Make("t");
	std::shared_ptr<node::Node> c0 = node::Integer::Make(1) - t; // 1 - t
	std::shared_ptr<node::Node> c1 = t;                          // t

	return Blend(t,
	             std::vector<std::shared_ptr<node::Node>>{c0, c1},
	             std::vector<ProductsOfLinearsBlock>{A, B});
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

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	bertini::Vec<dbl> result(1);
	blend.EvalInPlace<dbl>(result, x, dbl(0.5));

	BOOST_CHECK_CLOSE(result(0).real(), 0.5, 1e-11);
	BOOST_CHECK_SMALL(result(0).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	bertini::Mat<dbl> J(1, 2);
	blend.JacobianInPlace<dbl>(J, x, dbl(0.5));

	BOOST_CHECK_CLOSE(J(0, 0).real(), 1.0, 1e-11);
	BOOST_CHECK_CLOSE(J(0, 1).real(), 1.5, 1e-11);
}

BOOST_AUTO_TEST_CASE(time_derivative_double)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	bertini::Vec<dbl> dHdt(1);
	blend.TimeDerivInPlace<dbl>(dHdt, x, dbl(0.5));

	BOOST_CHECK_CLOSE(dHdt(0).real(), -7.0, 1e-11);
}

BOOST_AUTO_TEST_CASE(eval_and_time_derivative_mpfr)
{
	DefaultPrecision(30);
	auto blend = MakeTestBlend();
	blend.Precision(30);

	bertini::Vec<mpfr_complex> x(2); x << mpfr_complex(1), mpfr_complex(1);
	mpfr_complex t("0.5");

	bertini::Vec<mpfr_complex> result(1);
	blend.EvalInPlace<mpfr_complex>(result, x, t);
	BOOST_CHECK(abs(result(0) - mpfr_complex(1) / mpfr_complex(2)) < mpfr_float("1e-25"));

	bertini::Vec<mpfr_complex> dHdt(1);
	blend.TimeDerivInPlace<mpfr_complex>(dHdt, x, t);
	BOOST_CHECK(abs(dHdt(0) - mpfr_complex(-7)) < mpfr_float("1e-25"));
}

BOOST_AUTO_TEST_SUITE_END()
