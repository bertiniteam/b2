//This file is part of Bertini 2.
//
//system_blocks_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_blocks_test.cpp is distributed in the hope that it will be useful,
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
#include "bertini2/system/blocks/products_of_linears_block.hpp"

BOOST_AUTO_TEST_SUITE(system_blocks_suite)

using namespace bertini;
using Var = std::shared_ptr<node::Variable>;
using bertini::blocks::ProductsOfLinearsBlock;
using bertini::DefaultPrecision;

// A block-composed System with one affine variable group {x,y} (unhomogenized, no
// patch) whose single function is the product of linears (2x+3y+1)(x-y+4).
// At (x,y)=(1,1): value = 6*4 = 24; Jacobian = [14, 6].
static System MakeBlockSystem()
{
	System sys;
	Var x = node::Variable::Make("x"), y = node::Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});

	bertini::Mat<mpfr_complex> f(2, 3);
	f << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1),
	     mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);
	sys.AddBlock(ProductsOfLinearsBlock(2, std::vector<bertini::Mat<mpfr_complex>>{f}));
	return sys;
}

BOOST_AUTO_TEST_CASE(counts)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();
	BOOST_CHECK(sys.HasBlocks());
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 1u);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2u);
	BOOST_CHECK_EQUAL(sys.NumTotalFunctions(), 1u); // no patch
}

BOOST_AUTO_TEST_CASE(eval_double)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	auto v = sys.Eval(x);

	BOOST_CHECK_EQUAL(v.size(), 1);
	BOOST_CHECK_CLOSE(v(0).real(), 24.0, 1e-11);
	BOOST_CHECK_SMALL(v(0).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();

	bertini::Vec<dbl> x(2); x << dbl(1), dbl(1);
	sys.Eval(x);                       // sets the current variable values
	auto J = sys.Jacobian<dbl>();      // uses the current variable values

	BOOST_CHECK_EQUAL(J.rows(), 1);
	BOOST_CHECK_EQUAL(J.cols(), 2);
	BOOST_CHECK_CLOSE(J(0, 0).real(), 14.0, 1e-11);
	BOOST_CHECK_CLOSE(J(0, 1).real(),  6.0, 1e-11);
}

BOOST_AUTO_TEST_CASE(eval_mpfr)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();
	sys.precision(30);                 // propagate precision to variables and blocks

	bertini::Vec<mpfr_complex> x(2); x << mpfr_complex(1), mpfr_complex(1);
	auto v = sys.Eval(x);

	BOOST_CHECK(abs(v(0) - mpfr_complex(24)) < mpfr_float("1e-25"));
}

BOOST_AUTO_TEST_SUITE_END()
