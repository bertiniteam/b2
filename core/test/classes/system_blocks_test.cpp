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

#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdexcept>

#include "bertini2/system/system.hpp"
#include "bertini2/system/slice.hpp"
#include "bertini2/system/blocks/products_of_linears_block.hpp"
#include "bertini2/system/blocks/linear_forms_block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"

BOOST_AUTO_TEST_SUITE(system_blocks_suite)

using namespace bertini;
using Var = std::shared_ptr<node::Variable>;
using bertini::blocks::ProductsOfLinearsBlock;
using bertini::blocks::LinearFormsBlock;
using bertini::blocks::BlendBlock;
using bertini::DefaultPrecision;

// A block-composed System with one affine variable group {x,y} (unhomogenized, no
// patch) whose single function is the product of linears (2x+3y+1)(x-y+4).
// At (x,y)=(1,1): value = 6*4 = 24; Jacobian = [14, 6].
static System MakeBlockSystem()
{
	System sys;
	Var x = node::Variable::Make("x"), y = node::Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});

	bertini::Mat<complex_mp> f(2, 3);
	f << complex_mp(2), complex_mp(3),  complex_mp(1),
	     complex_mp(1), complex_mp(-1), complex_mp(4);
	sys.AddBlock(ProductsOfLinearsBlock(2, std::vector<bertini::Mat<complex_mp>>{f}));
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

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	auto v = sys.Eval(x);

	BOOST_CHECK_EQUAL(v.size(), 1);
	BOOST_CHECK_CLOSE(v(0).real(), 24.0, 1e-11);
	BOOST_CHECK_SMALL(v(0).imag(), 1e-11);
}

BOOST_AUTO_TEST_CASE(jacobian_double)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();

	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	sys.Eval(x);                       // sets the current variable values
	auto J = sys.Jacobian<complex_dbl>();      // uses the current variable values

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

	bertini::Vec<complex_mp> x(2); x << complex_mp(1), complex_mp(1);
	auto v = sys.Eval(x);

	BOOST_CHECK(abs(v(0) - complex_mp(24)) < real_mp("1e-25"));
}

// issue #263: Function(i) on a system with no PolynomialBlock used to dereference a null
// PolyBlockPtr() and segfault.  It now expands the structured block to a node on demand.
BOOST_AUTO_TEST_CASE(function_accessor_on_structured_block)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();                 // a single products-of-linears function, no PolynomialBlock

	BOOST_CHECK_EQUAL(sys.GetNaturalFunctions().size(), 1u);  // not empty, in sync with NumNaturalFunctions
	auto f = sys.Function(0);
	BOOST_CHECK(f != nullptr);                                // a real function-tree node
	BOOST_CHECK_THROW(sys.Function(1), std::out_of_range);    // past the end raises, no segfault
}

// issue #263: Slices() recovers the linear-form slices embedded in a system (inverse of AsSystem).
BOOST_AUTO_TEST_CASE(slices_round_trip_linear_forms)
{
	DefaultPrecision(30);
	Var x = node::Variable::Make("x"), y = node::Variable::Make("y");
	VariableGroup vars{x, y};

	bertini::Mat<complex_mp> M(2, 3);             // 2x+3y+1 ; x-y+4  (augmented, last col = constant)
	M << complex_mp(2), complex_mp(3),  complex_mp(1),
	     complex_mp(1), complex_mp(-1), complex_mp(4);

	System sys;
	sys.AddVariableGroup(vars);
	sys.AddBlock(LinearFormsBlock(2, M));

	auto slices = sys.Slices();
	BOOST_REQUIRE_EQUAL(slices.size(), 1u);
	auto const& C = slices[0].Coefficients();
	BOOST_CHECK_EQUAL(C.rows(), 2);
	BOOST_CHECK_EQUAL(C.cols(), 3);
	BOOST_CHECK(abs(C(0, 0) - complex_mp(2)) < real_mp("1e-25"));
	BOOST_CHECK(abs(C(1, 2) - complex_mp(4)) < real_mp("1e-25"));

	System empty;                                  // a system with no linear-forms block has no slices
	empty.AddVariableGroup(vars);
	empty.AddFunction(x * x - y);
	BOOST_CHECK(empty.Slices().empty());
}

// round-trip a System through a boost text archive, returning the deserialized copy.
static System RoundTrip(System const& sys)
{
	std::stringstream ss;
	{
		boost::archive::text_oarchive oa(ss);
		oa << sys;
	}
	System out;
	{
		boost::archive::text_iarchive ia(ss);
		ia >> out;
	}
	return out;
}

BOOST_AUTO_TEST_CASE(serialize_products_of_linears_block)
{
	DefaultPrecision(30);
	auto sys = MakeBlockSystem();
	bertini::Vec<complex_dbl> x(2); x << complex_dbl(1), complex_dbl(1);
	auto before = sys.Eval(x);

	auto sys2 = RoundTrip(sys);
	BOOST_REQUIRE(sys2.HasBlocks());
	auto after = sys2.Eval(x);

	BOOST_CHECK_EQUAL(after.size(), before.size());
	BOOST_CHECK_CLOSE(after(0).real(), before(0).real(), 1e-11);  // 24
}

BOOST_AUTO_TEST_CASE(serialize_linear_forms_block)
{
	DefaultPrecision(30);
	System sys;
	Var x = node::Variable::Make("x"), y = node::Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	// f0 = 2x + 3y + 1, f1 = x - y + 4   (augmented rows)
	bertini::Mat<complex_mp> M(2, 3);
	M << complex_mp(2), complex_mp(3),  complex_mp(1),
	     complex_mp(1), complex_mp(-1), complex_mp(4);
	sys.AddBlock(LinearFormsBlock(2, M));

	bertini::Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto before = sys.Eval(p);          // [6, 4]

	auto sys2 = RoundTrip(sys);
	BOOST_REQUIRE(sys2.HasBlocks());
	auto after = sys2.Eval(p);

	BOOST_REQUIRE_EQUAL(after.size(), 2);
	BOOST_CHECK_CLOSE(after(0).real(), before(0).real(), 1e-11);
	BOOST_CHECK_CLOSE(after(1).real(), before(1).real(), 1e-11);
}

BOOST_AUTO_TEST_CASE(serialize_blend_block)
{
	DefaultPrecision(30);
	// H(x,t) = (1-t)(x-2) + t(x-5) = x - 2 - 3t, blending two single-function systems.
	auto t = node::Variable::Make("t");
	auto x = node::Variable::Make("x");

	auto target = std::make_shared<System>();
	target->AddVariableGroup(VariableGroup{x});
	target->AddFunction(x - node::Integer::Make(2));

	auto start = std::make_shared<System>();
	start->AddVariableGroup(VariableGroup{x});
	start->AddFunction(x - node::Integer::Make(5));

	System H;
	H.AddVariableGroup(VariableGroup{x});
	H.AddPathVariable(t);
	std::vector<std::shared_ptr<node::Node>> coeffs{ node::Integer::Make(1) - t, t };
	std::vector<std::shared_ptr<const System>> operands{ target, start };
	H.AddBlock(BlendBlock<System>(t, std::move(coeffs), std::move(operands)));

	bertini::Vec<complex_dbl> p(1); p << complex_dbl(1);
	auto before = H.Eval(p, complex_dbl(0));        // x - 2 - 0 = -1

	auto H2 = RoundTrip(H);
	BOOST_REQUIRE(H2.HasBlocks());
	auto after = H2.Eval(p, complex_dbl(0));

	BOOST_REQUIRE_EQUAL(after.size(), 1);
	BOOST_CHECK_CLOSE(after(0).real(), before(0).real(), 1e-9);  // -1
}

BOOST_AUTO_TEST_SUITE_END()
