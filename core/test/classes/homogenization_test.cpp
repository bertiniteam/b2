//This file is part of Bertini 2.
//
//homogenization_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//homogenization_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with homogenization_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

//  homogenization_test.cpp



#include <boost/test/unit_test.hpp>



#include "bertini2/system/system.hpp"

#include "externs.hpp"


BOOST_AUTO_TEST_SUITE(homogenization)

using Variable = bertini::node::Variable;
using Complex = bertini::node::Complex;
using real_mp = bertini::real_mp;
using Var = std::shared_ptr<bertini::node::Variable>;

using Flt = std::shared_ptr<bertini::node::Complex>;
using VariableGroup = bertini::VariableGroup;
using complex_dbl = bertini::complex_dbl;
using mpfr = bertini::complex_mp;





BOOST_AUTO_TEST_CASE(no_homogenization_needed_x)
{
	Var x = Variable::Make("x");
	Var h = Variable::Make("h");

	std::shared_ptr<bertini::node::Node> f1 = x;  // widened to Node for functional Homogenized()
	
	BOOST_CHECK(f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK( f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));


	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_1)
{
	Var x = Variable::Make("x");
	Var t = Variable::Make("t");
	Var h = Variable::Make("h");

	auto f1 = x-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);

	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));


	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_1_minus_t_x_plus_t_1_minus_x)
{
	std::shared_ptr<bertini::node::Variable> x = Variable::Make("x");
	std::shared_ptr<bertini::node::Variable> t = Variable::Make("t");
	auto f1 = (1-t)*x + t*(1-x);

	BOOST_CHECK(!f1->IsHomogeneous());

	Var h = Variable::Make("h");
	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_t)
{
	Var x = Variable::Make("x");
	Var t = Variable::Make("t");
	Var h = Variable::Make("h");

	auto f1 = x-t;
	
	BOOST_CHECK(f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(no_homogenization_needed_x_minus_y_t)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");
	Var h = Variable::Make("h");

	auto f1 = x-y*t;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);
	vars.push_back(y);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}






BOOST_AUTO_TEST_CASE(homogenization_needed_sphere)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");
	Var h = Variable::Make("h");

	auto f1 = pow(x,2) + pow(y,2) + pow(z,2)-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(y));
	BOOST_CHECK(!f1->IsHomogeneous(z));
	BOOST_CHECK(!f1->IsHomogeneous(h));


	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadric)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");
	Var h = Variable::Make("h");

	auto f1 = x*y+x*z+y*z-1;
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(y));
	BOOST_CHECK(!f1->IsHomogeneous(z));
	BOOST_CHECK(!f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));


	BOOST_CHECK(f1->IsPolynomial());
}	





BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic)
{
	Var x = Variable::Make("x");
	Var h = Variable::Make("h");
	

	auto f1 = pow(x,2) + x + 1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant)
{
	Var x = Variable::Make("x");
	Var h = Variable::Make("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant_wrt_y)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var h = Variable::Make("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(y);

	f1 = f1->Homogenized(vars,h);  // functional: rebind to the homogenized copy
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));
	BOOST_CHECK( f1->IsHomogeneous(y));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}





BOOST_AUTO_TEST_CASE(nothomogeneous_sin_x)
{
	Var x = Variable::Make("x");

	auto f1 = sin(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_cos_x)
{
	Var x = Variable::Make("x");

	auto f1 = cos(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(nothomogeneous_tan_x)
{
	Var x = Variable::Make("x");

	auto f1 = tan(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_exp_x)
{
	Var x = Variable::Make("x");

	auto f1 = exp(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(nothomogeneous_sqrt_x)
{
	Var x = Variable::Make("x");

	auto f1 = sqrt(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(is_homogeneous_sin_0)
{
	Flt n = Complex::Make("1");

	auto f1 = sin(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_cos_1)
{
	Flt n = Complex::Make("1");

	auto f1 = cos(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_sin_1_plus_1)
{
	Flt n = Complex::Make("1");

	auto f1 = sin(n + n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(is_homogeneous_summands_homogeneous)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * real_mp("4.12331") * pow(x,1);
	
	auto f1 = a+b;
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsHomogeneous(x));
	BOOST_CHECK(f1->IsHomogeneous(y));


	BOOST_CHECK(f1->IsPolynomial());

}

BOOST_AUTO_TEST_CASE(not_homogeneous_summands_inhomogeneous)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * real_mp("4.12331");
	
	auto f1 = a+b;
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(y));

	BOOST_CHECK(f1->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(system_with_an_affine_products_of_linears_block_homogenizes_and_expands)
{
	// b2#376: the regeneration shape -- the products-of-linears start rows as a block beside a
	// polynomial slice -- must homogenize as a whole, expand to nodes with the block's values,
	// and patch.  The block's Homogenize was a no-op that reported success, so the expansion
	// threw "variable count mismatch" and the system could not be tracked projectively.
	bertini::DefaultPrecision(30);
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	VariableGroup vars{x, y};

	bertini::Mat<bertini::complex_mp> f0(2, 3);   // (x - 1)(x + 1)
	f0 << bertini::complex_mp(1), bertini::complex_mp(0), bertini::complex_mp(-1),
	      bertini::complex_mp(1), bertini::complex_mp(0), bertini::complex_mp(1);

	bertini::System sys;
	sys.AddVariableGroup(vars);
	sys.AddBlock(bertini::blocks::ProductsOfLinearsBlock(2, std::vector<bertini::Mat<bertini::complex_mp>>{f0}));
	sys.AddFunction(y - bertini::node::Rational::Make(1, 2, 0, 1));   // the static slice y = 1/2

	BOOST_CHECK(sys.IsPolynomial());
	BOOST_CHECK(!sys.IsHomogeneous());

	sys.Homogenize();
	BOOST_CHECK(sys.IsHomogeneous());
	BOOST_CHECK_EQUAL(sys.NumVariables(), 3u);

	// the node expansion agrees with the block evaluation at a generic projective point
	bertini::System expanded = sys.ExpandToFunctionTree();
	bertini::Vec<bertini::complex_mp> pt(3);
	pt << bertini::complex_mp("0.7", "0.2"), bertini::complex_mp("-0.4", "0.9"), bertini::complex_mp("1.3", "-0.5");
	auto a = sys.Eval(pt);
	auto b = expanded.Eval(pt);
	BOOST_REQUIRE_EQUAL(a.size(), b.size());
	for (Eigen::Index i = 0; i < a.size(); ++i)
		BOOST_CHECK(abs(a(i) - b(i)) < bertini::real_mp("1e-25"));

	BOOST_CHECK_NO_THROW(sys.AutoPatch());
	BOOST_CHECK(sys.IsPatched());
}



// ---- a variable-free subtree of ANY operator is a constant (b2#419) ----
//
// An algebraic constant written symbolically -- the golden ratio (1+sqrt(5))/2, or its
// Bertini 1 spelling (5^(1/2)+1)/2 -- is a polynomial of degree 0 and homogeneous, whatever
// operator builds it.  Every classifier must agree: Degree, IsPolynomial and IsHomogeneous
// are consulted by different callers, and when IsHomogeneous disagreed with Degree for a
// power with a non-integer exponent, Homogenize() found nothing to pad and AutoPatch()
// refused the result.  Table-driven over the operators, so a new one that breaks the
// agreement fails here by name.

namespace {

using Nd = std::shared_ptr<bertini::node::Node>;

struct VariableFreeCase
{
	std::string name;
	Nd node;
};

std::vector<VariableFreeCase> VariableFreeSubtrees()
{
	using namespace bertini::node;
	using bertini::mpq_rational;
	Nd one = Integer::Make(1), two = Integer::Make(2), three = Integer::Make(3), five = Integer::Make(5);
	return {
		{"sqrt(5)",        sqrt(five)},
		{"exp(1)",         exp(one)},
		{"log(2)",         log(two)},
		{"sin(1)",         sin(one)},
		{"cos(1)",         cos(one)},
		{"tan(1)",         tan(one)},
		{"5^(1/2)",        pow(five, mpq_rational(1, 2))},   // a power with a rational exponent: the #419 case
		{"5^(3/2)",        pow(five, mpq_rational(3, 2))},
		{"(1+sqrt(5))/2",  (sqrt(five) + one) / two},
		{"(5^(1/2)+1)/2",  (pow(five, mpq_rational(1, 2)) + one) / two},   // Bertini 1's spelling of phi
		{"sqrt(2)^3",      pow(sqrt(two), 3)},
		{"exp(sqrt(2))",   exp(sqrt(two))},
		{"-sqrt(3)",       -sqrt(three)},
		{"sqrt(pi)",       sqrt(Pi())},
	};
}

} // namespace


BOOST_AUTO_TEST_CASE(variable_free_subtree_of_every_operator_is_a_degree_zero_homogeneous_constant)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");
	VariableGroup vars{x, y, z};

	for (auto const& c : VariableFreeSubtrees())
	{
		BOOST_TEST_CONTEXT(c.name)
		{
			BOOST_CHECK_EQUAL(c.node->Degree(), 0);
			BOOST_CHECK_EQUAL(c.node->Degree(x), 0);
			BOOST_CHECK_EQUAL(c.node->Degree(vars), 0);
			BOOST_CHECK(c.node->IsPolynomial(vars));
			BOOST_CHECK(c.node->IsHomogeneous(x));
			BOOST_CHECK(c.node->IsHomogeneous(vars));

			// used as a coefficient, it must leave the polynomial's degree and homogeneity alone
			Nd f = c.node * pow(x, 2) - pow(y, 2) + pow(z, 2);
			BOOST_CHECK_EQUAL(f->Degree(vars), 2);
			BOOST_CHECK(f->IsPolynomial(vars));
			BOOST_CHECK(f->IsHomogeneous(vars));
		}
	}
}


BOOST_AUTO_TEST_CASE(system_with_a_symbolic_constant_coefficient_homogenizes_and_patches)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");
	VariableGroup vars{x, y, z};

	for (auto const& c : VariableFreeSubtrees())
	{
		BOOST_TEST_CONTEXT(c.name)
		{
			bertini::System sys;
			sys.AddVariableGroup(vars);
			sys.AddFunction(c.node * pow(x, 2) - pow(y, 2) + pow(z, 2));   // already homogeneous
			sys.AddFunction(c.node * pow(x, 2) - y);                        // needs a homogenizing variable on y

			BOOST_CHECK(sys.IsPolynomial());
			BOOST_CHECK(!sys.IsHomogeneous());

			sys.Homogenize();
			BOOST_CHECK(sys.IsHomogeneous());
			BOOST_CHECK_EQUAL(sys.NumVariables(), 4u);
			BOOST_CHECK_NO_THROW(sys.AutoPatch());
			BOOST_CHECK(sys.IsPatched());
		}
	}
}


BOOST_AUTO_TEST_CASE(barth_sextic_with_the_golden_ratio_as_a_node_homogenizes_and_patches)
{
	// the reproduction from b2#419, with phi in Bertini 1's spelling (5^(1/2)+1)/2 -- a power
	// with a rational exponent, which is how bertini_real's test/surface/barth6/input writes it
	using namespace bertini::node;
	using bertini::mpq_rational;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");
	VariableGroup vars{x, y, z};

	Nd one = Integer::Make(1), two = Integer::Make(2), four = Integer::Make(4), five = Integer::Make(5);
	Nd phi = (pow(five, mpq_rational(1, 2)) + one) / two;
	Nd phi2 = pow(phi, 2);
	Nd f = four * (phi2 * pow(x, 2) - pow(y, 2)) * (phi2 * pow(y, 2) - pow(z, 2)) * (phi2 * pow(z, 2) - pow(x, 2))
	       - (one + two * phi) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2) - one, 2);

	BOOST_CHECK_EQUAL(f->Degree(vars), 6);
	BOOST_CHECK(f->IsPolynomial(vars));

	bertini::System sys;
	sys.AddVariableGroup(vars);
	sys.AddFunction(f);
	sys.Homogenize();
	BOOST_CHECK(sys.IsHomogeneous());
	BOOST_CHECK_NO_THROW(sys.AutoPatch());
	BOOST_CHECK(sys.IsPatched());
}


BOOST_AUTO_TEST_CASE(variable_dependent_non_polynomial_subtrees_stay_non_polynomial)
{
	// the other side of the same line: once a variable is inside, these are not polynomials,
	// not homogeneous, and their degree is negative
	using namespace bertini::node;
	using bertini::mpq_rational;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	VariableGroup vars{x, y};

	std::vector<VariableFreeCase> cases = {
		{"sqrt(x)",   sqrt(x)},
		{"exp(x)",    exp(x)},
		{"log(x)",    log(x)},
		{"sin(x)",    sin(x)},
		{"x^(1/2)",   pow(x, mpq_rational(1, 2))},
		{"sqrt(x*y)", sqrt(x * y)},
		{"2^x",       pow(Integer::Make(2), x)},
	};
	for (auto const& c : cases)
	{
		BOOST_TEST_CONTEXT(c.name)
		{
			BOOST_CHECK_LT(c.node->Degree(vars), 0);
			BOOST_CHECK(!c.node->IsPolynomial(vars));
			BOOST_CHECK(!c.node->IsHomogeneous(vars));
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()





