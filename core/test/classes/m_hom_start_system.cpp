//This file is part of Bertini 2.
//
//m_hom_start_system.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//m_hom_start_system.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with m_hom_start_system.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Tim Hodges, Colorado State University

#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"
#include <iostream>
#include <iomanip>

using System = bertini::System;

using Variable = bertini::node::Variable;
using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;

using mpq_rational = bertini::mpq_rational;
using mpfr_float = bertini::mpfr_float;
using mpz_int = bertini::mpz_int;
using dbl = bertini::dbl;
using mpfr = bertini::mpfr_complex;
template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;

#include "externs.hpp"

using namespace bertini::start_system;

using bertini::DefaultPrecision;

BOOST_AUTO_TEST_SUITE(m_hom_system_class)

BOOST_AUTO_TEST_CASE(m_hom_system_preliminary_construction_small_example)
{

	/* 
		Test case to check if we are creating a degree matrix correctly. 
		This is not checking how homogenization or patching effects our MHomogeneous start system.

	   	f = x*y;
	   	g = x^2*y^2;  f, g are homogeneous w.r.t x and y.

		degree matrix: [1 1]
					   [2 2]
		valid partitions <0,1>, <1,0> these are the entries in each row to grab.
		
		num start points is 1*2 + 1*2 = 4.			   

	*/
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	// second projective coordinate of each homogeneous group: a hom group of size k is
	// P^{k-1}, so a size-1 group is the degenerate P^0.  Use size-2 groups (P^1, dimension
	// 1) so each group absorbs one function -- the degree matrix, partitions, and start-
	// point count are exactly as for the original (capacity = dimension = size - 1 = 1).
	auto x1 = Variable::Make("x1");
	auto y1 = Variable::Make("y1");

	System sys;

	VariableGroup v1{x, x1};
	VariableGroup v2{y, y1};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);

	sys.AddFunction(x*y);
	sys.AddFunction(pow(x,2)*pow(y,2));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	Vec<int> partition_1(2);
	partition_1 << 0, 1;

	Vec<int> partition_2(2);
	partition_2 << 1, 0;


	BOOST_CHECK(mhom_start_system.valid_partitions_[0] == partition_1);
	BOOST_CHECK(mhom_start_system.valid_partitions_[1] == partition_2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 2);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 4);

}


// Each generated start point must be an actual root of the (homogenized + patched)
// start system: evaluating the start system there is ~0.  This validates both the
// start-point linear solve and its homogenization onto the patch (the part the solve
// flow needs).
BOOST_AUTO_TEST_CASE(start_points_are_roots_of_the_start_system)
{
	DefaultPrecision(30);
	bertini::SetGlobalSeed(1u);

	System sys;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x*y - 1);
	sys.AddFunction(x + y);
	sys.Homogenize();
	sys.AutoPatch();

	auto mhom = MHomogeneous(sys);

	// Partition capacity respects each group's declared dimension even after homogenization
	// (VariableGroupSizes() would count the homogenizing variable and double it): the
	// m-homogeneous Bezout number here is 2, not 4.
	BOOST_CHECK_EQUAL(mhom.NumStartPoints(), 2ull);

	const auto n = mhom.NumStartPoints();
	for (unsigned long long i = 0; i < n; ++i)
	{
		auto sp = mhom.StartPoint<dbl>(i);
		BOOST_CHECK_EQUAL(static_cast<size_t>(sp.size()), mhom.NumVariables());
		auto v = mhom.Eval(sp);
		for (Eigen::Index j = 0; j < v.size(); ++j)
			BOOST_CHECK(std::abs(v(j)) < 1e-10);   // each start point is a root, on the patch
	}

	// The blend homotopy (built as FormHomotopy does) is correct: zero at the start points
	// at t=1, and its Jacobian and dH/dt agree with finite differences.
	auto t = Variable::Make("t");
	auto gamma = bertini::node::Rational::Make(bertini::node::Rational::Rand());
	System H = sys;                 // target's variable structure + patch
	H.ClearFunctions();             // the blend supplies the rows; don't also eval target's own functions
	H.AddPathVariable(t);
	std::vector<std::shared_ptr<bertini::node::Node>> coeffs{ 1 - t, gamma * t };
	std::vector<std::shared_ptr<const System>> operands{
		std::make_shared<System>(sys),
		std::make_shared<System>(mhom) };
	H.AddBlock(bertini::blocks::BlendBlock<System>(t, coeffs, operands));

	for (unsigned long long i = 0; i < n; ++i)
	{
		auto Hval = H.Eval(mhom.StartPoint<dbl>(i), dbl(1));
		for (Eigen::Index j = 0; j < Hval.size(); ++j)
			BOOST_CHECK(std::abs(Hval(j)) < 1e-9);
	}

	bertini::Vec<dbl> xq(H.NumVariables());
	for (Eigen::Index k = 0; k < xq.size(); ++k)
		xq(k) = dbl(0.37 * (k + 1) + 0.11, -0.19 * k + 0.07);
	const dbl tv(0.42, -0.13);
	const dbl hstep(1e-6, 0);

	auto J = H.Jacobian(xq, tv);
	auto f0 = H.Eval(xq, tv);
	for (Eigen::Index c = 0; c < xq.size(); ++c)
	{
		auto xp = xq; xp(c) += hstep;
		auto fp = H.Eval(xp, tv);
		for (Eigen::Index r = 0; r < f0.size(); ++r)
			BOOST_CHECK(std::abs((fp(r) - f0(r)) / hstep - J(r, c)) < 1e-5);
	}

	auto dHdt = H.TimeDerivative(xq, tv);
	auto fpt = H.Eval(xq, tv + hstep);
	for (Eigen::Index r = 0; r < f0.size(); ++r)
		BOOST_CHECK(std::abs((fpt(r) - f0(r)) / hstep - dHdt(r)) < 1e-5);
}




BOOST_AUTO_TEST_CASE(m_hom_system_preliminary_construction_larger_example)
{

	/* Test case to check if we are creating a degree matrix correctly. 
	   This is not checking how homogenization or patching effects our MHomogeneous start system.

	  	f = x*y;
	   	g = y^2*z^2;  
	   	h = x^3*z^3;  f, g, and h are homogeneous w.r.t x, y, and z.

		degree matrix: [1 1 0]
					   [0 2 2]
					   [3 0 3]
		valid partitions <0,1,2>, <1,2,0> these are the entries in each row to grab.
		
		num start points is 1*2*3 + 1*2*3 = 12.	
	*/
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	// second projective coordinate of each homogeneous group (see the small example above):
	// a hom group of size k is P^{k-1}, so use size-2 groups (P^1, dimension 1).  Degree
	// matrix, partitions, and start-point count are exactly as for size-1 groups under the
	// correct convention (capacity = dimension = size - 1 = 1).
	auto x1 = Variable::Make("x1");
	auto y1 = Variable::Make("y1");
	auto z1 = Variable::Make("z1");


	System sys;

	VariableGroup v1{x, x1};
	VariableGroup v2{y, y1};
	VariableGroup v3{z, z1};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);
	sys.AddHomVariableGroup(v3);

	sys.AddFunction(x*y);
	sys.AddFunction(pow(y,2)*pow(z,2));
	sys.AddFunction(pow(x,3)*pow(z,3));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);


	Vec<int> partition_1(3);
	partition_1 << 0, 1,2;

	Vec<int> partition_2(3);
	partition_2 << 1, 2, 0;

	BOOST_CHECK(mhom_start_system.valid_partitions_[0]== partition_1);
	BOOST_CHECK(mhom_start_system.valid_partitions_[1]== partition_2);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,2) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,2) == 2);

	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 3);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,1) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,2) == 3);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 12);

}

BOOST_AUTO_TEST_CASE(Two_var_groups_4_vars_4_fctns_example)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1,x2}, {x3,x4}
	// 		f1 = x1^2 + x4
	//  	f2 = x1*x2*x3
	//  	f3 = x4^2 
	//  	f4 = x3^2
	//  	Degree matrix 
	//  	[2 1]
	//		[2 1]
	//		[0 2]
	//		[0 2]
	//		Number of paths: 16
 	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");
	auto x3 = Variable::Make("x3");
	auto x4 = Variable::Make("x4");


	System sys;

	VariableGroup v1{x1,x2};
	VariableGroup v2{x3,x4};

	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);

	sys.AddFunction(pow(x1,2) + x4);
	sys.AddFunction(x1*x2*x3);
	sys.AddFunction(pow(x4,2));
	sys.AddFunction(pow(x3,2));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,0) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,1) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,1) == 2);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 16);

}

BOOST_AUTO_TEST_CASE(four_var_groups_4_vars_4_fctns_example)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1}, {x2}, {x3}, {x4}
	// 		f1 = x1^2 + x4
	//  	f2 = x1*x2*x3
	//  	f3 = x4^2 
	//  	f4 = x3^2
	//  	Degree matrix 
	//  	[2 0 0 1]
	//		[1 1 1 0]
	//		[0 0 0 2]
	//		[0 0 2 0]
	//		Number of paths: 8
	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");
	auto x3 = Variable::Make("x3");
	auto x4 = Variable::Make("x4");


	System sys;

	VariableGroup v1{x1};
	VariableGroup v2{x2};
	VariableGroup v3{x3};
	VariableGroup v4{x4};

	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);
	sys.AddVariableGroup(v3);
	sys.AddVariableGroup(v4);

	sys.AddFunction(pow(x1,2) + x4);
	sys.AddFunction(x1*x2*x3);
	sys.AddFunction(pow(x4,2));
	sys.AddFunction(pow(x3,2));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,0) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,1) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,1) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,2) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,2) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,2) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,2) == 2);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,3) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,3) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,3) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,3) == 0);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 8);

}

BOOST_AUTO_TEST_CASE(one_var_group_4_vars_4_fctns_example)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1,x2,x3,x4}
	// 		f1 = x1^2 + x4
	//  	f2 = x1*x2*x3
	//  	f3 = x4^2 
	//  	f4 = x3^2
	//  	Degree matrix 
	//  	[2]
	//		[3]
	//		[2]
	//		[2]
	//		Number of paths: 24
	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");
	auto x3 = Variable::Make("x3");
	auto x4 = Variable::Make("x4");


	System sys;

	VariableGroup v1{x1,x2,x3,x4};

	sys.AddVariableGroup(v1);

	sys.AddFunction(pow(x1,2) + x4);
	sys.AddFunction(x1*x2*x3);
	sys.AddFunction(pow(x4,2));
	sys.AddFunction(pow(x3,2));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 3);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(3,0) == 2);


	BOOST_CHECK(mhom_start_system.NumStartPoints() == 24);

}


BOOST_AUTO_TEST_CASE(zero_column_in_degree_matrix)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1,x2},{x3}
	// 		f1 = x1 + x2
	//  	f2 = x1
	//  	f3 = x2
	//  	Degree matrix 
	//  	[1 0]
	//		[1 0]
	//		[1 0]
	//		Should throw for having a zero column.
	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");
	auto x3 = Variable::Make("x3");



	System sys;

	VariableGroup v1{x1,x2};
	VariableGroup v2{x3};


	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);


	sys.AddFunction(x1 + x2);
	sys.AddFunction(x1);
	sys.AddFunction(x2);

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}

BOOST_AUTO_TEST_CASE(diagonal_degree_matrix)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1},{x2},{x3}
	// 		f1 = x1^2
	//  	f2 = x2^2
	//  	f3 = x3^2
	//  	Degree matrix 
	//  	[2 0 0]
	//		[0 2 0]
	//		[0 0 2]
	//		Number of paths: 8
	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");
	auto x3 = Variable::Make("x3");

	System sys;

	VariableGroup v1{x1};
	VariableGroup v2{x2};
	VariableGroup v3{x3};


	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);
	sys.AddVariableGroup(v3);

	sys.AddFunction(pow(x1,2));
	sys.AddFunction(pow(x2,2));
	sys.AddFunction(pow(x3,2));


	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,1) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(0,2) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,2) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,2) == 2);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 8);

}


BOOST_AUTO_TEST_CASE(zero_row_in_degree_matrix)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x1},{x2}
	// 		f1 = x1^2 + x2^2
	//  	f2 = 1
	//  	Degree matrix 
	//  	[2 2]
	//		[0 0]
	//		Should throw for having a zero row.
	auto x1 = Variable::Make("x1");
	auto x2 = Variable::Make("x2");

	System sys;

	VariableGroup v1{x1};
	VariableGroup v2{x2};


	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);

	sys.AddFunction(pow(x1,2) + pow(x2,2));
	sys.AddFunction(bertini::node::Integer::Make(1));

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}


BOOST_AUTO_TEST_CASE(non_square_target_system_for_mhom)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x}
	// 		f1 = x^2
	//  	f2 = x^2 - 1
	// Will throw for not being square
	auto x = Variable::Make("x");


	System sys;

	VariableGroup v1{x};

	sys.AddVariableGroup(v1);


	sys.AddFunction(pow(x,2));
	sys.AddFunction(pow(x,2) - 1);

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}


BOOST_AUTO_TEST_CASE(empty_variable_group)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x,y},{}
	// 		f1 = x^2
	//  	f2 = y^2 - 1
	// Should throw for having an empty variable group. 
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	System sys;

	VariableGroup v1{x,y};
	VariableGroup v2{};

	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);


	sys.AddFunction(pow(x,2));
	sys.AddFunction(pow(y,2) - 1);

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}


BOOST_AUTO_TEST_CASE(ungrouped_variable_for_target_system_in_mhom_construction)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x}, y ungrouped
	// 		f1 = x^2
	//  	f2 = y^2 - 1
	// Should throw for having an ungrouped variable.
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	System sys;

	VariableGroup v1{x};

	sys.AddVariableGroup(v1);


	sys.AddFunction(pow(x,2));
	sys.AddFunction(pow(y,2) - 1);

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}

BOOST_AUTO_TEST_CASE(variable_in_many_variable_groups_in_mhom_construction)
{
	//	Basic problem from Dan Bates 
	//  	variable groups {x,y},{y} 
	// 		f1 = x^2
	//  	f2 = y^2 - 1
	// Should throw for having the variable y in many variable groups. 
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	System sys;

	VariableGroup v1{x,y};
	VariableGroup v2{y};

	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);


	sys.AddFunction(pow(x,2));
	sys.AddFunction(pow(y,2) - 1);

	BOOST_CHECK_THROW(auto mhom_start_system = bertini::start_system::MHomogeneous(sys), std::runtime_error);

}

BOOST_AUTO_TEST_SUITE_END()




