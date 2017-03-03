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

using System = bertini::System;

using Variable = bertini::node::Variable;
using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;
using bertini::MakeVariable;
using mpq_rational = bertini::mpq_rational;
using mpfr_float = bertini::mpfr_float;
using mpz_int = bertini::mpz_int;
using dbl = bertini::dbl;
using mpfr = bertini::mpfr;
template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

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
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);

	sys.AddFunction(x*y);
	sys.AddFunction(pow(x,2)*pow(y,2));
	
	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	Vec<int> partition_1(2);
	partition_1 << 0, 1;

	Vec<int> partition_2(2);
	partition_2 << 1, 0;


	BOOST_CHECK(mhom_start_system.valid_partitions[0] == partition_1);
	BOOST_CHECK(mhom_start_system.valid_partitions[1] == partition_2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 2);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 4);

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
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");
	auto z = MakeVariable("z");


	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};
	VariableGroup v3{z};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);
	sys.AddHomVariableGroup(v3);

	sys.AddFunction(x*y);
	sys.AddFunction(pow(y,2)*pow(z,2));
	sys.AddFunction(pow(x,3)*pow(z,3));

	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	// for(int ii = 0; ii < mhom_start_system.valid_partitions.size(); ii++)
	// {
	// 	std::cout << "valid partition: " << std::endl <<  mhom_start_system.valid_partitions[ii] << std::endl; 
	// }

	Vec<int> partition_1(3);
	partition_1 << 0, 1,2;

	Vec<int> partition_2(3);
	partition_2 << 1, 2, 0;

	BOOST_CHECK(mhom_start_system.valid_partitions[0]== partition_1);
	BOOST_CHECK(mhom_start_system.valid_partitions[1]== partition_2);

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

BOOST_AUTO_TEST_SUITE_END()




