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


BOOST_AUTO_TEST_CASE(m_hom_start_system_construction_of_degree_matrix_num_start_pts_and_partitions)
{

	/* Test case to check if we are creating a degree matrix correctly. 
	   This is not checking how homogenization or patching effects our MHomogeneous start system.
	*/
	DefaultPrecision(30);

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");

	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);

	sys.AddFunction(x*y - 1);
	sys.AddFunction(pow(x,2) - 1);
	
	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	// std::cout << "deg mat is " << std::endl;
	// std::cout << mhom_start_system.degree_matrix_ << std::endl;

	// std::cout << "1st partition is " << std::endl;
	// std::cout << mhom_start_system.valid_partitions[0] << std::endl;



	Vec<int> partition_check(2);
	partition_check << 1, 0;

	BOOST_CHECK(mhom_start_system.valid_partitions.front() == partition_check);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 0);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 2);

}


BOOST_AUTO_TEST_CASE(m_hom_start_system_bigger_system_contruction_test)
{

	/* Test case to check if we are creating a degree matrix correctly. 
	   This is not checking how homogenization or patching effects our MHomogeneous start system.
	*/
	DefaultPrecision(30);

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var z = std::make_shared<Variable>("z");

	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};
	VariableGroup v3{z};	

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);
	sys.AddHomVariableGroup(v3);

	sys.AddFunction(x*y);
	sys.AddFunction(pow(x,2)*pow(y,2) + 2*z);
	sys.AddFunction(x*z);

	
	auto mhom_start_system = bertini::start_system::MHomogeneous(sys);

	Vec<int> partition_1(3);
	partition_1 << 0, 1, 2;

	Vec<int> partition_2(3);
	partition_2 << 1, 0, 2;

	Vec<int> partition_3(3);
	partition_3 << 1, 2, 0;



	BOOST_CHECK(mhom_start_system.valid_partitions[0] == partition_1);
	BOOST_CHECK(mhom_start_system.valid_partitions[1] == partition_2);
	BOOST_CHECK(mhom_start_system.valid_partitions[2] == partition_3);


	BOOST_CHECK(mhom_start_system.degree_matrix_(0,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,1) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(0,2) == 0);

	BOOST_CHECK(mhom_start_system.degree_matrix_(1,0) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,1) == 2);
	BOOST_CHECK(mhom_start_system.degree_matrix_(1,2) == 1);

	BOOST_CHECK(mhom_start_system.degree_matrix_(2,0) == 1);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,1) == 0);
	BOOST_CHECK(mhom_start_system.degree_matrix_(2,2) == 1);

	BOOST_CHECK(mhom_start_system.NumStartPoints() == 5);

}


BOOST_AUTO_TEST_SUITE_END()




