//This file is part of Bertini 2.0.
//
//node_serialization.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//node_serialization.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with node_serialization.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  node_serialization.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include <iostream>

#include <cstdlib>
#include <cmath>

#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>


using Variable = bertini::Variable;
using Node = bertini::Node;
using Float = bertini::Float;




extern double threshold_clearance_d;
extern unsigned FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS;
extern double threshold_clearance_mp;


BOOST_AUTO_TEST_SUITE(node_serialization)


BOOST_AUTO_TEST_CASE(serialize_variable)
{
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");


	{
		std::ofstream fout("serialization_test_node");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << x;
	}
	
	std::shared_ptr<Variable> x2;
	{
		std::ifstream fin("serialization_test_node");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> x2;
	}

	std::cout << *x2 << "\n";
}


BOOST_AUTO_TEST_CASE(serialize_float)
{

}

BOOST_AUTO_TEST_CASE(serialize_integer)
{

}


BOOST_AUTO_TEST_CASE(serialize_rational)
{

}

BOOST_AUTO_TEST_SUITE_END()




