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

#include "system.hpp"
#include "system_parsing.hpp"

using Variable = bertini::Variable;
using Node = bertini::Node;
using Float = bertini::Float;


using System = bertini::System;

extern double threshold_clearance_d;
extern unsigned FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS;
extern double threshold_clearance_mp;

template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

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

	BOOST_CHECK(x->name()==x2->name());
}


BOOST_AUTO_TEST_CASE(serialize_float)
{
	std::shared_ptr<Float> two_point_oh_four = std::make_shared<Float>(2.04);

	{
		std::ofstream fout("serialization_test_node");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << two_point_oh_four;
	}
	
	std::shared_ptr<Float> two_point_oh_four2;
	{
		std::ifstream fin("serialization_test_node");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> two_point_oh_four2;
	}

	BOOST_CHECK(two_point_oh_four->Eval<dbl>()==two_point_oh_four2->Eval<dbl>());
}

BOOST_AUTO_TEST_CASE(serialize_complicated_expression)
{
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");

	auto f = exp(sqrt(pow(pow(x*x+ (-x) -sin(x)+cos(x)+tan(x),x),3)))/x;

	{
		std::ofstream fout("serialization_test_node");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << x;
		oa << f;
	}
	
	std::shared_ptr<Node> f2;
	std::shared_ptr<Variable> x2;
	{
		std::ifstream fin("serialization_test_node");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> x2;
		ia >> f2;
	}

	BOOST_CHECK(x->name()==x2->name());

	x->set_current_value(dbl(1.2,0.9));
	x2->set_current_value(dbl(1.2,0.9));

	BOOST_CHECK(abs(f->Eval<dbl>() - f2->Eval<dbl>()) < threshold_clearance_d);

}




BOOST_AUTO_TEST_CASE(system_serialize)
{
	std::string str = "function f1, f2; variable_group x1, x2; y = x1*x2; f1 = y*y; f2 = x1*y; ";

	bertini::System sys;
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bertini::SystemParser<std::string::const_iterator> S;
	phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);

	

	
	{
		std::ofstream fout("serialization_test_node");
		
		boost::archive::text_oarchive oa(fout);
		
		// write class instance to archive
		oa << sys;
	}
	
	bertini::System sys2;
	{
		std::ifstream fin("serialization_test_node");
		
		boost::archive::text_iarchive ia(fin);
		// read class state from archive
		ia >> sys2;
	}

	

	Vec<dbl> values(2);

	values(0) = dbl(2.0);
	values(1) = dbl(3.0);

	Vec<dbl> v = sys2.Eval(values);

	BOOST_CHECK_EQUAL(v.size(),2);

	BOOST_CHECK_EQUAL(v(0), 36.0);
	BOOST_CHECK_EQUAL(v(1), 12.0);
	

}


BOOST_AUTO_TEST_SUITE_END()




