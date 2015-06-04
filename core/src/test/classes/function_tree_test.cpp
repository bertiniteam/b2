//This file is part of Bertini 2.0.
//
//function_tree_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_tree_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_tree_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  function_tree_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

#include <cmath>
#include <iostream>
#include <vector>

#include "bertini.hpp"
#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>



#include <boost/test/unit_test.hpp>


using Variable = bertini::Variable;
using Node = bertini::Node;


BOOST_AUTO_TEST_SUITE(function_tree_class)



BOOST_AUTO_TEST_CASE(manual_construction_x_squared){
	
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
	
	std::shared_ptr<Node> N = x;

	N *= N;
    
	x->set_current_value<dbl>(3.5);
	

	BOOST_CHECK_EQUAL(N->Eval<dbl>() , 12.25 );
	
}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_x){
	
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
	
	std::shared_ptr<Node> N = pow(x, 1.0/2);
	
    x->set_current_value<dbl>(std::complex<double>(3.0,0.0));
	
    BOOST_CHECK(abs(N->Eval<dbl>().real() - 1.73205080756887729352744634151) < 1e-15);  // root 3 to 30 digits
	BOOST_CHECK(abs(N->Eval<dbl>().imag() - 0.0) < 1e-15);
}


BOOST_AUTO_TEST_CASE(manual_construction_x_plus_y){
    
    std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
    
    std::shared_ptr<Node> N = pow(x, 1.0/2);
    
    x->set_current_value<dbl>(std::complex<double>(3.0,0.0));
    
    BOOST_CHECK(abs(N->Eval<dbl>().real() - 1.73205080756887729352744634151) < 1e-15);  // root 3 to 30 digits
    BOOST_CHECK(abs(N->Eval<dbl>().imag() - 0.0) < 1e-15);
}




BOOST_AUTO_TEST_SUITE_END()
